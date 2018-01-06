#include "Physics.h"
#include <iostream>
#include "FieldOptimization.h"
#include <Eigen/Dense>

using namespace Eigen;

Matrix3d crossMatrix(const Vector3d &n)
{
    Matrix3d ret;
    ret << 0, -n[2], n[1],
        n[2], 0, -n[0],
        -n[1], n[0], 0;
    return ret;
}

void precomputeFundamentalForms(
    const Eigen::MatrixXd &V,          // vertices, |V| x 3
    const Eigen::Vector3i &faceVerts,  // the triangle's vertices, in orientation order
    const Eigen::Vector3i &oppVerts,   // the triangle's opposite vertices, in "one left out" order (q0 is the vertex *opposite* p0 across the common edge)
    Eigen::Matrix2d &a,        // as above 
    Eigen::Matrix2d &b,
    Eigen::MatrixXd &da,
    Eigen::MatrixXd &db
)
{
    da.resize(3, 9);
    db.resize(3, 18);
    da.setZero();
    db.setZero();

    Vector3d p[3];
    Vector3d q[3];
    for (int i = 0; i < 3; i++)
    {
        p[i] = V.row(faceVerts[i]).transpose();
        q[i] = V.row(oppVerts[i]).transpose();
    }

    a << (p[1] - p[0]).dot(p[1] - p[0]), (p[1] - p[0]).dot(p[2] - p[0]),
        (p[1] - p[0]).dot(p[2] - p[0]), (p[2] - p[0]).dot(p[2] - p[0]);

    da.block<1, 3>(0, 0) -= 2.0 * (p[1] - p[0]).transpose();
    da.block<1, 3>(0, 3) += 2.0 * (p[1] - p[0]).transpose();

    da.block<1, 3>(1, 0) -= (p[2] - p[0]).transpose();
    da.block<1, 3>(1, 0) -= (p[1] - p[0]).transpose();
    da.block<1, 3>(1, 3) += (p[2] - p[0]).transpose();
    da.block<1, 3>(1, 6) += (p[1] - p[0]).transpose();

    da.block<1, 3>(2, 0) -= 2.0*(p[2] - p[0]).transpose();
    da.block<1, 3>(2, 6) += 2.0*(p[2] - p[0]).transpose();

    Vector3d n[3];
    Vector3d tilden[3];
    Matrix3d I;
    I.setIdentity();

    for (int i = 0; i < 3; i++)
    {
        int ip2 = (i + 2) % 3;
        int ip1 = (i + 1) % 3;
        n[i] = (p[ip1] - p[i]).cross(p[ip2] - p[i]);
        tilden[i] = (p[ip2] - q[i]).cross(p[ip1] - q[i]);
    }

    // mid-edge normals
    Vector3d m[3];
    for (int i = 0; i < 3; i++)
    {
        m[i] = n[i] / n[i].norm() + tilden[i] / tilden[i].norm();
    }

    Vector3d mhat[3];
    for (int i = 0; i < 3; i++)
    {
        mhat[i] = m[i] / m[i].norm();
    }

    b << 0.5 * (mhat[0] - mhat[1]).dot(p[1] - p[0]), 0.5 * (mhat[0] - mhat[1]).dot(p[2] - p[0]),
        0.5 * (mhat[0] - mhat[2]).dot(p[1] - p[0]), 0.5*(mhat[0] - mhat[2]).dot(p[2] - p[0]);

    // b is symmetric, though this is not obvious from the definition
    if (fabs(b(0, 1) - b(1, 0)) > 1e-8)
    {
        std::cerr << "Bug in b code" << std::endl;
        exit(-1);
    }

    db.block<1, 3>(0, 0) -= 0.5 * (mhat[0] - mhat[1]).transpose();
    db.block<1, 3>(0, 3) += 0.5 * (mhat[0] - mhat[1]).transpose();
    db.block<1, 3>(1, 0) -= 0.5 * (mhat[0] - mhat[1]).transpose();
    db.block<1, 3>(1, 6) += 0.5 * (mhat[0] - mhat[1]).transpose();
    db.block<1, 3>(2, 0) -= 0.5 * (mhat[0] - mhat[2]).transpose();
    db.block<1, 3>(2, 6) += 0.5 * (mhat[0] - mhat[2]).transpose();

    Matrix3d dmiterm1[3];
    Matrix3d dmiterm2[3];
    Matrix3d dmiterm3[3];
    Matrix3d dmiterm4[3];
    for (int i = 0; i < 3; i++)
    {
        int ip1 = (i + 1) % 3;
        int ip2 = (i + 2) % 3;
        dmiterm1[i] = 1.0 / m[i].norm() * (I - m[i] * m[i].transpose() / m[i].squaredNorm()) * 1.0 / n[i].norm() * (I - n[i] * n[i].transpose() / n[i].squaredNorm()) * crossMatrix(p[ip1] - p[i]);
        dmiterm2[i] = -1.0 / m[i].norm() * (I - m[i] * m[i].transpose() / m[i].squaredNorm()) * 1.0 / n[i].norm() * (I - n[i] * n[i].transpose() / n[i].squaredNorm()) * crossMatrix(p[ip2] - p[i]);
        dmiterm3[i] = 1.0 / m[i].norm() * (I - m[i] * m[i].transpose() / m[i].squaredNorm()) * 1.0 / tilden[i].norm() * (I - tilden[i] * tilden[i].transpose() / tilden[i].squaredNorm()) * crossMatrix(p[ip2] - q[i]);
        dmiterm4[i] = -1.0 / m[i].norm() * (I - m[i] * m[i].transpose() / m[i].squaredNorm()) * 1.0 / tilden[i].norm() * (I - tilden[i] * tilden[i].transpose() / tilden[i].squaredNorm()) * crossMatrix(p[ip1] - q[i]);
    }

    // db_11
    db.block<1, 3>(0, 0) -= 0.5 * (p[1] - p[0]).transpose() * dmiterm1[0];
    db.block<1, 3>(0, 6) += 0.5 * (p[1] - p[0]).transpose() * dmiterm1[0];
    db.block<1, 3>(0, 0) -= 0.5 * (p[1] - p[0]).transpose() * dmiterm2[0];
    db.block<1, 3>(0, 3) += 0.5 * (p[1] - p[0]).transpose() * dmiterm2[0];
    db.block<1, 3>(0, 9) -= 0.5 * (p[1] - p[0]).transpose() * dmiterm3[0];
    db.block<1, 3>(0, 3) += 0.5 * (p[1] - p[0]).transpose() * dmiterm3[0];
    db.block<1, 3>(0, 9) -= 0.5 * (p[1] - p[0]).transpose() * dmiterm4[0];
    db.block<1, 3>(0, 6) += 0.5 * (p[1] - p[0]).transpose() * dmiterm4[0];

    db.block<1, 3>(0, 3) += 0.5 * (p[1] - p[0]).transpose() * dmiterm1[1];
    db.block<1, 3>(0, 0) -= 0.5 * (p[1] - p[0]).transpose() * dmiterm1[1];
    db.block<1, 3>(0, 3) += 0.5 * (p[1] - p[0]).transpose() * dmiterm2[1];
    db.block<1, 3>(0, 6) -= 0.5 * (p[1] - p[0]).transpose() * dmiterm2[1];
    db.block<1, 3>(0, 12) += 0.5 * (p[1] - p[0]).transpose() * dmiterm3[1];
    db.block<1, 3>(0, 6) -= 0.5 * (p[1] - p[0]).transpose() * dmiterm3[1];
    db.block<1, 3>(0, 12) += 0.5 * (p[1] - p[0]).transpose() * dmiterm4[1];
    db.block<1, 3>(0, 0) -= 0.5 * (p[1] - p[0]).transpose() * dmiterm4[1];

    // db_12
    db.block<1, 3>(1, 0) -= 0.5 * (p[2] - p[0]).transpose() * dmiterm1[0];
    db.block<1, 3>(1, 6) += 0.5 * (p[2] - p[0]).transpose() * dmiterm1[0];
    db.block<1, 3>(1, 0) -= 0.5 * (p[2] - p[0]).transpose() * dmiterm2[0];
    db.block<1, 3>(1, 3) += 0.5 * (p[2] - p[0]).transpose() * dmiterm2[0];
    db.block<1, 3>(1, 9) -= 0.5 * (p[2] - p[0]).transpose() * dmiterm3[0];
    db.block<1, 3>(1, 3) += 0.5 * (p[2] - p[0]).transpose() * dmiterm3[0];
    db.block<1, 3>(1, 9) -= 0.5 * (p[2] - p[0]).transpose() * dmiterm4[0];
    db.block<1, 3>(1, 6) += 0.5 * (p[2] - p[0]).transpose() * dmiterm4[0];

    db.block<1, 3>(1, 3) += 0.5 * (p[2] - p[0]).transpose() * dmiterm1[1];
    db.block<1, 3>(1, 0) -= 0.5 * (p[2] - p[0]).transpose() * dmiterm1[1];
    db.block<1, 3>(1, 3) += 0.5 * (p[2] - p[0]).transpose() * dmiterm2[1];
    db.block<1, 3>(1, 6) -= 0.5 * (p[2] - p[0]).transpose() * dmiterm2[1];
    db.block<1, 3>(1, 12) += 0.5 * (p[2] - p[0]).transpose() * dmiterm3[1];
    db.block<1, 3>(1, 6) -= 0.5 * (p[2] - p[0]).transpose() * dmiterm3[1];
    db.block<1, 3>(1, 12) += 0.5 * (p[2] - p[0]).transpose() * dmiterm4[1];
    db.block<1, 3>(1, 0) -= 0.5 * (p[2] - p[0]).transpose() * dmiterm4[1];

    // db_22
    db.block<1, 3>(2, 0) -= 0.5 * (p[2] - p[0]).transpose() * dmiterm1[0];
    db.block<1, 3>(2, 6) += 0.5 * (p[2] - p[0]).transpose() * dmiterm1[0];
    db.block<1, 3>(2, 0) -= 0.5 * (p[2] - p[0]).transpose() * dmiterm2[0];
    db.block<1, 3>(2, 3) += 0.5 * (p[2] - p[0]).transpose() * dmiterm2[0];
    db.block<1, 3>(2, 9) -= 0.5 * (p[2] - p[0]).transpose() * dmiterm3[0];
    db.block<1, 3>(2, 3) += 0.5 * (p[2] - p[0]).transpose() * dmiterm3[0];
    db.block<1, 3>(2, 9) -= 0.5 * (p[2] - p[0]).transpose() * dmiterm4[0];
    db.block<1, 3>(2, 6) += 0.5 * (p[2] - p[0]).transpose() * dmiterm4[0];

    db.block<1, 3>(2, 6) += 0.5 * (p[2] - p[0]).transpose() * dmiterm1[2];
    db.block<1, 3>(2, 3) -= 0.5 * (p[2] - p[0]).transpose() * dmiterm1[2];
    db.block<1, 3>(2, 6) += 0.5 * (p[2] - p[0]).transpose() * dmiterm2[2];
    db.block<1, 3>(2, 0) -= 0.5 * (p[2] - p[0]).transpose() * dmiterm2[2];
    db.block<1, 3>(2, 15) += 0.5 * (p[2] - p[0]).transpose() * dmiterm3[2];
    db.block<1, 3>(2, 0) -= 0.5 * (p[2] - p[0]).transpose() * dmiterm3[2];
    db.block<1, 3>(2, 15) += 0.5 * (p[2] - p[0]).transpose() * dmiterm4[2];
    db.block<1, 3>(2, 3) -= 0.5 * (p[2] - p[0]).transpose() * dmiterm4[2];
}

void physicsDataFromMesh(const MeshData &mesh, PhysicsData &phys)
{
    int nfaces = mesh.F.rows();
    
    phys.as.resize(nfaces);
    phys.bs.resize(nfaces);
    phys.das.resize(nfaces);
    phys.dbs.resize(nfaces);

    for (int i = 0; i < nfaces; i++)
    {
        precomputeFundamentalForms(mesh.V, mesh.F.row(i), mesh.faceWings.row(i), phys.as[i], phys.bs[i], phys.das[i], phys.dbs[i]);
    }
}

double physicalEnergy(const MeshData &mesh, const PhysicsData &fundamentalForms, const Eigen::VectorXd &widthFractions, const Eigen::VectorXd &vectorField, double kb, double kt)
{
    double result = 0;

    int nfaces = mesh.F.rows();
    for (int i = 0; i < nfaces; i++)
    {
        Matrix2d ainv = fundamentalForms.as[i].inverse();
        Vector2d v = vectorField.segment<2>(2 * i);
        Vector2d vperp = (ainv * v) * (v.transpose() * fundamentalForms.as[i] * v) / (v.transpose() * ainv * v);
        double af = sqrt(fundamentalForms.as[i].determinant()) * 0.5;
        Vector3i faceverts = mesh.F.row(i);
        Vector3i flapverts = mesh.faceWings.row(i);
        double w = widthFractions[i];

        double vTbv = v.transpose() * fundamentalForms.bs[i] * v;
        double JvTbv = vperp.transpose() * fundamentalForms.bs[i] * v;
        double vTav = v.transpose() * fundamentalForms.as[i] * v;

        result += af * w * 0.5 * kb * (vTbv * vTbv) / (vTav * vTav);
        result += af * w * 0.5 * kt * (JvTbv * JvTbv) / (vTav * vTav);
    }
    return result;
}

void physicalForces(const MeshData &mesh, const PhysicsData &fundamentalForms, const Eigen::VectorXd &widthFractions, const Eigen::VectorXd &vectorField, double kb, double kt,
    Eigen::MatrixXd &forceField)
{
    int nfaces = mesh.F.rows();
    int nverts = mesh.V.rows();
    forceField.resize(nverts, 3);
    forceField.setZero();
    // energy is
    // \sum_faces A_f w_f kb/2 ( v_f^T b v_f )^2 / (v_f^T a v_f)^2 + \sum_faces A_f w_f kt/2 (Jv_f^T b v_f)^2 / (v_f^T a v_f)^2

    for (int i = 0; i < nfaces; i++)
    {
        Matrix2d ainv = fundamentalForms.as[i].inverse();
        Vector2d v = vectorField.segment<2>(2 * i);
        Vector2d vperp = (ainv * v) * (v.transpose() * fundamentalForms.as[i] * v) / (v.transpose() * ainv * v);
        double af = sqrt(fundamentalForms.as[i].determinant()) * 0.5;
        Vector3i faceverts = mesh.F.row(i);
        Vector3i flapverts = mesh.faceWings.row(i);
        double w = widthFractions[i];

        double vTbv = v.transpose() * fundamentalForms.bs[i] * v;
        double JvTbv = vperp.transpose() * fundamentalForms.bs[i] * v;
        double vTav = v.transpose() * fundamentalForms.as[i] * v;

        // derivative of b terms
        double bendcoeff1 = af * w * kb * vTbv / vTav / vTav;
        for (int j = 0; j < 3; j++)
        {
            forceField.row(faceverts[j]) += bendcoeff1 * v[0] * v[0] * fundamentalForms.dbs[i].block<1, 3>(0, 3 * j);
            forceField.row(faceverts[j]) += bendcoeff1 * v[0] * v[1] * fundamentalForms.dbs[i].block<1, 3>(1, 3 * j);
            forceField.row(faceverts[j]) += bendcoeff1 * v[1] * v[0] * fundamentalForms.dbs[i].block<1, 3>(1, 3 * j);
            forceField.row(faceverts[j]) += bendcoeff1 * v[1] * v[1] * fundamentalForms.dbs[i].block<1, 3>(2, 3 * j);

            forceField.row(flapverts[j]) += bendcoeff1 * v[0] * v[0] * fundamentalForms.dbs[i].block<1, 3>(0, 9 + 3 * j);
            forceField.row(flapverts[j]) += bendcoeff1 * v[0] * v[1] * fundamentalForms.dbs[i].block<1, 3>(1, 9 + 3 * j);
            forceField.row(flapverts[j]) += bendcoeff1 * v[1] * v[0] * fundamentalForms.dbs[i].block<1, 3>(1, 9 + 3 * j);
            forceField.row(flapverts[j]) += bendcoeff1 * v[1] * v[1] * fundamentalForms.dbs[i].block<1, 3>(2, 9 + 3 * j);
        }        
        double twistcoeff1 = af * w * kt * JvTbv / vTav / vTav;
        for (int j = 0; j < 3; j++)
        {
            forceField.row(faceverts[j]) += twistcoeff1 * vperp[0] * v[0] * fundamentalForms.dbs[i].block<1, 3>(0, 3 * j);
            forceField.row(faceverts[j]) += twistcoeff1 * vperp[0] * v[1] * fundamentalForms.dbs[i].block<1, 3>(1, 3 * j);
            forceField.row(faceverts[j]) += twistcoeff1 * vperp[1] * v[0] * fundamentalForms.dbs[i].block<1, 3>(1, 3 * j);
            forceField.row(faceverts[j]) += twistcoeff1 * vperp[1] * v[1] * fundamentalForms.dbs[i].block<1, 3>(2, 3 * j);

            forceField.row(flapverts[j]) += twistcoeff1 * vperp[0] * v[0] * fundamentalForms.dbs[i].block<1, 3>(0, 9 + 3 * j);
            forceField.row(flapverts[j]) += twistcoeff1 * vperp[0] * v[1] * fundamentalForms.dbs[i].block<1, 3>(1, 9 + 3 * j);
            forceField.row(flapverts[j]) += twistcoeff1 * vperp[1] * v[0] * fundamentalForms.dbs[i].block<1, 3>(1, 9 + 3 * j);
            forceField.row(flapverts[j]) += twistcoeff1 * vperp[1] * v[1] * fundamentalForms.dbs[i].block<1, 3>(2, 9 + 3 * j);
        }        

        // derivative of a terms
        double bendcoeff2 = -af * w * kb * vTbv * vTbv / vTav / vTav / vTav;
        for (int j = 0; j < 3; j++)
        {
            forceField.row(faceverts[j]) += bendcoeff2 * v[0] * v[0] * fundamentalForms.das[i].block<1, 3>(0, 3 * j);
            forceField.row(faceverts[j]) += bendcoeff2 * v[0] * v[1] * fundamentalForms.das[i].block<1, 3>(1, 3 * j);
            forceField.row(faceverts[j]) += bendcoeff2 * v[1] * v[0] * fundamentalForms.das[i].block<1, 3>(1, 3 * j);
            forceField.row(faceverts[j]) += bendcoeff2 * v[1] * v[1] * fundamentalForms.das[i].block<1, 3>(2, 3 * j);
        }

        double twistcoeff2 = -af * w * kb * JvTbv * JvTbv / vTav / vTav / vTav;
        for (int j = 0; j < 3; j++)
        {
            forceField.row(faceverts[j]) += twistcoeff2 * v[0] * v[0] * fundamentalForms.das[i].block<1, 3>(0, 3 * j);
            forceField.row(faceverts[j]) += twistcoeff2 * v[0] * v[1] * fundamentalForms.das[i].block<1, 3>(1, 3 * j);
            forceField.row(faceverts[j]) += twistcoeff2 * v[1] * v[0] * fundamentalForms.das[i].block<1, 3>(1, 3 * j);
            forceField.row(faceverts[j]) += twistcoeff2 * v[1] * v[1] * fundamentalForms.das[i].block<1, 3>(2, 3 * j);
        }
    }

    // force is *negative* derivative of energy
    forceField *= -1.0;
}

void testForces(const MeshData &mesh, const PhysicsData &fundamentalForms, const Eigen::VectorXd &widthFractions, const Eigen::VectorXd &vectorField, double kb, double kt, int vertex, int coord)
{
    double origenergy = physicalEnergy(mesh, fundamentalForms, widthFractions, vectorField, kb, kt);
    Eigen::MatrixXd forceField;
    physicalForces(mesh, fundamentalForms, widthFractions, vectorField, kb, kt, forceField);

    MeshData newmesh = mesh;
    newmesh.V(vertex, coord) += 1e-6;
    PhysicsData newdata;
    physicsDataFromMesh(newmesh, newdata);
    double newenergy = physicalEnergy(newmesh, newdata, widthFractions, vectorField, kb, kt);
    double findiff = (newenergy - origenergy) / 1e-6;
    std::cout << findiff << " vs " << forceField(vertex, coord) << std::endl;
}