#include "Physics.h"
#include <iostream>
#include "FieldOptimization.h"

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
    MatrixXi faceWings(nfaces, 3);
    for (int i = 0; i < nfaces; i++) {
        for (int j = 0; j < 3; j++) {
            int result = -1;
            int p1 = mesh.F(i, (j + 1) % 3);
            int p2 = mesh.F(i, (j + 2) % 3);
            for (int k = 0; k < nfaces; k++) {
                for (int l = 0; l < 3; l++) {
                    if (mesh.F(k, (l + 1) % 3) == p2 && mesh.F(k, (l + 2) % 3) == p1) {
                        result = mesh.F(k, l);
                    }
                }
            }
            faceWings(i, j) = result;
        }
    }

    phys.as.resize(nfaces);
    phys.bs.resize(nfaces);
    phys.das.resize(nfaces);
    phys.dbs.resize(nfaces);

    for (int i = 0; i < nfaces; i++)
    {
        precomputeFundamentalForms(mesh.V, mesh.F.row(i), faceWings.row(i), phys.as[i], phys.bs[i], phys.das[i], phys.dbs[i]);
    }
}