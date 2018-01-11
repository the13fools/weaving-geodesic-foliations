#include "Weave.h"
#include <igl/read_triangle_mesh.h>
#include <map>
#include <Eigen/Dense>
#include "Colors.h"

Weave::Weave(const std::string &objname, int m)
{
    Eigen::MatrixXd Vtmp;
    if (!igl::read_triangle_mesh(objname, Vtmp, F))
    {
        std::cerr << "Couldn't load mesh " << objname << std::endl;
        exit(-1);
    }
    if (Vtmp.cols() < 3)
    {
        std::cerr << "Mesh must 3D" << std::endl;
        exit(-1);
    }
    V.resize(Vtmp.rows(), 3);
    //wtf
    for (int i = 0; i < 3; i++)
    {
        V.col(i) = Vtmp.col(i);
    }

    buildConnectivityStructures();
    buildGeometricStructures();

    // initialize vector fields
    nFields_ = m;
    vectorFields.resize(5*F.rows()*m);
    vectorFields.setZero();
    vectorFields.segment(0, 2 * F.rows()*m).setRandom();
    normalizeFields();

    // initialize permutation matrices
    int nedges = nEdges();
    Ps.resize(nedges);
    for (int i = 0; i < nedges; i++)
    {
        Ps[i].resize(m, m);
        Ps[i].setIdentity();
    }
}

Weave::~Weave()
{    
}

void Weave::buildConnectivityStructures()
{
    std::map<std::pair<int, int>, Eigen::Vector2i > edgemap;
    int nfaces = nFaces();
    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int nextj = (j + 1) % 3;
            int v1 = F(i, j);
            int v2 = F(i, nextj);
            int idx = 0;
            if (v1 > v2)
            {
                idx = 1;
                std::swap(v1, v2);
            }
            edgemap[std::pair<int, int>(v1, v2)][idx] = i;
        }
    }

    E.resize(edgemap.size(), 2);
    faceEdges.resize(nfaces, 3);
    edgeVerts.resize(edgemap.size(), 2);
    int idx = 0;
    for (std::map<std::pair<int, int>, Eigen::Vector2i >::iterator it = edgemap.begin(); it != edgemap.end(); ++it)
    {
        if (it->second.size() != 2)
        {
            std::cerr << "Code only supports watertight manifold meshes without boundary" << std::endl;
            exit(-1);
        }
        E(idx, 0) = it->second[0];
        E(idx, 1) = it->second[1];
        edgeVerts(idx, 0) = it->first.first;
        edgeVerts(idx, 1) = it->first.second;
        idx++;
    }
    faceEdges.resize(nfaces, 3);
    faceNeighbors.resize(nfaces, 3);
    faceWings.resize(nfaces, 3);
    int nedges = nEdges();
    for (int edge = 0; edge < nedges; edge++)
    {
        Eigen::Vector3i face1 = F.row(E(edge, 0));
        Eigen::Vector3i face2 = F.row(E(edge, 1));
        int idx1 = -1;
        for (int i = 0; i < 3; i++)
        {
            bool ok = true;
            for (int j = 0; j < 3; j++)
            {
                if (face2[j] == face1[i])
                    ok = false;
            }
            if (ok)
                idx1 = i;
        }
        int idx2 = -1;
        for (int i = 0; i < 3; i++)
        {
            bool ok = true;
            for (int j = 0; j < 3; j++)
            {
                if (face1[j] == face2[i])
                    ok = false;
            }
            if (ok)
                idx2 = i;
        }
        faceEdges(E(edge,0), idx1) = edge;
        faceEdges(E(edge,1), idx2) = edge;
        faceNeighbors(E(edge, 0), idx1) = E(edge, 1);
        faceNeighbors(E(edge, 1), idx2) = E(edge, 0);
        faceWings(E(edge, 0), idx1) = face2[idx2];
        faceWings(E(edge, 1), idx2) = face1[idx1];
    }
}

Eigen::Vector3d Weave::faceNormal(int face)
{
    Eigen::Vector3d e1 = (V.row(F(face, 1)) - V.row(F(face, 0)));
    Eigen::Vector3d e2 = (V.row(F(face, 2)) - V.row(F(face, 0)));
    Eigen::Vector3d result = e1.cross(e2);
    result /= result.norm();
    return result;
}

void Weave::buildGeometricStructures()
{
    // compute barycentric matrices and Js
    int nfaces = nFaces();
    Bs.resize(nfaces);
    Js.resize(2 * nfaces, 2);

    averageEdgeLength = 0;

    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector3d v0 = V.row(F(i, 0));
        Eigen::Vector3d v1 = V.row(F(i, 1));
        Eigen::Vector3d v2 = V.row(F(i, 2));
        Bs[i].col(0) = v1 - v0;
        Bs[i].col(1) = v2 - v0;

        averageEdgeLength += (v2 - v1).norm();
        averageEdgeLength += (v2 - v0).norm();
        averageEdgeLength += (v1 - v0).norm();

        Eigen::Vector3d n = (v1 - v0).cross(v2 - v0);
        n /= n.norm();

        Eigen::Matrix2d BTB = Bs[i].transpose() * Bs[i];
        Eigen::Matrix<double, 3, 2> ncrossB;
        ncrossB.col(0) = n.cross(v1 - v0);
        ncrossB.col(1) = n.cross(v2 - v0);
        Js.block<2, 2>(2 * i, 0) = BTB.inverse() * Bs[i].transpose() * ncrossB;
    }

    averageEdgeLength /= (3.0 * nfaces);

    // compute cDiffs and transition matrices
    int nedges = nEdges();
    cDiffs.resize(2 * nedges, 2);
    Ts.resize(2 * nedges, 4);
    for (int edgeidx = 0; edgeidx < nedges; edgeidx++)
    {
        //collect neighboring face indices
        int face1 = E(edgeidx, 0);
        int face2 = E(edgeidx, 1);

        int v1 = edgeVerts(edgeidx, 0);
        int v2 = edgeVerts(edgeidx, 1);

        Eigen::Vector3d n1 = faceNormal(face1);
        Eigen::Vector3d n2 = faceNormal(face2);
        // collect: (1) the midpoint of the common edge, (2) unit vector in direction of common edge,
        // (3) the face normals, (4) the centroid of neighboring faces

        Eigen::Vector3d midpt = 0.5 * (V.row(v1) + V.row(v2));
        Eigen::Vector3d commone = V.row(v2) - V.row(v1);
        commone /= commone.norm();
        Eigen::Vector3d centroids[2];
        centroids[0].setZero();
        centroids[1].setZero();

        for (int i = 0; i < 3; i++)
        {
            centroids[0] += V.row(F(face1, i));
            centroids[1] += V.row(F(face2, i));
        }

        centroids[0] /= 3.0;
        centroids[1] /= 3.0;

        //rotate each centroid into the plane of the opposite triangle and compute ci minus c

        Eigen::Vector3d t1 = n1.cross(commone);
        Eigen::Vector3d t2 = n2.cross(commone);
        Eigen::Vector3d diff2 = centroids[1] - midpt;
        double alpha = commone.dot(diff2);
        double beta = t2.dot(diff2);

        Eigen::Matrix2d BTB1 = Bs[face1].transpose() * Bs[face1];
        Eigen::Matrix2d BTB2 = Bs[face2].transpose() * Bs[face2];

        cDiffs.row(2 * edgeidx) = BTB1.inverse() * Bs[face1].transpose() * (midpt + alpha * commone + beta * t1 - centroids[0]);
        Eigen::Vector3d diff1 = centroids[0] - midpt;
        alpha = commone.dot(diff1);
        beta = t1.dot(diff1);
        cDiffs.row(2 * edgeidx + 1) = BTB2.inverse() * Bs[face2].transpose() * (midpt + alpha*commone + beta * t2 - centroids[1]);

        Eigen::Vector3d e1 = V.row(F(face1, 1)) - V.row(F(face1, 0));
        Eigen::Vector3d e2 = V.row(F(face1, 2)) - V.row(F(face1, 0));

        double alpha1 = commone.dot(e1);
        double beta1 = t1.dot(e1);
        Eigen::Vector3d newe1 = alpha1*commone + beta1 * t2;
        Ts.block<2, 1>(2 * edgeidx, 0) = BTB2.inverse() * Bs[face2].transpose() * newe1;

        double alpha2 = commone.dot(e2);
        double beta2 = t1.dot(e2);
        Eigen::Vector3d newe2 = alpha2*commone + beta2*t2;
        Ts.block<2, 1>(2 * edgeidx, 1) = BTB2.inverse() * Bs[face2].transpose() * newe2;

        e1 = V.row(F(face2, 1)) - V.row(F(face2, 0));
        e2 = V.row(F(face2, 2)) - V.row(F(face2, 0));

        alpha1 = commone.dot(e1);
        beta1 = t2.dot(e1);
        newe1 = alpha1 * commone + beta1 * t1;
        Ts.block<2, 1>(2 * edgeidx, 2) = BTB1.inverse() * Bs[face1].transpose() * newe1;

        alpha2 = commone.dot(e2);
        beta2 = t2.dot(e2);
        newe2 = alpha2*commone + beta2*t1;
        Ts.block<2, 1>(2 * edgeidx, 3) = BTB1.inverse() * Bs[face1].transpose() * newe2;
    }
}

int Weave::vidx(int face, int field) const
{
    return (2 * nFields()*face + 2 * field);
}

Eigen::Vector2d Weave::v(int face, int field) const
{
    return vectorFields.segment<2>(vidx(face,field));
}

int Weave::betaidx(int face, int field) const
{
    return 2 * nFields()*nFaces() + 2 * nFields()*face + 2 * field;
}
Eigen::Vector2d Weave::beta(int face, int field) const
{
    return vectorFields.segment<2>(betaidx(face,field));
}

int Weave::alphaidx(int face, int field) const
{
    return 4 * nFields()*nFaces() + nFields()*face + field;
}

double Weave::alpha(int face, int field) const
{
    return vectorFields[alphaidx(face,field)];
}

void Weave::createVisualizationEdges(Eigen::MatrixXd &edgePts, Eigen::MatrixXd &edgeVecs, Eigen::MatrixXi &edgeSegs, Eigen::MatrixXd &colors)
{
    int nfaces = nFaces();
    int m = nFields();
    int nhandles = nHandles();
    edgePts.resize(m*nfaces + nhandles, 3);
    edgeVecs.resize(m*nfaces + nhandles, 3);
    edgeVecs.setZero();
    edgeSegs.resize(m*nfaces + nhandles, 2);
    colors.resize(m*nfaces + nhandles, 3);
    
    Eigen::MatrixXd fcolors(m, 3);
    for (int i = 0; i < m; i++)
        fcolors.row(i).setZero();//heatmap(double(i), 0.0, double(m-1));

    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector3d centroid;
        centroid.setZero();
        for (int j = 0; j < 3; j++)
            centroid += V.row(F(i, j));
        centroid /= 3.0;

        for (int j = 0; j < m; j++)
        {
            edgePts.row(m*i + j) = centroid;
            edgeVecs.row(m*i + j) = Bs[i] * v(i, j);
            edgeSegs(m*i + j, 0) = 2 * (m*i + j);
            edgeSegs(m*i + j, 1) = 2 * (m*i + j) + 1;
            colors.row(m*i + j) = fcolors.row(j);
        }
    }

    for (int i = 0; i < nhandles; i++)
    {
        Eigen::Vector3d centroid;
        centroid.setZero();
        for (int j = 0; j < 3; j++)
            centroid += V.row(F(handles[i].face, j));
        centroid /= 3.0;

        Eigen::Vector3d white(1, 1, 1);
        edgePts.row(m*nfaces + i) = centroid;
        edgeVecs.row(m*nfaces + i) = Bs[handles[i].face] * handles[i].dir;
        edgeSegs(m*nfaces + i, 0) = 2 * m*nfaces + 2 * i;
        edgeSegs(m*nfaces + i, 1) = 2 * m*nfaces + 2 * i + 1;
        colors.row(m*nfaces + i).setConstant(1.0);
    }
}

bool Weave::addHandle(Handle h)
{
    if (h.face < 0 || h.face > nFaces())
        return false;
    if(h.field < 0 || h.field > nFields())
        return false;

    Eigen::Vector3d extrinsic = Bs[h.face] * h.dir;
    double mag = extrinsic.norm();
    h.dir /= mag;
    handles.push_back(h);
    return true;
}

void Weave::normalizeFields()
{
    int nfaces = nFaces();
    int m = nFields();
    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < m; j++)
        {
            Eigen::Vector2d vif = v(i, j);
            double norm = sqrt(vif.transpose() * Bs[i].transpose() * Bs[i] * vif);
            vectorFields.segment<2>(vidx(i, j)) /= norm;
        }
    }
}