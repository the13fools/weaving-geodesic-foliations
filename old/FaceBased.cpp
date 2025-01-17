#include "FaceBased.h"
#include <Eigen/Geometry>
#include "VectorUtils.h"


static void oneFaceGradientMatrix(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const Eigen::MatrixXi &F_edges, int edgeidx, Eigen::Vector3d &cdiff1, Eigen::Vector3d &cdiff2)
{
    cdiff1.setZero();
    cdiff2.setZero();
    //collect neighboring face indices
    int face1 = E(edgeidx, 2);
    int face2 = E(edgeidx, 3);

    //ignore boundary edges
    if (face1 == -1 || face2 == -1)
        return;

    int v1 = E(edgeidx, 0);
    int v2 = E(edgeidx, 1);

    Eigen::Vector3d n1 = faceNormal(F, V, face1);
    Eigen::Vector3d n2 = faceNormal(F, V, face2);
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
    
    cdiff1 = midpt + alpha * commone + beta * t1 - centroids[0];
    Eigen::Vector3d diff1 = centroids[0] - midpt;
    alpha = commone.dot(diff1);
    beta = t1.dot(diff1);
    cdiff2 = midpt + alpha*commone + beta * t2 - centroids[1];            
}

void computeGradientMatrices(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const Eigen::MatrixXi &F_edges, Eigen::MatrixXd &Ms)
{
    int nedges = (int)E.rows();
    Ms.resize(2 * nedges, 3);
    for (int i = 0; i < nedges; i++)
    {
        Eigen::Vector3d c1, c2;
        oneFaceGradientMatrix(F, V, E, F_edges, i, c1, c2);
        Ms.row(2 * i) = c1;
        Ms.row(2 * i + 1) = c2;
    }
}

Eigen::Vector3d projectOntoFace(const Eigen::Vector3d &v, const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, int face)
{
    Eigen::Vector3d n = faceNormal(F, V, face);
    return v - (v.dot(n))*n;
}

void computeBarycentricOperators(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, Eigen::MatrixXd &B, Eigen::SparseMatrix<double> &Bmat)
{
    int nfaces = (int)F.rows();
    B.resize(3 * nfaces, 2);
    B.setZero();
    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector3d p0 = V.row(F(i, 0));
        Eigen::Vector3d p1 = V.row(F(i, 1));
        Eigen::Vector3d p2 = V.row(F(i, 2));
        Eigen::Vector3d e1 = p1 - p0;
        Eigen::Vector3d e2 = p2 - p0;
        for (int j = 0; j < 3; j++)
        {
            B(3 * i + j, 0) = e1[j];
            B(3 * i + j, 1) = e2[j];
        }
    }

    Bmat.resize(3 * nfaces, 2 * nfaces);
    std::vector<Eigen::Triplet<double> > Bcoeffs;
    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Bcoeffs.push_back(Eigen::Triplet<double>(3 * i + j, 2 * i, B(3 * i + j, 0)));
            Bcoeffs.push_back(Eigen::Triplet<double>(3 * i + j, 2 * i + 1, B(3 * i + j, 1)));
        }
    }
    Bmat.setFromTriplets(Bcoeffs.begin(), Bcoeffs.end());
}

Eigen::Vector2d projectOntoBarycentric(const Eigen::Vector3d &v, const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, const Eigen::MatrixXd &B, int face)
{
    Eigen::MatrixXd bblock = B.block<3, 2>(3 * face, 0);
    Eigen::Matrix2d BTB = bblock.transpose()*bblock;
    return BTB.inverse() * bblock.transpose() * v;
}
