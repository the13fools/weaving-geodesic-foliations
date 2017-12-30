#include "FaceBased.h"
#include <Eigen/Geometry>

Eigen::Vector3d faceNormal(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, int faceidx)
{
    Eigen::Vector3d p0 = V.row(F(faceidx, 0));
    Eigen::Vector3d p1 = V.row(F(faceidx, 1));
    Eigen::Vector3d p2 = V.row(F(faceidx, 2));
    Eigen::Vector3d n = (p1 - p0).cross(p2 - p0);
    n /= n.norm();
    return n;
}

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

