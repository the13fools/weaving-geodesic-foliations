

#include "VectorUtils.h"
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


// This function computes the discrete connection. 
// Conceptually, consider that this unfolds two adjacent triangles, 
// shifts a vector across the common edge, and then refolds the edge.
Eigen::Vector3d mapVectorToAdjacentFace(const Eigen::MatrixXi &F, 
                                        const Eigen::MatrixXd &V,
                                        const Eigen::MatrixXi &E, 
                                        int edgeIdx, int init_fidx, int neighbor_fidx,
                                        const Eigen::Vector3d &inp)
{ 
    Eigen::Vector3d n1 = faceNormal(F, V, init_fidx);
    Eigen::Vector3d n2 = faceNormal(F, V, neighbor_fidx);
    Eigen::Vector3d commone = V.row( E(edgeIdx, 0) ) - V.row( E(edgeIdx, 1) );
    commone.normalize();

    Eigen::Vector3d t1 = n1.cross(commone);
    Eigen::Vector3d t2 = n2.cross(commone);

    double alpha = commone.dot( inp );
    double beta  = t1.dot( inp );

    return  alpha * commone + beta * t2;
}
