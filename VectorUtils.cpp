

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
