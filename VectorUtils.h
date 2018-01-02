#ifndef VECTOR_UTILS
#define VECTOR_UTILS

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>


Eigen::Vector3d faceNormal(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, int faceidx);


#endif
