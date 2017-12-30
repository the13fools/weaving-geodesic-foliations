#ifndef FACE_BASED
#define FACE_BASED

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

void computeGradientMatrices(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const Eigen::MatrixXi &F_edges, Eigen::MatrixXd &Ms);

#endif