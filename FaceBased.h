#ifndef FACE_BASED
#define FACE_BASED

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

void computeGradientMatrices(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const Eigen::MatrixXi &F_edges, Eigen::MatrixXd &Ms);
void computeBarycentricOperators(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, Eigen::MatrixXd &B, Eigen::SparseMatrix<double> &Bmat);

Eigen::Vector3d projectOntoFace(const Eigen::Vector3d &v, const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, int face);
Eigen::Vector2d projectOntoBarycentric(const Eigen::Vector3d &v, const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, const Eigen::MatrixXd &B, int face);

#endif