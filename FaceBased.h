#ifndef FACE_BASED
#define FACE_BASED

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

/*
Computes a set of 3 x F matrices that, when multiplied by a discrete face-based function, returns the per-face gradient (in world coordinates) of that function.
For now, ignores boundary triangles (the matrix returned for these is just the zero matrix)
*/
void computeGradientMatrices(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const Eigen::MatrixXi &F_edges, std::vector<Eigen::SparseMatrix<double> > &Ms);

#endif