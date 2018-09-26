#ifndef GAUSSNEWTON_H
#define GAUSSNEWTON_H

#include <Eigen/Core>
#include <Eigen/Sparse>

class Weave;

struct SolverParams
{
    double lambdacompat; // weight of compatibility term
    double lambdareg;    // Tilhonov regularization
    double curlreg; // Weight on the curl component of v
    int eigenvector; // Choose the eigenvector of s to display (for debugging - delete me!)
    Eigen::VectorXd edgeWeights;
};

void GNmetric(const Weave &weave, Eigen::SparseMatrix<double> &M);
void GNEnergy(const Weave &weave, SolverParams params, Eigen::VectorXd &E);
void GNGradient(const Weave &weave, SolverParams params, Eigen::SparseMatrix<double> &J);

void GNtestFiniteDifferences(Weave &weave, SolverParams params);
double lineSearch(Weave &weave, SolverParams params, const Eigen::VectorXd &update);
void oneStep(Weave &weave, SolverParams params);

/*
 * Computes |F| x m matrix of face energies due to vector field derivative incompatibility, contributed by each vector field on each face.
 */
void faceEnergies(const Weave &weave, SolverParams params, Eigen::MatrixXd &E);
#endif
