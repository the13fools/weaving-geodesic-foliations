#ifndef GAUSSNEWTON_H
#define GAUSSNEWTON_H

#include <Eigen/Core>
#include <Eigen/Sparse>

class Weave;

struct SolverParams
{
    double lambdacompat; // weight of compatibility term
    double lambdareg;    // Tilhonov regularization
};

void GNmetric(const Weave &weave, Eigen::SparseMatrix<double> &M);
void GNEnergy(const Weave &weave, SolverParams params, Eigen::VectorXd &E);
void GNGradient(const Weave &weave, SolverParams params, Eigen::SparseMatrix<double> &J);

void GNtestFiniteDifferences(Weave &weave, SolverParams params);

void oneStep(Weave &weave, SolverParams params);

/*
 * Computes |F| x m matrix of face energies due to vector field derivative incompatibility, contributed by each vector field on each face.
 */
void faceEnergies(const Weave &weave, SolverParams params, Eigen::MatrixXd &E);
#endif
