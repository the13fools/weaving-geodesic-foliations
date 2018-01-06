#ifndef GAUSSNEWTON_H
#define GAUSSNEWTON_H

#include <Eigen/Core>
#include <Eigen/Sparse>

class Weave;

struct SolverParams
{
    double lambdacompat; // weight of compatibility term
};

double energy(const Weave &weave, SolverParams params);
void trueGradient(const Weave &weave, SolverParams params, Eigen::VectorXd &dE);

void GNEnergy(const Weave &weave, SolverParams params, Eigen::VectorXd &E);
void GNGradient(const Weave &weave, SolverParams params, Eigen::SparseMatrix<double> &J);

void testFiniteDifferences(Weave &weave, SolverParams params);
void GNtestFiniteDifferences(Weave &weave, SolverParams params);
#endif
