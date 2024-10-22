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
    double handleScale;
    bool softHandleConstraint; // in the Linear Solver, controls the handle mode
    double vizVectorCurl; // in field surface, vizualization variable
    double vizCorrectionCurl; // in field surface, vizualization variable
    bool vizNormalizeVecs;
    bool vizShowCurlSign; // round the viz to see the sign.
    bool disableCurlConstraint;
    int rosyN;
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
