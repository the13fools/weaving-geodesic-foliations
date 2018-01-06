#include "WeaveHook.h"
#include "GaussNewton.h"
#include <iostream>

using namespace std;

bool WeaveHook::simulateOneStep()
{
    SolverParams params;
    params.lambdacompat = 100.0;
    //GNtestFiniteDifferences(*weave, params);
    //exit(-1);
    std::cout << "original energy: " << energy(*weave, params) << std::endl;;
    Eigen::VectorXd r;
    Eigen::SparseMatrix<double> J;
    GNEnergy(*weave, params, r);
    GNGradient(*weave, params, J);
    Eigen::SparseMatrix<double> M = J.transpose() * J;
    Eigen::SparseMatrix<double> I(weave->vectorFields.size(), weave->vectorFields.size());
    I.setIdentity();
    M += 1e-3 * I;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    solver.compute(M);
    Eigen::VectorXd rhs = J.transpose() * r;
    Eigen::VectorXd update = solver.solve(rhs);
    weave->vectorFields -= update;
    std::cout << "new energy: " << energy(*weave, params) << std::endl;;
    return false;
}