#ifndef LinearSolver_H
#define LinearSolver_H

#include <Eigen/Core>
#include <vector>
#include <Eigen/Sparse>

class Weave;
struct SolverParams;

struct Handle;

class LinearSolver
{
public:
    void addHandle(const Handle &h);
    void clearHandles();
    const std::vector<Handle> &getHandles() { return handles; }

    void handleConstraintOperator(const Weave &weave, SolverParams params, Eigen::VectorXd &primalVars, Eigen::SparseMatrix<double> &H, Eigen::VectorXd &h0);
    
    void curlOperator(const Weave &weave, SolverParams params, Eigen::SparseMatrix<double> &curlOp);
    void differentialOperator(const Weave &weave, SolverParams params, Eigen::SparseMatrix<double> &D);
    void differentialOperator_rosy(const Weave &weave, SolverParams params, Eigen::SparseMatrix<double> &D);

    void computeEnergy(const Weave &weave, SolverParams params, const Eigen::VectorXd &primalVars, const Eigen::VectorXd &dualVars, bool isRoSy );
    void updatePrimalVars(const Weave &weave, SolverParams params, Eigen::VectorXd &primalVars, Eigen::VectorXd &dualVars, bool isRoSy);
    void updateDualVars_new(const Weave &weave, SolverParams params, Eigen::VectorXd &primalVars, Eigen::VectorXd &dualVars, bool isRoSy);

    std::vector<Handle> handles;

private:
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > * dualSolver;
};

#endif
