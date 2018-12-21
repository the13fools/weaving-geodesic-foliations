#ifndef BOMMESFIELDINTEGRATION_H
#define BOMMESFIELDINTEGRATION_H

#include "FieldIntegration.h"
#include <Eigen/Sparse>

class BommesFieldIntegration : public FieldIntegration
{
public:
    BommesFieldIntegration(double anisotropy, double smoothnessReg, double globalScale) :
        aniso_(anisotropy),
        smoothreg_(smoothnessReg),
        globalScale_(globalScale)
    {}

    virtual void integrateOneComponent(const Surface &surf, const Eigen::MatrixXd &v, Eigen::VectorXd &theta);

private:
    void ComisoWrapper(const Eigen::SparseMatrix<double> &constraints,
        const Eigen::SparseMatrix<double> &A,
        Eigen::VectorXd &result,
        const Eigen::VectorXd &rhs,
        const Eigen::VectorXi &toRound,
        double reg);

    double aniso_;
    double smoothreg_;
    double globalScale_;
};

#endif