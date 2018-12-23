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
    double aniso_;
    double smoothreg_;
    double globalScale_;
};

#endif