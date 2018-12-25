#ifndef MIGLOBALINTEGRATION_H
#define MIGLOBALINTEGRATION_H

#include "FieldIntegration.h"
#include <Eigen/Sparse>

class MIGlobalIntegration : public GlobalFieldIntegration
{
public:
    MIGlobalIntegration(double anisotropy, double smoothnessReg) :
        aniso_(anisotropy),
        smoothreg_(smoothnessReg)
    {}

    virtual void globallyIntegrateOneComponent(const Surface &surf, const Eigen::MatrixXd &v, Eigen::VectorXd &s, Eigen::VectorXd &theta);

private:
    double aniso_;
    double smoothreg_;
};

#endif