#ifndef GNGLOBALINTEGRATION_H
#define GNGLOBALINTEGRATION_H

#include "FieldIntegration.h"

// Our method, based on Gauss-Newton optimization of theta and s

class GNGlobalIntegration : public  GlobalFieldIntegration
{
public:
    GNGlobalIntegration() {}

    void globallyIntegrateOneComponent(const Surface &surf, const Eigen::MatrixXd &v, Eigen::VectorXd &scales, Eigen::VectorXd &theta);
};

#endif