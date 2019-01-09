#ifndef GNGLOBALINTEGRATION_H
#define GNGLOBALINTEGRATION_H

#include "FieldIntegration.h"

// Our method, based on Gauss-Newton optimization of theta and s

class GNGlobalIntegration : public  GlobalFieldIntegration
{
public:
    GNGlobalIntegration(int alternationIters, int powerIters) : outerIters_(alternationIters), powerIters_(powerIters) {}

    void globallyIntegrateOneComponent(const Surface &surf, const Eigen::MatrixXd &v, Eigen::VectorXd &scales, Eigen::VectorXd &theta);
    
private:
    int outerIters_;
    int powerIters_;
};

#endif
