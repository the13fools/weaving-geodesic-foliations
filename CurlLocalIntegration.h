#ifndef CURLLOCALINTEGRATION_H
#define CURLLOCALINTEGRATION_H

#include "FieldIntegration.h"

// Compute s to minimize curl (Ray et al)
class CurlLocalIntegration : public LocalFieldIntegration
{
public:
    CurlLocalIntegration(double sSmoothnessReg) : sreg_(sSmoothnessReg) {}

    void locallyIntegrateOneComponent(const Surface &surf, const Eigen::MatrixXd &v, Eigen::VectorXd &s);

private:
    double sreg_;
};

#endif