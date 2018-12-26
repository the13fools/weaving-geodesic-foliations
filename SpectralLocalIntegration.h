#ifndef SPECTRALLOCALINTEGRATION_H
#define SPECTRALLOCALINTEGRATION_H

#include "FieldIntegration.h"

// Our method that estimates s using eigenvalue problem
class SpectralLocalIntegration : public LocalFieldIntegration
{
public:
    SpectralLocalIntegration(double sSmoothnessReg) : sreg_(sSmoothnessReg) {}

    void locallyIntegrateOneComponent(const Surface &surf, const Eigen::MatrixXd &v, Eigen::VectorXd &s);

private:
    double sreg_;
};

#endif