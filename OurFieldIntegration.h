#ifndef OURFIELDINTEGRATION_H
#define OURFIELDINTEGRATION_H

#include "FieldIntegration.h"

class OurFieldIntegration : public FieldIntegration
{
public:
    OurFieldIntegration(double sSmoothnessReg, double globalScale) : 
        sreg_(sSmoothnessReg),
        globalScale_(globalScale)
    {}

    virtual void integrateOneComponent(const Surface &surf, const Eigen::MatrixXd &v, Eigen::VectorXd &theta);

private:
    void initializeS(const Surface &surf, const Eigen::MatrixXd &v, Eigen::VectorXd &s);
    void globalSThetaSolve(const Surface &surf, const Eigen::MatrixXd &v, Eigen::VectorXd &s, Eigen::VectorXd &theta);

    double sreg_;
    double globalScale_;
};

#endif