#ifndef FIELDINTEGRATION_H
#define FIELDINTEGRATION_H

#include "Surface.h"
#include <Eigen/Core>

class FieldIntegration
{
public:
    // Integrates a given vector field, assuming:
    // - the surface s has one connected component
    // - the vector field v has no singularities
    // Result is a periodic function (values in [0, 2pi)) on the vertices of s.
    virtual void integrateOneComponent(const Surface &surf, const Eigen::MatrixXd &v, Eigen::VectorXd &theta) = 0;
};

#endif