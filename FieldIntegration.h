#ifndef FIELDINTEGRATION_H
#define FIELDINTEGRATION_H

#include "Surface.h"
#include <Eigen/Core>

class LocalFieldIntegration
{
public:
    // Locally integrates a given vector field, assuming:
    // - the surface surf has one connected component
    // - the vector field v has no singularities
    // Result is a rescaling s on the faces of surf.
    virtual void locallyIntegrateOneComponent(const Surface &surf, const Eigen::MatrixXd &v, Eigen::VectorXd &s) = 0;
};

class GlobalFieldIntegration
{
public:
    // Integrates a given vector field, assuming:
    // - the surface surf has one connected component
    // - the vector field v has no singularities
    // Result is a periodic function (values in [0, 2pi)) on the vertices of surf.
    virtual void globallyIntegrateOneComponent(const Surface &surf, const Eigen::MatrixXd &v, Eigen::VectorXd &s, Eigen::VectorXd &theta) = 0;
};

// does nothing except normalize the vector field
class TrivialLocalIntegration : public LocalFieldIntegration
{
public:
    void locallyIntegrateOneComponent(const Surface &surf, const Eigen::MatrixXd &v, Eigen::VectorXd &s)
    {
        int nfaces = surf.nFaces();
        s.resize(nfaces);
        for (int i = 0; i < nfaces; i++)
        {
            s[i] = 1.0;
        }
    }
};

#endif