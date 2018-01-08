#ifndef TRACE_H
#define TRACE_H

#include <Eigen/Core>
#include <Eigen/Sparse>

class Weave;

class Trace
{
    void traceCurve(const Weave &md,
                const Eigen::Vector3d dir, int faceId, 
                std::vector<Eigen::Vector3d> &curve,
                std::vector<Eigen::Vector3d> &normal);
};


#endif
