#ifndef TRACE_H
#define TRACE_H

#include <Eigen/Core>
#include <Eigen/Sparse>

class Weave;

struct Collision {
    Collision(int rod1, int rod2, int seg1, int seg2);
    int rod1;
    int rod2; // can be the same rod 
    int seg1; 
    int seg2;
};

class Trace
{
public:
    Trace();
    ~Trace();

    // eeeeeew
    std::vector< Eigen::MatrixXd > curves;
    std::vector< Eigen::MatrixXd > normals;
        
    void traceCurve(const Weave &vw, const Eigen::Vector3d dir, int faceId, int steps);
    void popLastCurve();
};


#endif
