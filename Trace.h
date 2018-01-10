#ifndef TRACE_H
#define TRACE_H

#include <Eigen/Core>
#include <Eigen/Sparse>

class Weave;

enum Trace_Mode {
    GEODESIC = 0,
    FIELD    
};

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

    std::vector< Eigen::MatrixXd > curves;
    std::vector< Eigen::MatrixXd > normals;    
    std::vector< Trace_Mode > modes;

    void traceCurve(const Weave &vw, const Trace_Mode trace_state, const Eigen::Vector3d dir, int faceId, int steps);
    void popLastCurve();
    void logRibbonsToFile(std::string foldername, std::string filename);
    void computeIntersections(int curveIdx1, int curveIdx2, std::vector<Collision> &collisions);

};


#endif
