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
    std::vector< Eigen::VectorXd > bending;

    void loadSampledCurves(const std::string &filename);

    void traceCurve(const Weave &vw, const Trace_Mode trace_state, int traceIdx, int sign, int faceId, int steps);
    void popLastCurve();
    void logRibbonsToFile(std::string foldername, std::string filename);
    void computeIntersections(int curveIdx1, int curveIdx2, std::vector<Collision> &collisions);
    void save(const std::string &filename);
    void load(const std::string &filename);
};


#endif
