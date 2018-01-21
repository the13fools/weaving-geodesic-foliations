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

struct TracePoint {
    int face_id; // current face to trace
    int edge_id; // current edge that contains point
    int op_v_id;
    Eigen::Vector3d n; // Face normal
    Eigen::Vector3d point; // current point
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
    void logRibbonsToFile(std::string foldername, std::string filename, const Weave &wv);
    void computeIntersections(int curveIdx1, int curveIdx2, std::vector<Collision> &collisions, std::vector<Eigen::MatrixXd> &curves);
    void save(const std::string &filename);
    void load(const std::string &filename);

    void getNextTracePoint(const Weave &wv, 
                           int curr_face_id, int curr_edge_id, 
                           const Eigen::Vector3d prev_point, int op_v_id,  
                           const Eigen::Vector3d curr_dir, TracePoint &nextTrace);

    void startTraceFromPoint(const Weave &wv,
                             int curr_face_id, 
                             const Eigen::Vector3d startpoint, 
                             const Eigen::Vector3d curr_dir, TracePoint &startPoint);


};


#endif
