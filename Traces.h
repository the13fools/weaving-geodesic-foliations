#ifndef TRACES_H
#define TRACES_H

#include <vector>
#include <Eigen/Core>
class FieldSurface;

enum Trace_Mode {
    GEODESIC = 0,
    FIELD    
};

struct Collision {
    Collision(int trace1, int trace2, int seg1, int seg2);
    int trace1;
    int trace2; // can be the same rod 
    int seg1; 
    int seg2;
};

// one piece of a trace
struct TraceSegment
{
    int face; // index into the faces of the CoverMesh
    int side[2]; // integer [0,3]; the first endpoint of the segment lies on edge fs->data().faceEdges(side[0]), etc
    double bary[2]; // barycentric coordinates of segment endpoint along each edge
};

// trace in its exact form (as segments interpolating points on a surface's edges)
class Trace
{
public:
    Trace(const FieldSurface *parent, Trace_Mode type) : parent_(parent), type_(type) {}

    std::vector<TraceSegment> segs;    
    Trace_Mode type_;
    const FieldSurface *parent_;
};

// trace in its "rationalized" form (loose points and normals in space)
struct RationalizedTrace
{
    Eigen::MatrixXd pts;
    Eigen::MatrixXd normals;
};

// a collection of traces on a surface
class TraceSet
{
public:
    int nTraces() const { return traces_.size(); }
    const Trace &trace(int id) const { return traces_[id]; }
    
    void addTrace(const Trace &tr);
    void traceCurve(const FieldSurface &parent, const Trace_Mode trace_state, int traceIdx, int sign, int faceId, int steps);

    // delete all traces
    void clear();
    // delete only traces belonging to specific surface
    void purgeTraces(const FieldSurface *surface);

    void popLastCurve();



    // converts a trace to a set of points and normals; does *not* do any cleanup (just converts segments as-they-are)
    void renderTrace(int traceid, std::vector<Eigen::Vector3d> &verts, std::vector<Eigen::Vector3d> &normals) const;
private:
    void findTraceStart(const FieldSurface &parent,
        int curr_face_id,
        const Eigen::Vector3d startpoint,
        const Eigen::Vector3d curr_dir, int &startEdge, double &startBery);

    void getNextTracePoint(const FieldSurface &parent,
        int curr_face_id,
        int curr_edge_id,
        double curr_bary,
        const Eigen::Vector3d curr_dir,
        int &next_edge_id,
        double &next_bary,
        int &opp_face_id,
        int &opp_edge_id);

    Eigen::Vector3d pointFromBary(const FieldSurface &parent, int faceId, int faceEdge, double bary);

    std::vector<Trace> traces_;    
};

#endif
