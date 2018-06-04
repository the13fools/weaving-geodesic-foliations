#ifndef TRACES_H
#define TRACES_H

#include <vector>
#include <Eigen/Core>
#include <set>

class FieldSurface;

enum Trace_Mode {
    GEODESIC = 0,
    FIELD    
};

struct Collision {    
    int rod1;
    int rod2; // can be the same rod 
    int seg1; 
    int seg2;
    double bary1;
    double bary2;
};

// one piece of a trace
struct TraceSegment
{
    int face; // index into the faces of the CoverMesh
    int side[2]; // integer [0,3]; the first endpoint of the segment lies on edge fs->data().faceEdges(side[0]), etc
    double bary[2]; // barycentric coordinates of segment endpoint along each edge
};

struct TraceCollision
{
    int seg1, seg2;
    double bary1, bary2;
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
    
    int nRationalizedTraces() const { return rattraces_.size(); }
    const RationalizedTrace &rationalizedTrace(int id) const { return rattraces_[id]; }

    int nCollisions() const { return collisions_.size(); }
    
    void addTrace(const Trace &tr);
    void traceCurve(const FieldSurface &parent, const Trace_Mode trace_state, int traceIdx, int sign, int faceId, int steps);

    // delete all traces
    void clear();
    // delete only traces belonging to specific surface
    void purgeTraces(const FieldSurface *surface);

    void popLastCurve();

    void rationalizeTraces(double maxcurvature, double extenddist, double seglen, double minlen);

    // converts a trace to a set of points and normals; does *not* do any cleanup (just converts segments as-they-are)
    void renderTrace(int traceid, std::vector<Eigen::Vector3d> &verts, std::vector<Eigen::Vector3d> &normals) const;

    void collisionPoint(int collision, Eigen::Vector3d &pt0, Eigen::Vector3d &pt1) const;

    void exportRodFile(const char*filename);
    void exportForRendering(const char *filename);

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
        int &opp_edge_id) const;

    void findCurvedVerts(const Trace &tr, double maxcurvature, std::set<int> &badverts) const;
    void splitTrace(const Trace &tr, const std::set<int> &badvverts, std::vector<Trace> &splittr) const;
    // attempts to extend a trace in both directions by the given distance. Returns the amount it actually succeeded in extending (could be smaller or larger)
    void extendTrace(Trace &tr, double extbeginning, double extend, double &actualextbeginning, double &actualextend) const;

    double arclength(const Trace &tr) const;
    void sampleTrace(const Trace &tr, double start, double end, int nsegs, RationalizedTrace &rattrace, std::vector<double> &samples);
    void findPointOnTrace(const Trace &tr, double s, int &seg, double &bary);
    void computeIntersections(const Trace &tr1, const Trace &tr2, bool selfcollision,
        std::vector<TraceCollision> &collisions) const;

    Eigen::Vector3d pointFromBary(const FieldSurface &parent, int faceId, int faceEdge, double bary) const;
    double geodesicCurvature(const Eigen::Vector3d &p0, const Eigen::Vector3d &p1, const Eigen::Vector3d &p2, const Eigen::Vector3d &n) const;

    std::vector<Trace> traces_;    

    std::vector<RationalizedTrace> rattraces_;
    std::vector<Collision> collisions_;
};

#endif
