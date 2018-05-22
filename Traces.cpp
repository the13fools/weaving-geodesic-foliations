#include "Traces.h"
#include "FieldSurface.h"
#include "Distance.h"
#include <algorithm>

void TraceSet::addTrace(const Trace &tr)
{
    traces_.push_back(tr);
}

void TraceSet::renderTrace(int traceid, std::vector<Eigen::Vector3d> &verts, std::vector<Eigen::Vector3d> &normals) const
{
    int nsegs = traces_[traceid].segs.size();
    int nverts = nsegs + 1;
    verts.resize(nverts);
    for (int i = 0; i < nverts; i++)
        verts[i].setZero();
    normals.resize(nsegs);
    const FieldSurface &parent = *traces_[traceid].parent_;
    for (int i = 0; i < nsegs; i++)
    {
        int v0 = parent.data().F(traces_[traceid].segs[i].face, (traces_[traceid].segs[i].side[0]+1)%3);
        int v1 = parent.data().F(traces_[traceid].segs[i].face, (traces_[traceid].segs[i].side[0]+2)%3);
        Eigen::Vector3d pos = (1.0 - traces_[traceid].segs[i].bary[0])*parent.data().V.row(v0).transpose() + traces_[traceid].segs[i].bary[0] * parent.data().V.row(v1).transpose();
        verts[i] += pos;

        v0 = parent.data().F(traces_[traceid].segs[i].face, (traces_[traceid].segs[i].side[1]+1)%3);
        v1 = parent.data().F(traces_[traceid].segs[i].face, (traces_[traceid].segs[i].side[1]+2)%3);
        pos = (1.0 - traces_[traceid].segs[i].bary[1])*parent.data().V.row(v0).transpose() + traces_[traceid].segs[i].bary[1] * parent.data().V.row(v1).transpose();
        verts[i + 1] += pos;
    }
    for (int i = 1; i < nverts - 1; i++)
        verts[i] /= 2.0;

    for (int i = 0; i < nsegs; i++)
    {
        normals[i] = parent.faceNormal(traces_[traceid].segs[i].face);
    }
}

void TraceSet::clear()
{
    traces_.clear();
}

void TraceSet::purgeTraces(const FieldSurface *surface)
{
    std::vector<Trace> newtraces;
    for (auto it : traces_)
    {
        if (it.parent_ != surface)
            newtraces.push_back(it);
    }
    traces_ = newtraces;
}

void TraceSet::popLastCurve()
{
    if (!traces_.empty())
        traces_.pop_back();
}


void TraceSet::findTraceStart(const FieldSurface &parent,
    int curr_face_id, 
    const Eigen::Vector3d startpoint, 
    const Eigen::Vector3d curr_dir, int &startEdge, double &startBary)
{

    // Project backwards to initialize.  Kinda hacky.
    double min_dist = std::numeric_limits<double>::infinity();

    startEdge = -1;    

    for (int i = 0; i < 3; i++)
    {        
        Eigen::Vector3d op_v1 = parent.data().V.row(parent.data().F(curr_face_id, (i + 1) % 3)).transpose();
        Eigen::Vector3d op_v2 = parent.data().V.row(parent.data().F(curr_face_id, (i + 2) % 3)).transpose();
        
        double p0bary, p1bary, q0bary, q1bary;
        Eigen::Vector3d dist = Distance::edgeEdgeDistance(startpoint,
            startpoint - curr_dir,
            op_v1, op_v2,
            p0bary, p1bary, q0bary, q1bary);
        if (dist.norm() < min_dist)
        {
            min_dist = dist.norm();
            startBary = q1bary;
            startEdge = i;
        }
    }
}

void TraceSet::getNextTracePoint(const FieldSurface &parent,
    int curr_face_id,
    int curr_edge_id,
    double curr_bary,
    const Eigen::Vector3d curr_dir,
    int &next_edge_id,
    double &next_bary,
    int &opp_face_id,
    int &opp_edge_id) const
{
    Eigen::Vector3d op_vertex = parent.data().V.row(parent.data().F(curr_face_id,curr_edge_id)).transpose();
    Eigen::Vector3d prev_point = pointFromBary(parent, curr_face_id, curr_edge_id, curr_bary);
    Eigen::Vector3d split = op_vertex - prev_point;

    double split_len = split.norm();
    split = split / split_len;
    Eigen::Vector3d n = parent.faceNormal(curr_face_id);
    Eigen::Vector3d perp = split.cross(n);    

    int op_edge_id = -1;
    for (int j = 0; j < 3; j++)
    {
        if (j == curr_edge_id)
            continue;
        int v0 = parent.data().F(curr_face_id, (j + 1) % 3);
        int v1 = parent.data().F(curr_face_id, (j + 2) % 3);
            
        Eigen::Vector3d e_test = (parent.data().V.row(v0) + parent.data().V.row(v1)) * .5;
        e_test -= prev_point;
        if (e_test.dot(perp) * curr_dir.dot(perp) > 0.)
        {
            op_edge_id = j;
            break;
        }
    }
    // stop if we hit vertex
    if (op_edge_id == -1)
    {
        next_edge_id = -1;
        next_bary = 0.0;
        opp_face_id = -1;
        opp_edge_id = -1;
        return;
    }

    // Find intersection point.
    next_edge_id = op_edge_id;
    Eigen::Vector3d op_v1 = parent.data().V.row(parent.data().F(curr_face_id, (next_edge_id + 1) % 3));
    Eigen::Vector3d op_v2 = parent.data().V.row(parent.data().F(curr_face_id, (next_edge_id + 2) % 3));
    
    double p0bary, p1bary, q0bary, q1bary;
    Eigen::Vector3d dist = Distance::edgeEdgeDistance(prev_point,
        prev_point + curr_dir,
        op_v1, op_v2,
        p0bary, p1bary, q0bary, q1bary);

    if (dist.norm() > 1e-6)
    {
        next_edge_id = -1;
        next_bary = 0.0;
        opp_face_id = -1;
        opp_edge_id = -1;
        return;
    }
    next_bary = q1bary;

    opp_face_id = parent.data().E(parent.data().faceEdges(curr_face_id, next_edge_id), 0);
    if (opp_face_id == curr_face_id)
    {
        opp_face_id = parent.data().E(parent.data().faceEdges(curr_face_id, next_edge_id), 1);
    }    

    if (opp_face_id != -1)
    {
        for (int i = 0; i < 3; i++)
        {
            if (parent.data().faceEdges(curr_face_id, next_edge_id) == parent.data().faceEdges(opp_face_id, i))
                opp_edge_id = i;
        }
    }    
}


void TraceSet::traceCurve(const FieldSurface &parent, const Trace_Mode trace_state,
    int traceIdx, int sign, int faceId, int steps)
{
    int curr_dir_idx = abs(traceIdx);
    double coeff_dir = sign > 0 ? 1 : -1;
    Eigen::Vector3d curr_dir = coeff_dir * parent.data().Bs[faceId] * parent.v(faceId, curr_dir_idx);

    // assumes roughly delaunay    
    curr_dir = curr_dir.normalized() * parent.data().averageEdgeLength * 1000.;

    int curr_face_id = faceId;

    Eigen::Vector3d startpoint;
    startpoint.setZero();
    for (int i = 0; i < 3; i++)
    {
        startpoint += parent.data().V.row(parent.data().F(faceId, i));
    }
    startpoint /= 3.0;


    int curr_edge_id;
    double curr_bary;
    findTraceStart(parent, curr_face_id, startpoint, curr_dir, curr_edge_id, curr_bary); 
    
    Trace t(&parent, trace_state);

    assert(curr_edge_id > -1);
    for (int i = 0; i < steps; i++)
    {
        int next_edge_id;
        double next_bary;
        int opp_face_id;
        int opp_edge_id;
        getNextTracePoint(parent, curr_face_id, curr_edge_id, curr_bary, curr_dir, next_edge_id, next_bary, opp_face_id, opp_edge_id);  
        
        if (next_edge_id == -1)
            break;
        if (opp_face_id == -1)
            break;

        TraceSegment seg;
        seg.face = curr_face_id;
        seg.side[0] = curr_edge_id;
        seg.side[1] = next_edge_id;
        seg.bary[0] = curr_bary;
        seg.bary[1] = next_bary;
        t.segs.push_back(seg);

        switch (trace_state)
        {
        case GEODESIC:
        {
            // projet vector into barycentric coordinates
            Eigen::Matrix<double, 3, 2> B = parent.data().Bs[curr_face_id];
            Eigen::Vector2d barys = (B.transpose()*B).inverse() * B.transpose() * curr_dir;
            int edgeid = parent.data().faceEdges(curr_face_id, next_edge_id);
            int side = 0;
            if (parent.data().E(edgeid, side) != curr_face_id)
                side = 1;
            Eigen::Vector2d nextbarys = parent.data().Ts.block<2, 2>(2 * edgeid, 2 * side) * barys;
            curr_dir = parent.data().Bs[opp_face_id] * nextbarys;
        }
        break;
        case FIELD:
        {
            int edgeid = parent.data().faceEdges(curr_face_id, next_edge_id);
            Eigen::MatrixXi perm = parent.Ps(edgeid);
            if (opp_face_id == parent.data().E(edgeid, 1))
            {
                perm.transposeInPlace();
            }
            Eigen::VectorXi curr_vec(parent.nFields());
            curr_vec.setZero();
            curr_vec(curr_dir_idx) = 1;
            Eigen::VectorXi next_vec = perm * curr_vec;            
            for (int idx = 0; idx < parent.nFields(); idx++)
            {
                if (next_vec(idx) != 0)
                {
                    if (next_vec(idx) < 0)
                    {
                        coeff_dir *= -1.;
                    }
                    curr_dir = parent.data().Bs[opp_face_id] * parent.v(opp_face_id, idx);
                    curr_dir_idx = idx;
                    break;
                }
            }
            curr_dir = coeff_dir * curr_dir.normalized() * parent.data().averageEdgeLength * 1000.;
            
        }
            break;
        }

        curr_face_id = opp_face_id;
        curr_edge_id = opp_edge_id;
        curr_bary = 1.0 - next_bary;                
    }

    addTrace(t);    
}

Eigen::Vector3d TraceSet::pointFromBary(const FieldSurface &parent, int faceId, int faceEdge, double bary) const
{
    Eigen::Vector3d v0 = parent.data().V.row(parent.data().F(faceId, (faceEdge + 1) % 3)).transpose();
    Eigen::Vector3d v1 = parent.data().V.row(parent.data().F(faceId, (faceEdge + 2) % 3)).transpose();
    return (1.0 - bary)*v0 + bary*v1;
}

static double angle(Eigen::Vector3d v, Eigen::Vector3d w, Eigen::Vector3d n)
{
    return 2.0 * atan2( v.cross(w).dot(n), v.norm() * w.norm() + v.dot(w));
}

void TraceSet::findCurvedVerts(const Trace &tr, double maxcurvature, std::set<int> &badverts) const
{
    badverts.clear();
    int nsegs = tr.segs.size();
    for (int i = 0; i < nsegs - 1; i++)
    {
        Eigen::Vector3d v0 = pointFromBary(*tr.parent_, tr.segs[i].face, tr.segs[i].side[0], tr.segs[i].bary[0]);
        Eigen::Vector3d v1 = pointFromBary(*tr.parent_, tr.segs[i].face, tr.segs[i].side[1], tr.segs[i].bary[1]);
        v1 += pointFromBary(*tr.parent_, tr.segs[i + 1].face, tr.segs[i + 1].side[0], tr.segs[i + 1].bary[0]);
        v1 /= 2.0;
        Eigen::Vector3d v2 = pointFromBary(*tr.parent_, tr.segs[i + 1].face, tr.segs[i + 1].side[1], tr.segs[i + 1].bary[1]);
        Eigen::Vector3d n = tr.parent_->faceNormal(tr.segs[i].face);
        n += tr.parent_->faceNormal(tr.segs[i + 1].face);
        n.normalize();
        double theta = fabs(angle(v1 - v0, v2 - v1, n));
        if (theta > maxcurvature)
            badverts.insert(i+1);
    }
}

void TraceSet::splitTrace(const Trace &tr, const std::set<int> &badverts, std::vector<Trace> &splittr) const
{
    splittr.clear();
    std::vector<int> splitsegs;
    for (auto it : badverts)
        splitsegs.push_back(it-1);
    std::sort(splitsegs.begin(), splitsegs.end());
    splittr.push_back(Trace(tr.parent_, tr.type_));
    int nsegs = tr.segs.size();
    splitsegs.push_back(nsegs);
    int idx = 0;
    for (int i = 0; i < nsegs; i++)
    {
        splittr.back().segs.push_back(tr.segs[i]);
        if (i == splitsegs[idx])
        {
            splittr.push_back(Trace(tr.parent_, tr.type_));
            idx++;
        }
    }
}

void TraceSet::extendTrace(Trace &tr, double extbeginning, double extend, double &actualextbeginning, double &actualextend) const
{
    actualextbeginning = 0;
    actualextend = 0;
    if (tr.segs.size() == 0)
        return;

    // extend beginning
    std::vector<TraceSegment> toprepend;
    int begface = tr.segs[0].face;
    int begedge = tr.segs[0].side[1];
    double begbary = tr.segs[0].bary[1];
    Eigen::Vector3d v0 = pointFromBary(*tr.parent_, begface, begedge, begbary);
    Eigen::Vector3d v1 = pointFromBary(*tr.parent_, begface, tr.segs[0].side[0], tr.segs[0].bary[0]);
    Eigen::Vector3d dir = v1 - v0;
    dir *= 1e3 / dir.norm();
    int nextedge, oppface, oppedge;
    double nextbary;
    getNextTracePoint(*tr.parent_, begface, begedge, begbary, dir, nextedge, nextbary, oppface, oppedge);    
    while (extbeginning > 0)
    {
        if (nextedge == -1 || oppface == -1)
        {
            // failed to extend
            break;
        }

        // projet vector into barycentric coordinates
        Eigen::Matrix<double, 3, 2> B = tr.parent_->data().Bs[begface];
        Eigen::Vector2d barys = (B.transpose()*B).inverse() * B.transpose() * dir;
        int edgeid = tr.parent_->data().faceEdges(begface, nextedge);
        int side = 0;
        if (tr.parent_->data().E(edgeid, side) != begface)
            side = 1;
        Eigen::Vector2d nextbarys = tr.parent_->data().Ts.block<2, 2>(2 * edgeid, 2 * side) * barys;
        dir = tr.parent_->data().Bs[oppface] * nextbarys;

        begface = oppface;
        begedge = oppedge;
        begbary = 1.0 - nextbary;

        getNextTracePoint(*tr.parent_, begface, begedge, begbary, dir, nextedge, nextbary, oppface, oppedge);

        if (nextedge == -1)
        {
            // failed to extend
            break;
        }
        TraceSegment seg;
        seg.face = begface;
        seg.side[0] = nextedge;
        seg.bary[0] = nextbary;
        seg.side[1] = begedge;
        seg.bary[1] = begbary;
        toprepend.push_back(seg);

        Eigen::Vector3d pt0 = pointFromBary(*tr.parent_, begface, seg.side[0], seg.bary[0]);
        Eigen::Vector3d pt1 = pointFromBary(*tr.parent_, begface, seg.side[1], seg.bary[1]);
        double dist = (pt1 - pt0).norm();
        extbeginning -= dist;
        actualextbeginning += dist;
    }

    // extend end
    std::vector<TraceSegment> toappend;
    int endface = tr.segs.back().face;
    int endedge = tr.segs.back().side[0];
    double endbary = tr.segs.back().bary[0];
    v0 = pointFromBary(*tr.parent_, endface, endedge, endbary);
    v1 = pointFromBary(*tr.parent_, endface, tr.segs.back().side[1], tr.segs.back().bary[1]);
    dir = v1 - v0;
    dir *= 1e3 / dir.norm();
    getNextTracePoint(*tr.parent_, endface, endedge, endbary, dir, nextedge, nextbary, oppface, oppedge);    

    while (extend > 0)
    {
        if (nextedge == -1 || oppface == -1)
        {
            // failed to extend
            break;
        }

        // project vector into barycentric coordinates
        Eigen::Matrix<double, 3, 2> B = tr.parent_->data().Bs[endface];
        Eigen::Vector2d barys = (B.transpose()*B).inverse() * B.transpose() * dir;
        int edgeid = tr.parent_->data().faceEdges(endface, nextedge);
        int side = 0;
        if (tr.parent_->data().E(edgeid, side) != endface)
            side = 1;
        Eigen::Vector2d nextbarys = tr.parent_->data().Ts.block<2, 2>(2 * edgeid, 2 * side) * barys;
        dir = tr.parent_->data().Bs[oppface] * nextbarys;

        endface = oppface;
        endedge = oppedge;
        endbary = 1.0 - nextbary;

        getNextTracePoint(*tr.parent_, endface, endedge, endbary, dir, nextedge, nextbary, oppface, oppedge);

        if (nextedge == -1)
        {
            // failed to extend
            break;
        }
        TraceSegment seg;
        seg.face = endface;
        seg.side[0] = endedge;
        seg.bary[0] = endbary;
        seg.side[1] = nextedge;
        seg.bary[1] = nextbary;
        toappend.push_back(seg);

        Eigen::Vector3d pt0 = pointFromBary(*tr.parent_, endface, seg.side[0], seg.bary[0]);
        Eigen::Vector3d pt1 = pointFromBary(*tr.parent_, endface, seg.side[1], seg.bary[1]);
        double dist = (pt1 - pt0).norm();
        extend -= dist;
        actualextend += dist;
    }

    std::reverse(toprepend.begin(), toprepend.end());
    for (auto &it : tr.segs)
        toprepend.push_back(it);
    for (auto &it : toappend)
        toprepend.push_back(it);
    tr.segs = toprepend;
}

double TraceSet::arclength(const Trace &tr) const
{
    double result = 0;
    int nsegs = tr.segs.size();
    for (int i = 0; i < nsegs; i++)
    {
        Eigen::Vector3d pt0 = pointFromBary(*tr.parent_, tr.segs[i].face, tr.segs[i].side[0], tr.segs[i].bary[0]);
        Eigen::Vector3d pt1 = pointFromBary(*tr.parent_, tr.segs[i].face, tr.segs[i].side[1], tr.segs[i].bary[1]);
        result += (pt1 - pt0).norm();
    }
    return result;
}


Eigen::Vector3d parallelTransport(const Eigen::Vector3d &v, const Eigen::Vector3d &e1, const Eigen::Vector3d &e2)
{
    Eigen::Vector3d t1 = e1 / e1.norm();
    Eigen::Vector3d t2 = e2 / e2.norm();
    Eigen::Vector3d n = t1.cross(t2);
    if (n.norm() < 1e-8)
        return v;
    n /= n.norm();
    Eigen::Vector3d p1 = n.cross(t1);
    Eigen::Vector3d p2 = n.cross(t2);
    return v.dot(n)*n + v.dot(t1)*t2 + v.dot(p1)*p2;
}

void TraceSet::findPointOnTrace(const Trace &tr, double s, int &seg, double &bary)
{
    // this should be rewritten to be O(log n) using binary search

    if (s < 0)
    {
        seg = 0;
        bary = 0;
        return;
    }

    double curs = 0;
    int nsegs = tr.segs.size();
    for (int i = 0; i < nsegs; i++)
    {
        Eigen::Vector3d pt0 = pointFromBary(*tr.parent_, tr.segs[i].face, tr.segs[i].side[0], tr.segs[i].bary[0]);
        Eigen::Vector3d pt1 = pointFromBary(*tr.parent_, tr.segs[i].face, tr.segs[i].side[1], tr.segs[i].bary[1]);
        double dist = (pt1 - pt0).norm();
        if (curs <= s && s < curs + dist)
        {
            bary = (s - curs) / dist;
            seg = i;
            return;
        }
        curs += dist;
    }
    // past end of rod
    seg = nsegs - 1;
    bary = 1.0;
}

void TraceSet::sampleTrace(const Trace &tr, double start, double end, int nsegs, RationalizedTrace &rattrace, std::vector<double> &samples)
{
    samples.clear();
    rattrace.normals.resize(nsegs, 3);
    rattrace.pts.resize(nsegs + 1, 3);    
    for (int i = 0; i < nsegs+1; i++)
    {
        double s = start + (end - start)*double(i) / double(nsegs);
        samples.push_back(s);
        int seg;
        double bary;
        findPointOnTrace(tr, s, seg, bary);
        Eigen::Vector3d pt0 = pointFromBary(*tr.parent_, tr.segs[seg].face, tr.segs[seg].side[0], tr.segs[seg].bary[0]);
        Eigen::Vector3d pt1 = pointFromBary(*tr.parent_, tr.segs[seg].face, tr.segs[seg].side[1], tr.segs[seg].bary[1]);
        Eigen::Vector3d interp = (1.0 - bary)*pt0 + bary*pt1;
        rattrace.pts.row(i) = interp.transpose();
    }

    for (int i = 0; i < nsegs; i++)
    {
        double s = start + (end - start)*double(2*i+1) / (2.0*nsegs);
        int seg;
        double bary;
        findPointOnTrace(tr, s, seg, bary);
        Eigen::Vector3d oldpt0 = pointFromBary(*tr.parent_, tr.segs[seg].face, tr.segs[seg].side[0], tr.segs[seg].bary[0]);
        Eigen::Vector3d oldpt1 = pointFromBary(*tr.parent_, tr.segs[seg].face, tr.segs[seg].side[1], tr.segs[seg].bary[1]);
        Eigen::Vector3d normal = tr.parent_->faceNormal(tr.segs[seg].face);

        Eigen::Vector3d newpt0 = rattrace.pts.row(i);
        Eigen::Vector3d newpt1 = rattrace.pts.row(i+1);
        Eigen::Vector3d newnormal = parallelTransport(normal, oldpt1 - oldpt0, newpt1 - newpt0);
        rattrace.normals.row(i) = newnormal.transpose();
    }

}

void TraceSet::rationalizeTraces(double maxcurvature, double extenddist, double seglen, double minlen)
{
    rattraces_.clear();
    collisions_.clear();

    std::vector<Trace> cleanedtraces;
    for (int i = 0; i < traces_.size(); i++)
    {
        std::set<int> badverts;
        findCurvedVerts(traces_[i], maxcurvature, badverts);        
        std::vector<Trace> splittr;
        splitTrace(traces_[i], badverts, splittr);
        for (auto &it : splittr)
        {
            if (minlen < arclength(it))
                cleanedtraces.push_back(it);
        }
    }
    
    // extend traces if desired
    std::vector<double> starts;
    std::vector<double> ends;
    for (auto &it : cleanedtraces)
    {
        double actualextbeginning;
        double actualextend;
        extendTrace(it, extenddist, extenddist, actualextbeginning, actualextend);
        double len = arclength(it);
        double start = actualextbeginning - extenddist;
        start = std::max(start, 0.0);
        start = std::min(start, len);
        starts.push_back(start);
        double end = len + extenddist - actualextend;
        end = std::max(end, 0.0);
        end = std::min(end, len);
        ends.push_back(end);
    }

    // find collisions
    int ntraces = cleanedtraces.size();

    std::map<std::pair<int, int>, std::vector<TraceCollision> > cols;
    for (int i = 0; i < ntraces; i++)
    {
        for (int j = i; j < ntraces; j++)
        {
            bool self = i == j;
            std::vector<TraceCollision> paircols;
            computeIntersections(cleanedtraces[i], cleanedtraces[j], self, paircols);
            cols[std::pair<int, int>(i, j)] = paircols;
        }
    }

    // convert collisions to pairs of arclength values
    std::vector<std::vector<double> > svals;
    svals.resize(ntraces);
    for (int i = 0; i < ntraces; i++)
    {
        svals[i].push_back(0);
        double s = 0;
        int nsegs = cleanedtraces[i].segs.size();
        for (int j = 0; j < nsegs; j++)
        {
            Eigen::Vector3d pt0 = pointFromBary(*cleanedtraces[i].parent_, cleanedtraces[i].segs[j].face, cleanedtraces[i].segs[j].side[0], cleanedtraces[i].segs[j].bary[0]);
            Eigen::Vector3d pt1 = pointFromBary(*cleanedtraces[i].parent_, cleanedtraces[i].segs[j].face, cleanedtraces[i].segs[j].side[1], cleanedtraces[i].segs[j].bary[1]);
            double dist = (pt1 - pt0).norm();
            s += dist;
            svals[i].push_back(s);
        }
    }
    struct ArcCollision
    {
        int rod1, rod2;
        double s1, s2;
    };
    std::vector<ArcCollision> arccols;
    for (auto &it : cols)
    {
        int ncols = it.second.size();
        for (int i = 0; i < ncols; i++)
        {
            ArcCollision ac;
            ac.rod1 = it.first.first;
            ac.rod2 = it.first.second;
            ac.s1 = (1.0 - it.second[i].bary1)*svals[ac.rod1][it.second[i].seg1] + it.second[i].bary1 * svals[ac.rod1][it.second[i].seg1 + 1];
            ac.s2 = (1.0 - it.second[i].bary2)*svals[ac.rod2][it.second[i].seg2] + it.second[i].bary2 * svals[ac.rod2][it.second[i].seg2 + 1];
            arccols.push_back(ac);
        }
    }

    // sample traces into rod segments
    std::vector<std::vector<double> > samples;
    samples.resize(ntraces);
    for (int i = 0; i < ntraces; i++)
    {
        int nsegs = 1 + int( (ends[i] - starts[i]) / seglen);
        RationalizedTrace rat;
        sampleTrace(cleanedtraces[i], starts[i], ends[i], nsegs, rat, samples[i]);
        rattraces_.push_back(rat);
    }

    // compute collisions on sampled rod segments
    // should be rewritten to be O(m log n) with binary search
    for (auto &it : arccols)
    {
        int seg1 = -1;
        double bary1 = 0;
        double curs = samples[it.rod1][0];
        if (it.s1 < curs)
            continue;
        for (int i = 0; i < samples[it.rod1].size() - 1; i++)
        {
            double nexts = samples[it.rod1][i + 1];
            if (curs <= it.s1 && it.s1 < nexts)
            {
                seg1 = i;
                bary1 = (it.s1 - curs) / (nexts - curs);
                break;
            }
            curs = nexts;
        }
        if (seg1 == -1)
            continue;

        int seg2 = -1;
        double bary2 = 0;
        curs = samples[it.rod2][0];
        if (it.s2 < curs)
            continue;
        for (int i = 0; i < samples[it.rod2].size() - 1; i++)
        {
            double nexts = samples[it.rod2][i + 1];
            if (curs <= it.s2 && it.s2 < nexts)
            {
                seg2 = i;
                bary2 = (it.s2 - curs) / (nexts - curs);
                break;
            }
            curs = nexts;
        }
        if (seg2 == -1)
            continue;

        Collision col;
        col.rod1 = it.rod1;
        col.seg1 = seg1;
        col.bary1 = bary1;
        col.rod2 = it.rod2;
        col.seg2 = seg2;
        col.bary2 = bary2;
        collisions_.push_back(col);
    }

    // Project along a geodesic curve at the end point of each rod that ends at a singularity.  
    /*
    Eigen::VectorXd sqrD;
    Eigen::VectorXi nearFace;
    Eigen::MatrixXd C;
    Eigen::MatrixXd point = Eigen::MatrixXd::Zero(1,3); 
    igl::AABB<Eigen::MatrixXd,3> tree; 

    tree.init(wv.fs->data().V,wv.fs->data().F);


    // need to deal with boundary case
    int stepstoextend = 0;
    int backoff = 0; // The very ends of curves seem generally bad.  Taking a few steps to back off

    for (int i = 0; i < curves.size(); i++)
    { 
        int len = curves[i].rows();
        int cols = curves[i].cols();
        curves[i].conservativeResize(len + stepstoextend - backoff, cols);
        normals[i].conservativeResize(len + stepstoextend - backoff, cols);

        Eigen::Vector3d endpoint = curves[i].row(len-1 - backoff);
        Eigen::Vector3d prevpoint = curves[i].row(len-2 - backoff);
        Eigen::Vector3d curr_dir = endpoint - prevpoint;
        // hack-y way of figuring out the edge of the end-point 
        point.row(0) = prevpoint + curr_dir * .01;

        curr_dir *= 1000.0; // Assumes mesh is roughly deluanay

        tree.squared_distance(wv.fs->data().V,wv.fs->data().F,point,sqrD,nearFace,C);
        int curr_face_id = nearFace(0);

        TracePoint tp;
        startTraceFromPoint(wv, curr_face_id, point.row(0), curr_dir, tp);
        //        std::cout << point.row(0) - prevpoint << " should be same \n";
        int curr_edge_id = tp.edge_id;
        getNextTracePoint(wv, curr_face_id, tp.edge_id, prevpoint, tp.op_v_id, curr_dir, tp); 
        //        std::cout << tp.point - endpoint << " endpoints should be same \n";

        for (int j = 0; j < stepstoextend; j++)
        {
            curr_dir = mapVectorToAdjacentFace(wv.fs->data().F, wv.fs->data().V, wv.fs->data().edgeVerts,
                tp.edge_id, curr_face_id, tp.face_id, curr_dir);
            curr_face_id = tp.face_id;  
            curves[i].row(len + j - backoff)  = tp.point;
            normals[i].row(len + j - backoff) = tp.n;
            std::cout << normals[i] << "\n";

            getNextTracePoint(wv, curr_face_id, tp.edge_id, tp.point, tp.op_v_id, curr_dir, tp); 
        }

        Eigen::MatrixXd revcurve = curves[i];
        Eigen::MatrixXd revnorm = normals[i];
        int c_len = curves[i].rows();
        for (int r = 0; r < c_len; r++)
        {
            revcurve.row(r) = curves[i].row(c_len - r - 1);
            revnorm.row(r) = normals[i].row(c_len - r - 1);
        }
        curves[i] = revcurve;
        normals[i] = revnorm;
    }   

    // Cut strips at areas of high curvature.

    std::vector<Eigen::MatrixXd> splitcurves;
    std::vector<Eigen::MatrixXd> splitnormals;

    double maxcurvature = 2000; // angle - currently disabled
    double minrodlen = 2;

    for(int i=0; i<curves.size(); i++)
    {
        if(curves[i].rows() < 3)
            continue;
        std::vector<Eigen::Vector3d> curpts;
        std::vector<Eigen::Vector3d> curnormals;
        curpts.push_back(curves[i].row(0));
        curpts.push_back(curves[i].row(1));
        curnormals.push_back(normals[i].row(0));
        curnormals.push_back(normals[i].row(1));
        Eigen::Vector3d prevedge = curves[i].row(1) - curves[i].row(0);
        for(int j=2; j<curves[i].rows(); j++)
        {
            Eigen::Vector3d nextedge = curves[i].row(j).transpose() - curpts.back();
            Eigen::Vector3d prevproj = prevedge - prevedge.dot(normals[i].row(j)) * normals[i].row(j).transpose();
            if(fabs(angle(prevproj, nextedge, normals[i].row(j)))/(prevedge.norm() + nextedge.norm()) > maxcurvature)
            {
                // cut
                if(curpts.size() > minrodlen)
                {
                    Eigen::MatrixXd newcurve(curpts.size(), 3);
                    Eigen::MatrixXd newnormal(curpts.size(), 3);
                    for(int j=0; j<curpts.size(); j++)
                    {
                        newcurve.row(j) = curpts[j].transpose();
                        newnormal.row(j) = curnormals[j].transpose();
                    }
                    splitcurves.push_back(newcurve);
                    splitnormals.push_back(newnormal);
                }
                curpts.clear();
                curnormals.clear();
                curpts.push_back(curves[i].row(j-1));
                curpts.push_back(curves[i].row(j));
                curnormals.push_back(normals[i].row(j-1));
                curnormals.push_back(normals[i].row(j));
            }
            else
            {
                curpts.push_back(curves[i].row(j));
                curnormals.push_back(normals[i].row(j));
            }
            prevedge = nextedge;
        }
        if(curpts.size() > minrodlen)
        {
            Eigen::MatrixXd newcurve(curpts.size(), 3);
            Eigen::MatrixXd newnormal(curpts.size(), 3);
            for(int j=0; j<curpts.size(); j++)
            {
                newcurve.row(j) = curpts[j].transpose();
                newnormal.row(j) = curnormals[j].transpose();
            }
            splitcurves.push_back(newcurve);
            splitnormals.push_back(newnormal);
        }             
    }









    // Find collisions between rods
    std::vector<Collision> collisions;
    for (int i = 0; i < splitcurves.size(); i++)
    {
        for (int j = i; j < splitcurves.size(); j++)
        {
            computeIntersections(i, j, splitcurves, collisions);
        }
    }

    // Decimate and log rods
    std::vector<Eigen::VectorXd> desc_maps;
    std::vector< std::vector<Eigen::Vector3d> > desc_curves; //eeew
    std::vector< std::vector<Eigen::Vector3d> > desc_normals;

    double decimation_factor = 1.; // make this bigger and add intelligent subdivision

    for (int curveId = 0; curveId < splitcurves.size(); curveId++)
    {
        Eigen::MatrixXd curve = splitcurves[curveId];
        Eigen::MatrixXd curveNormals = splitnormals[curveId];

        double max_length = 0.;
        for (int i = 0; i < curve.rows() - 1; i++)
        {
            double seg_length = (curve.row(i) - curve.row(i + 1)).norm();
            if (seg_length > max_length)
            {
                max_length = seg_length;
                //        max_length = 0.;
            }
        }

        // Decimate 
        std::vector<Eigen::Vector3d> cnew;
        std::vector<Eigen::Vector3d> nnew;
        cnew.push_back(curve.row(0));
        int seg_counter = 0;
        Eigen::VectorXd desc_mapping = Eigen::VectorXd::Zero(curve.rows());
        Eigen::Vector3d prev_point = cnew.back();
        for (int i = 1; i < curve.rows(); i++)
        {
            Eigen::Vector3d curr_point = curve.row(i);
            double seg_length = (prev_point - curr_point).norm();
            desc_mapping(i-1) = seg_counter;
            if (seg_length > max_length / decimation_factor)
            {
                seg_counter++;
                cnew.push_back(curve.row(i));
                Eigen::Vector3d currEdge = curve.row(i) - curve.row(i - 1);
                Eigen::Vector3d targEdge = cnew[seg_counter] - cnew[seg_counter - 1];
                nnew.push_back(parallelTransport(curveNormals.row(i), currEdge, targEdge));
                prev_point = cnew.back();
            }
        }

        desc_maps.push_back(desc_mapping);
        desc_curves.push_back(cnew);
        desc_normals.push_back(nnew);
    }


    std::vector<Collision> desc_collisions;
    for (int i = 0; i < collisions.size(); i++)
    {
        Collision col = collisions[i];
        std::vector<Eigen::Vector3d> &c1 = desc_curves[col.rod1];
        std::vector<Eigen::Vector3d> &c2 = desc_curves[col.rod2];
        int idx1 = desc_maps[col.rod1](col.seg1);
        int idx2 = desc_maps[col.rod2](col.seg2);
        if(idx1 < c1.size() - 1 && idx2 < c2.size() - 1)
        {
            Collision newcol(col.rod1, col.rod2, idx1, idx2);
            desc_collisions.push_back(newcol);
        }
    }

    // Write Header 
    myfile << desc_curves.size() << std::endl;;
    myfile << desc_collisions.size() << std::endl;;
    //  myfile << 0 << std::endl;;
    myfile << "0.001"  << std::endl;;
    myfile << "1e+08"  << std::endl;;
    myfile << "1"  << std::endl  << std::endl  << std::endl;;


    for (int curveId = 0; curveId < desc_curves.size(); curveId++)
    {
        std::vector<Eigen::Vector3d> &cnew = desc_curves[curveId];
        std::vector<Eigen::Vector3d> &nnew = desc_normals[curveId];
        myfile << cnew.size() << "\n";
        myfile << "0\n";
        for (int i = 0; i < cnew.size(); i++)
        {
            myfile << cnew[i](0) << " " << cnew[i](1) << " " << cnew[i](2) << " ";
        }
        myfile << "\n";

        for (int i = 0; i < nnew.size(); i++)
        {
            myfile << nnew[i](0) << " " << nnew[i](1) << " " << nnew[i](2) << " ";
        }
        myfile << "\n";

        for (int i = 0; i < cnew.size()-1; i++)
        {
            myfile << " 0.02";
        }
        myfile << "\n";

    }
    for (int i = 0; i < desc_collisions.size(); i++)
    {
        Collision col = desc_collisions[i];
        std::vector<Eigen::Vector3d> c1 = desc_curves[col.rod1];
        std::vector<Eigen::Vector3d> c2 = desc_curves[col.rod2];
        int idx1 = col.seg1;
        int idx2 = col.seg2;

        double p0bary, p1bary, q0bary, q1bary;

        if(idx1 >= c1.size()-1) exit(-1);

        Eigen::Vector3d dist = Distance::edgeEdgeDistance(c1[idx1],
            c1[idx1 + 1],
            c2[idx2],
            c2[idx2 + 1],
            p0bary, p1bary, q0bary, q1bary);

        myfile << col.rod1 << " " << col.rod2 << " " << idx1 << " " << idx2
            << " " << p1bary << " " << q1bary << " " << 1000. << "\n";

    }


    myfile.close();*/
}

void TraceSet::computeIntersections(const Trace &tr1, const Trace &tr2, bool selfcollision, 
    std::vector<TraceCollision> &collisions) const
{
    collisions.clear();
    int nsegs1 = tr1.segs.size();
    int nsegs2 = tr2.segs.size();

    for (int i = 0; i < nsegs1; i++)
    {
        for (int j = 0; j < nsegs2; j++)
        {
            if (selfcollision && (i - j) < 2) { continue; }
            else
            {
                Eigen::Vector3d p0 = pointFromBary(*tr1.parent_, tr1.segs[i].face, tr1.segs[i].side[0], tr1.segs[i].bary[0]);
                Eigen::Vector3d p1 = pointFromBary(*tr1.parent_, tr1.segs[i].face, tr1.segs[i].side[1], tr1.segs[i].bary[1]);
                Eigen::Vector3d q0 = pointFromBary(*tr2.parent_, tr2.segs[j].face, tr2.segs[j].side[0], tr2.segs[j].bary[0]);
                Eigen::Vector3d q1 = pointFromBary(*tr2.parent_, tr2.segs[j].face, tr2.segs[j].side[1], tr2.segs[j].bary[1]);
                double p0bary, p1bary, q0bary, q1bary;
                Eigen::Vector3d dist = Distance::edgeEdgeDistance(p0, p1, q0, q1,
                    p0bary, p1bary, q0bary, q1bary);
                if (dist.norm() < 1e-6 && p0bary != 0 && p0bary != 1.0 && q0bary != 0 && q0bary != 1.0)
                {
                    TraceCollision tc;
                    tc.seg1 = i;
                    tc.seg2 = j;
                    tc.bary1 = p1bary;
                    tc.bary2 = q1bary;
                    collisions.push_back(tc);
                }
            }
        }
    }
}

void TraceSet::collisionPoint(int collision, Eigen::Vector3d &pt0, Eigen::Vector3d &pt1) const
{
    const Collision &col = collisions_[collision];
    pt0 = (1.0 - col.bary1) * rattraces_[col.rod1].pts.row(col.seg1).transpose() + col.bary1 * rattraces_[col.rod1].pts.row(col.seg1 + 1).transpose();
    pt1 = (1.0 - col.bary2) * rattraces_[col.rod2].pts.row(col.seg2).transpose() + col.bary2 * rattraces_[col.rod2].pts.row(col.seg2 + 1).transpose();
}