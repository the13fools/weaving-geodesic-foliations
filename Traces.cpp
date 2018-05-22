#include "Traces.h"
#include "FieldSurface.h"
#include "Distance.h"
#include <algorithm>
#include <fstream>

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
}
    

void TraceSet::exportRodFile(const char*filename)
{
    std::ofstream myfile(filename);
    // Write Header 
    myfile << -217 << std::endl;
    myfile << 1 << std::endl;
    myfile << rattraces_.size() << std::endl;;
    myfile << collisions_.size() << std::endl;;
    myfile << "0.001"  << std::endl;;
    myfile << "1e+08"  << std::endl;;
    myfile << "1"  << std::endl  << std::endl  << std::endl;;

    int color=0;
    for (auto &it : rattraces_)
    {
        int nverts = it.pts.rows();
        int nsegs = it.normals.rows();        
        myfile << nverts << std::endl;
        myfile << 0 << std::endl;
        myfile << 1 << std::endl;
        myfile << color << std::endl;
        color = (color+1)%7;
        
        for (int i = 0; i < nverts; i++)
        {
            myfile << it.pts(i,0) << " " << it.pts(i,1) << " " << it.pts(i,2) << " ";
        }
        myfile << std::endl;

        for (int i = 0; i < nsegs; i++)
        {
            myfile << it.normals(i,0) << " " << it.normals(i,1) << " " << it.normals(i,2) << " ";
        }
        myfile << std::endl;

        for (int i = 0; i < nsegs; i++)
        {
            myfile << "0 ";
        }
        myfile << std::endl;
        
        for(int i=0; i<nsegs; i++)
        {
            myfile << "0.02 " << std::endl;
        }
        myfile << std::endl;
    }
    
    int ncollisions = collisions_.size();
    for (int i = 0; i < ncollisions; i++)
    {
        Collision &col = collisions_[i];
        myfile << col.rod1 << std::endl;
        myfile << col.rod2 << std::endl;
        myfile << col.seg1 << std::endl;
        myfile << col.seg2 << std::endl;
        myfile << col.bary1 << std::endl;
        myfile << col.bary2 << std::endl;
        myfile << "1000." << std::endl;

    }
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
