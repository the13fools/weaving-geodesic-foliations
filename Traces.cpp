#include "Traces.h"
#include "FieldSurface.h"
#include "Distance.h"

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
    int &opp_edge_id)
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

Eigen::Vector3d TraceSet::pointFromBary(const FieldSurface &parent, int faceId, int faceEdge, double bary)
{
    Eigen::Vector3d v0 = parent.data().V.row(parent.data().F(faceId, (faceEdge + 1) % 3)).transpose();
    Eigen::Vector3d v1 = parent.data().V.row(parent.data().F(faceId, (faceEdge + 2) % 3)).transpose();
    return (1.0 - bary)*v0 + bary*v1;
}