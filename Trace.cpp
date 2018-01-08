#include <igl/avg_edge_length.h>
#include <igl/viewer/Viewer.h>

#include "Trace.h"
#include "VectorUtils.h"
#include "Distance.h"

double energy_OP = 0.;
//void updateView(const MeshData &curMesh, igl::viewer::Viewer &viewer){}


Eigen::MatrixXd colorField;


/*
// This: Recieves a point on an edge of a triangle
//       Returns the vertex id opposing that edge in the triangle
int getCurrEdge(const MeshData &md, const Eigen::Vector3d prev_point, int faceId)
{
    for (int j = 0; j < 3; j++)
    {
    Eigen::VectorXi e = md.E.row(md.F_edges(faceId, j));
    Eigen::Vector3d e_test =  md.V.row( e(0) ) - md.V.row( e(1) );
    e_test.normalize();
        Eigen::Vector3d point_test = md.V.row( e(0) );
        point_test = ( prev_point - point_test );
    point_test.normalize();
    if ( std::abs( e_test.dot( point_test )) > .9 )
    {
        return md.F_edges(faceId, j);
    } 
    }
    assert(false);
    return -1;
}


// This: Recieves a point on an edge of a triangle
//       Returns the vertex id opposing that edge in the triangle
int getOpVId(const MeshData &md, const Eigen::Vector3d prev_point, int faceId)
{
    for (int j = 0; j < 3; j++)
    {
    Eigen::VectorXi e = md.E.row(md.F_edges(faceId, j));
    Eigen::Vector3d e_test =  md.V.row( e(0) ) - md.V.row( e(1) );
    e_test.normalize();
        Eigen::Vector3d point_test = md.V.row( e(0) );
        point_test = ( prev_point - point_test );
    point_test.normalize();
    if ( e_test.dot( point_test ) > .99 )
    {
        Eigen::VectorXi e_next = md.E.row(md.F_edges(faceId, (j + 1) % 3 ));
        if (e_next(0) == e(0) || e_next(0) == e(1) )
        {
                return e_next(1);
        }
        else 
        {
                return e_next(0);
        }
    } 
    }
    return -1;
}

int getOpVIdFromEdge(const MeshData &md, int curr_edge, int faceId)
{
    Eigen::VectorXi e = md.E.row(curr_edge);   
    for (int i = 0; i < 3; i++) 
    {
        if ( md.F(faceId, i) != e(0) && md.F(faceId, i) != e(1) )
    {
            return md.F(faceId, i);
    }
    }
    assert(false);
    return -1;
}

void computeSelfIntersections(const std::vector<Eigen::Vector3d> &curve, int idx, std::vector<Collision> &collisions)
{
    double minSegmentLength = (curve[0] - curve[1]).norm();
    for (int i = 0; i < curve.size() - 1; i++)
    {
    double segLength = (curve[i] - curve[i + 1]).norm();
        if ( minSegmentLength > segLength && segLength > 0.)
    {
            minSegmentLength = segLength;
    }
    }
    minSegmentLength = minSegmentLength / 2.;
    std::cout << minSegmentLength << "\n";
    for (int i = 0; i < curve.size() - 3; i++) 
    {
        for (int j = i + 3; j < curve.size(); j++)
    {
         
        double p0bary, p1bary, q0bary, q1bary;
        Eigen::Vector3d dist = Distance::edgeEdgeDistance(curve[i], 
                                                      curve[i + 1],
                                                              curve[j],
                                  curve[j + 1],
                                  p0bary, p1bary, q0bary, q1bary);
        if ( dist.norm() < minSegmentLength )
        {
             //   Collision* c = new Collision(idx, idx, i, j);
        collisions.push_back( Collision(idx, idx, i, j)  );
        std::cout << idx << " " << i << " " << j << " "  << dist.norm() << "\n";
        }
    } 
    }
}


void traceCurve(const Weave &md,
                const Eigen::Vector3d dir, int faceId, 
                std::vector<Eigen::Vector3d> &curve,
                std::vector<Eigen::Vector3d> &normal)
{
    curve.clear();
    normal.clear();
    // This is hacky, make it possible to shoot from point in bary-centric coordinates
    curve.push_back( .5 * md.V.row(md.F(faceId, 0)) + .5 * md.V.row(md.F(faceId, 1)) ); 
    // check that vector is pointing in.  
    // TODO
    // find next edge
    int curr_face_id = faceId;
    Eigen::Vector3d curr_dir = -dir;
    
    int steps = 1000;

    int curr_edge_id = getCurrEdge(md, curve.back(), curr_face_id);
    for (int i = 0; i < steps; i++)
    {
        Eigen::Vector3d prev_point = curve.back();
    int op_v_id = getOpVIdFromEdge(md, curr_edge_id, curr_face_id);
        
        Eigen::Vector3d op_vertex = md.V.row(op_v_id);    
        Eigen::Vector3d split = op_vertex - prev_point;
    
    double split_len = split.norm();
    split = split / split_len;
    Eigen::Vector3d n = faceNormal(md.F, md.V, curr_face_id);
    Eigen::Vector3d perp = split.cross(n);
        normal.push_back(-n);

    int op_edge_id = 0;
        for (int j = 0; j < 3; j++)
    {
        Eigen::VectorXi e = md.E.row(md.F_edges(curr_face_id, j));
        
        if ( e(0) == op_v_id || e(1) == op_v_id )
            {
                Eigen::VectorXd e_test = (md.V.row(e(0)) + md.V.row(e(1))) * .5;
        e_test -= prev_point;
                if ( e_test.dot(perp) * curr_dir.dot(perp) > 0.) 
                {
                    op_edge_id = j;
            break;
                }               
            }
        }
        
    // Find intersection point.
        int next_edge_id = md.F_edges(curr_face_id, op_edge_id);
        Eigen::Vector3d op_v1 = md.V.row(md.E(next_edge_id, 0));
        Eigen::Vector3d op_v2 = md.V.row(md.E(next_edge_id, 1));

    // This 100 is hacky, figure out the appropriate scaling...
    double p0bary, p1bary, q0bary, q1bary;
    Eigen::Vector3d dist = Distance::edgeEdgeDistance(prev_point - curr_dir * split_len * 100, 
                                                  prev_point + curr_dir * split_len * 100,
                                                          op_v1, op_v2,
                              p0bary, p1bary, q0bary, q1bary);
    
    curve.push_back( op_v1 * q0bary + op_v2 * q1bary);
        int next_face_id = md.E(next_edge_id, 2);
        if ( next_face_id == curr_face_id )
    {
            next_face_id = md.E(next_edge_id, 3);
    }
        if ( next_face_id == -1) { break; }    
    
    curr_dir = mapVectorToAdjacentFace(md.F, md.V, md.E, 
              next_edge_id, curr_face_id, next_face_id, curr_dir);
        curr_face_id = next_face_id;
    curr_edge_id = next_edge_id;
    }
    return; 
}
*/
