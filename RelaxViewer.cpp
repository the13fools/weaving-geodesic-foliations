#include <igl/avg_edge_length.h>
#include <igl/viewer/Viewer.h>

#include "RelaxViewer.h"
#include "VectorUtils.h"
#include "Distance.h"

double energy_OP = 0.;
//void updateView(const MeshData &curMesh, igl::viewer::Viewer &viewer){}


Eigen::MatrixXd colorField;

shading_enum shading_enum_state = INIT_MAGNITUDE;


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
	if ( e_test.dot( point_test ) > .99 )
	{
	    return md.F_edges(faceId, j);
	} 
    }
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

int getOpVId_PrevEdge(const MeshData &md, int curr_edge, int faceId)
{
    Eigen::VectorXi e = md.E.row(curr_edge);   
    for (int i = 0; i < 3; i++) 
    {
        if ( md.F(faceId, i) != e(0) && md.F(faceId, i) != e(1) )
	{
            return md.F(faceId, i);
	}
    }

}


void traceCurve(const MeshData &md, const Eigen::Vector3d dir, int faceId, VisualizationState &vs)
{
    vs.curve.clear();
    vs.curve.push_back(.5 * md.V.row(md.F(faceId, 0)) + .5 * md.V.row(md.F(faceId, 1)) ); 
    // check that vector is pointing in.  
    // TODO
    // find next edge
    int curr_face_id = faceId;
    Eigen::Vector3d curr_dir = dir;
    
    std::cout << curr_dir << "\n";
    int steps = 20;
    vs.curve_dirs = Eigen::MatrixXd::Zero(steps, 3);

    int curr_edge_id  = getCurrEdge(md, vs.curve.back(), curr_face_id);
    for (int i = 0; i < steps; i++)
    {
	vs.curve_dirs.row(i) = dir;

        Eigen::Vector3d prev_point = vs.curve.back();
	int op_v_id = getOpVId_PrevEdge(md, curr_edge_id, curr_face_id);
        
        Eigen::Vector3d op_vertex = md.V.row(op_v_id);	
        Eigen::Vector3d split = prev_point - op_vertex;
        double split_len = split.norm();
	split = split / split_len;
	Eigen::Vector3d n = faceNormal(md.F, md.V, curr_face_id);
	Eigen::Vector3d perp = split.cross(n);
        
    //    std::cout << "Iteration i: " << i << " "  <<  split_len << " split len" << "\n";

	int op_edge_id = 0;
        for (int j = 0; j < 3; j++)
	{
	    Eigen::VectorXi e = md.E.row(md.F_edges(curr_face_id, j));
//	    std::cout << e << " " << j  << " edge \n";
	    if ( e(0) == op_v_id || e(1) == op_v_id )
            {
                Eigen::VectorXd e_test = md.V.row(e(0)) - md.V.row(e(1));
                if ( e_test.dot(perp) * dir.dot(perp) > 0.) 
                {
                    op_edge_id = j;
		    break;
                }               
            }
        }
    //    std::cout << op_edge_id << "op edge id"  << std::endl << std::endl; 
        // Find intersection point.
        int next_edge_id = md.F_edges(curr_face_id, op_edge_id);
    	Eigen::Vector3d op_v1 = md.V.row(md.E(next_edge_id, 0));
        Eigen::Vector3d op_v2 = md.V.row(md.E(next_edge_id, 1));
    //    std::cout << op_v1 << "op_v1"<< std::endl << op_v2 << "<-- op_v2"<<  std::endl;    

	double p0bary, p1bary, q0bary, q1bary;
	Eigen::Vector3d dist = Distance::edgeEdgeDistance(prev_point, 
		                                          prev_point + dir * split_len,
                                                          op_v1, op_v2,
							  p0bary, p1bary, q0bary, q1bary);
  //      std::cout << dist.norm() << " <- dist norm \n  " << prev_point << " <-- prev point " << std::endl;
	
//	std::cout <<  p0bary << " " <<   p1bary << " " << q0bary << " " << q1bary << " barys\n";
	vs.curve.push_back( op_v1 * q0bary + op_v2 * q1bary);
        int next_face_id = md.E(next_edge_id, 2);
        if ( next_face_id == curr_face_id )
	{
            next_face_id = md.E(next_edge_id, 3);
	}
        if ( next_face_id == -1) { break; }	
	
    std::cout << curr_dir << "\n";

  //      std::cout << next_edge_id << " " << curr_face_id << " " << next_face_id << "edge face face\n";
//	curr_dir = mapVectorToAdjacentFace(md.F, md.V, md.E, 
//		      next_edge_id, curr_face_id, next_face_id, curr_dir);
        curr_face_id = next_face_id;
	curr_edge_id = next_edge_id;
    }
    return; 
}



void updateView(const MeshData *curMesh, igl::viewer::Viewer *viewer)
{
    int nFaces = curMesh->F.rows();

    // Set mesh colors and log operator state
    Eigen::VectorXd Z(nFaces);
    Eigen::Vector3d testDir(1,0,0);
    energy_OP = 0.;

    double max_error = 0.;
    for (int i = 0; i < nFaces; i++)
    {
        switch (shading_enum_state)
        {
            case OP_ERROR:
     //           Z(i) = log(del_W_F.row(i).norm() + .000005);
    //            Z(i) = Op_Grad.row(i).norm();
//              std::cout << Z(i) << "\n";
  //              Z(i) = (Op_Grad.row(i) + Op_Grad2.row(i) - Op_Grad_fd.row(i)).norm() / ( Op_Grad.row(i) + Op_Grad2.row(i) ).norm() * 2;
                if (Z(i) > max_error)
                {
                    std::cout << i << " " << Z(i) << "\n";
                    max_error = Z(i);
                }
                break;
            case INIT_DIRECTION:
                Z(i) = (curMesh->optVars.W_opt-curMesh->v0).row(i).normalized().dot(testDir);
//              Z(i) = (Op_Grad).row(i).normalized().dot(testDir) + .000005;
//              Z(i) = (Op_Grad).row(i).norm();
                break;
            case INIT_MAGNITUDE:
              //  Z(i) = log( (curMesh->optVars.W_opt).row(i).squaredNorm() );
             //   Z(i) = log( (curMesh->optVars.W_opt-curMesh->v0).row(i).squaredNorm() );
                Z(i) = curMesh->vs.energy(i);
//                std::cout << curMesh->vs.energy;
//              Eigen::Vector3d test_vec(-Op_Grad(i,1), Op_Grad(i,0), 0);
//              Z(i) = (Op_Grad_fd).row(i).normalized()
//                         .dot(test_vec.normalized()) + .000005;
        //        std::cout << Z(i) << "\n";
//              Z(1) = 1;
//                Z(i) = (Op_Grad2).row(i).norm();
                break;
        }

  //      energy_OP += del_W_F.row(i).squaredNorm();
    }
//    std::cout << energy_OP << " Operator Energy\n";

    // Average edge length for sizing
    const double avg = igl::avg_edge_length(curMesh->V,curMesh->F);
    colorField.resize(nFaces, 3);

    //  igl::jet(Z,true,colorField);

    igl::ColorMapType viz_color = igl::COLOR_MAP_TYPE_JET;

    switch (shading_enum_state)
    {
        case OP_ERROR:
            igl::colormap(viz_color,Z, true, colorField);
            break;
        case INIT_DIRECTION:
            igl::colormap(viz_color,Z, true, colorField);
            break;
        case INIT_MAGNITUDE:
            igl::colormap(viz_color,Z, true, colorField);
            break; // MAGMA, JET
    }


    // Plot the mesh
    viewer->data.clear();
    viewer->data.set_mesh(curMesh->V, curMesh->F);
    viewer->data.set_face_based(true);

    viewer->data.set_colors(colorField);

    Eigen::MatrixXd eps = Eigen::MatrixXd::Constant(nFaces,3,.001);

    const Eigen::RowVector3d red(0.9,.1,.1),green(0.1,0.9,0.2),blue(0.1,0.2,0.8);
  //  viewer->data.add_edges(centroids_f  + del_w_f*avg/2, centroids_f, blue);

//    viewer->data.add_edges(centroids_F  + (W - W_init)*avg/2*operator_scale, centroids_F, red);
//    viewer->data.add_edges(centroids_F  + (Op_Grad*operator_scale - Op_Grad_fd)*avg/2, centroids_F, red);
   

  //  Eigen::MatrixXd curveStarts = Eigen::MatrixXd::Zero(curMesh->vs.curve.size() - 1, 3);
    Eigen::MatrixXd curveStarts = Eigen::MatrixXd::Zero(19, 3);
    Eigen::MatrixXd curveStarts_all = Eigen::MatrixXd::Zero(curMesh->vs.curve.size(), 3);
    Eigen::MatrixXd curveEnds = Eigen::MatrixXd::Zero(curMesh->vs.curve.size() - 1, 3);
    for(std::vector<Eigen::VectorXd>::size_type i = 0; i != curMesh->vs.curve.size() - 1; i++) {
      //  viewer->data.add_edges( curMesh->vs.curve[i], curMesh->vs.curve[i + 1].transpose(), red);
	curveStarts.row(i) = curMesh->vs.curve[i];
        curveStarts_all.row(i) = curMesh->vs.curve[i];
	curveEnds.row(i)  =  curMesh->vs.curve[i + 1];
//        std::cout << curMesh->vs.curve[i] << std::endl;
    }

    viewer->data.add_edges( curveStarts, curveEnds, red);
    viewer->data.add_points( curveStarts ,green);
    
    viewer->data.add_edges( curveStarts, curveStarts +  curMesh->vs.curve_dirs, blue);
    std::cout << curMesh->vs.curve_dirs << std::endl << std::endl;
    std::cout << curveStarts << std::endl <<  std::endl;
    std::cout << curveStarts + curMesh->vs.curve_dirs << std::endl;

    if (curMesh->vs.normFaceVectors)
    {
        Eigen::MatrixXd field;
        field.resize(nFaces, 3);
        for (int i = 0; i < nFaces; i++)
        {
            field.row(i) = curMesh->optVars.W_opt.row(i).normalized();
        }

  //      viewer->data.add_edges(curMesh->centroids_F + field*avg / 2,
  //          curMesh->centroids_F, green);
    }
    else
    {
//        viewer->data.add_edges(curMesh->centroids_F + curMesh->optVars.W_opt*avg / 2,
//            curMesh->centroids_F, green);
    }

//    viewer->data.add_edges(curMesh->centroids_F  + curMesh->optVars.W_opt.normalized()*avg/2 * curMesh->vs.physicsEnergyArrowScale, curMesh->centroids_F, red);

 //   viewer->data.add_edges(centroids_F  + (Op_Grad)*avg/2*operator_scale, centroids_F, green);
//    viewer->data.add_edges(curMesh->centroids_F  + curMesh->v0*avg/2, curMesh->centroids_F, blue);

    extern Eigen::MatrixXd forceField;

//    viewer->data.add_edges(curMesh->V  + (forceField) * avg * curMesh->vs.physicsEnergyArrowScale, curMesh->V, red);
} 



