#include <igl/avg_edge_length.h>
#include <igl/viewer/Viewer.h>

#include "RelaxViewer.h"
#include "VectorUtils.h"

double energy_OP = 0.;
//void updateView(const MeshData &curMesh, igl::viewer::Viewer &viewer){}


Eigen::MatrixXd colorField;

shading_enum shading_enum_state = INIT_MAGNITUDE;

/*
// This: Recieves a point on an edge of a triangle
//       Returns the vertex id opposing that edge in the triangle
int getOpVId(const MeshData &md, const Eigen::Vector3d prev_point, int faceId)
{
    for (int j = 0; j < 3; j++)
    {
	Eigen::VectorXi e = md.E.row(md.F_edges(faceId, j));
	Eigen::Vector3d e_test =  md.V.row( e(0) ) - md.V.row( e(1) );
	e_test.normalize();
	if ( e_test.dot( (prev_point - e_test ).normalized() ) > .9 )
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

}
 

void traceCurve(const MeshData &md, const Eigen::Vector3d dir, int faceId, VisualizationState &vs)
{
    vs.curve.clear();
    vs.curve.push_back(.5 * md.V.row(md.F(faceId, 0)) + .5 * md.V.row(md.F(faceId, 1)) ); 
    // check that vector is pointing in.  
    // TODO
    // find next edge
    int curr_face_id = faceId;
    Eigen::Vector3d cur_dir = dir;
    for (int i = 0; i < 200; i++)
    {
        Eigen::Vector3d prev_point = vs.curve.back();
	int op_v_id = getOpVId(md, prev_point, curr_face_id);
            
        Eigen::Vector3d split = prev_point - md.V.row(op_v_id);
        split.normalize();
	Eigen::Vector3d n = faceNormal(md.F, md.V, curr_face_id);
	Eigen::Vector3d perp = split.cross(n);

	int op_edge_id = 0;
        for (int j = 0; j < 3; j++)
	{
	    Eigen::VectorXi e = md.E.row(md.F_edges(faceId, j));
	    if ( e(0) == op_v_id || e(1) == op_v_id )
            {
                Eigen::VectorXd e_test = md.V.row(e(0)) - md.V.row(e(1));
                if ( e_test.dot(perp) * cur_dir.dot(perp) > 0.) 
                {
                    op_edge_id = j;
                    break;
                }               
            }
        }
        
        // Find intersection point.
        

     }
}
*/


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
                Z(i) = log( (curMesh->optVars.W_opt-curMesh->v0).row(i).squaredNorm() );
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

    igl::ColorMapType viz_color = igl::COLOR_MAP_TYPE_INFERNO;

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
    if (curMesh->vs.normFaceVectors)
    {
	Eigen::MatrixXd field;
	field.resize(nFaces, 3);
	for (int i = 0; i < nFaces; i++)
	{
	    field.row(i) = curMesh->optVars.W_opt.row(i).normalized();
	}

	viewer->data.add_edges(curMesh->centroids_F + field*avg/2, 
		    curMesh->centroids_F, blue);
    }
    else
    {
        viewer->data.add_edges(curMesh->centroids_F + curMesh->optVars.W_opt*avg/2, 
		    curMesh->centroids_F, blue);
    }

    viewer->data.add_edges(curMesh->centroids_F  + curMesh->optVars.W_opt.normalized()*avg/2, curMesh->centroids_F, red);

 //   viewer->data.add_edges(centroids_F  + (Op_Grad)*avg/2*operator_scale, centroids_F, green);
    viewer->data.add_edges(curMesh->centroids_F  + curMesh->v0*avg/2, curMesh->centroids_F, green);
} 



