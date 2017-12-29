#include <igl/viewer/Viewer.h>
#include <igl/avg_edge_length.h>

#include "RelaxViewer.h"

double energy_OP = 0.;
Eigen::MatrixXd colorField;

shading_enum shading_enum_state = INIT_DIRECTION; // Make this better...

void updateView(const MeshData &data, const OptState &state, igl::viewer::Viewer *viewer)
{
    int nFaces = data.F.rows(); 
    
    Eigen::MatrixXd W = state.W;
    Eigen::MatrixXd W_init = data.W_init;
 //   logToFile(W, folderName, std::to_string(step)); 


//    if (!is_fd_set) 
//    {
//	is_fd_set = true;
 //       computeOperatorGradient(Ms, del_W_F, W, Op_Grad_fd);
//        loadFiniteDifference();
//	std::cout << Op_Grad_fd.rows() << " " << Op_Grad_fd.cols() << "\n";
//    }

    Eigen::MatrixXd Op_Grad;
//    dvEnergy(data, W, W, Op_Grad);
//    logToFile(Op_Grad, folderName, "op_grad"); 

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
//		std::cout << Z(i) << "\n";
                if (Z(i) > max_error) 
		{
		    std::cout << i << " " << Z(i) << "\n";
		    max_error = Z(i);
		}
		break;
            case INIT_DIRECTION:
		Z(i) = (W-W_init).row(i).normalized().dot(testDir);
//		Z(i) = (Op_Grad).row(i).normalized().dot(testDir) + .000005;
//		Z(i) = (Op_Grad).row(i).norm();
		break;
            case INIT_MAGNITUDE:
//		Z(i) = log( (W-W_init).row(i).squaredNorm() );
//              Eigen::Vector3d test_vec(-Op_Grad(i,1), Op_Grad(i,0), 0);
//              Z(i) = (Op_Grad_fd).row(i).normalized()
//                         .dot(test_vec.normalized()) + .000005;
        //        std::cout << Z(i) << "\n";
//              Z(1) = 1;
                break;
        }

//        energy_OP += del_W_F.row(i).squaredNorm();
    }
//    std::cout << energy_OP << " Operator Energy\n";

    // Average edge length for sizing
    const double avg = igl::avg_edge_length(data.V,data.F);
    colorField.resize(nFaces, 3);

    //  igl::jet(Z,true,colorField);

    switch (shading_enum_state)
    {
        case OP_ERROR:
            igl::colormap(igl::COLOR_MAP_TYPE_INFERNO,Z, true, colorField);
            break;
        case INIT_DIRECTION:
            igl::colormap(igl::COLOR_MAP_TYPE_INFERNO,Z, true, colorField);
            break;
        case INIT_MAGNITUDE:
            igl::colormap(igl::COLOR_MAP_TYPE_INFERNO,Z, true, colorField);
            break; // MAGMA, JET
    }


    // Plot the mesh
    viewer->data.clear();
    viewer->data.set_mesh(data.V, data.F);
    viewer->data.set_face_based(true);

    viewer->data.set_colors(colorField);

    Eigen::MatrixXd eps = Eigen::MatrixXd::Constant(nFaces,3,.001);
/*
    Eigen::MatrixXd field;
    field.resize(nFaces, 3);
    for (int i = 0; i < nFaces; i++)
    {
        field.row(i) = W.row(i).normalized();
    }
*/
    const Eigen::RowVector3d red(0.9,.1,.1),green(0.1,0.9,0.2),blue(0.1,0.2,0.8);
  //  viewer->data.add_edges(centroids_f  + del_w_f*avg/2, centroids_f, blue);

//    viewer->data.add_edges(centroids_F  + (W - W_init)*avg/2*operator_scale, centroids_F, red);
//    viewer->data.add_edges(centroids_F  + (Op_Grad*operator_scale - Op_Grad_fd)*avg/2, centroids_F, red);
    viewer->data.add_edges(data.centroids_F  + W*avg/2, data.centroids_F, blue);
//    viewer->data.add_edges(data.centroids_F  + (Op_Grad)*avg/2*operator_scale, data.centroids_F, green);
}















