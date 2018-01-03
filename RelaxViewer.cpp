#include <igl/avg_edge_length.h>
#include <igl/viewer/Viewer.h>

#include "RelaxViewer.h"


double energy_OP = 0.;
//void updateView(const MeshData &curMesh, igl::viewer::Viewer &viewer){}


Eigen::MatrixXd colorField;

shading_enum shading_enum_state = INIT_MAGNITUDE;

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

        viewer->data.add_edges(curMesh->centroids_F + field*avg / 2,
            curMesh->centroids_F, blue);
    }
    else
    {
        viewer->data.add_edges(curMesh->centroids_F + curMesh->optVars.W_opt*avg / 2,
            curMesh->centroids_F, blue);
    }

    viewer->data.add_edges(curMesh->centroids_F  + curMesh->optVars.W_opt.normalized()*avg/2, curMesh->centroids_F, red);

 //   viewer->data.add_edges(centroids_F  + (Op_Grad)*avg/2*operator_scale, centroids_F, green);
    viewer->data.add_edges(curMesh->centroids_F  + curMesh->v0*avg/2, curMesh->centroids_F, green);

    extern Eigen::MatrixXd forceField;

    viewer->data.add_edges(curMesh->V  + forceField, curMesh->V, red);
} 



