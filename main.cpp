#include <igl/viewer/Viewer.h>
#include <igl/avg_edge_length.h>

#include <math.h>
#include <iostream>
#include <fstream>

#include "DataLoad.h"
#include "Covariant.h"
#include "FaceBased.h"


Eigen::MatrixXd V;
Eigen::MatrixXi F, E, F_edges;

igl::viewer::Viewer *viewer;

double px = 0;
double py = 0;

void computeCovariantOperator(const Eigen::VectorXd &scalar_F,
    const Eigen::MatrixXi &F,
    const Eigen::MatrixXi &F_edges,
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &E,
    Eigen::MatrixXd &DM_local)
{
    Eigen::VectorXd scalar_E;
    computeEdgeWeights(scalar_F, V, E, scalar_E); 
    
    int nfaces = F.rows();

    DM_local.resize(nfaces, 2);

    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector3d u = V.row(F(i, 1)) - V.row(F(i, 0));
        Eigen::Vector3d v = V.row(F(i, 2)) - V.row(F(i, 0));

        double u_weight = 2 * (scalar_E(F_edges(i, 0)) - scalar_E(F_edges(i, 1)));
        double v_weight = 2 * (scalar_E(F_edges(i, 0)) - scalar_E(F_edges(i, 2)));

        DM_local(i, 0) = u_weight;
	DM_local(i, 1) = v_weight; 
    }
}

void evaluateOperator(const Eigen::MatrixXd &DM_local,
	const Eigen::MatrixXd &W_local, 
        Eigen::MatrixXd &del_W_F, int idx)
{
    Eigen::Vector3d component(0, 0, 0);
    component(idx) = 1;

// This is super slow because it does a bunch of irrelevant multiplications...    
//    Eigen::MatrixXd weights = DM_local * W_local.transpose();

    for (int i = 0; i < DM_local.rows(); i++)
    {
 
 	del_W_F.row(i) += component * (DM_local(i, 0) * W_local(i, 0) + DM_local(i,1) * W_local(i, 1));
 
    }
}

void computeCovariantOperatorNew(const Eigen::MatrixXi &F, 
	const Eigen::MatrixXd &V, 
	const Eigen::MatrixXi &E, 
	const Eigen::MatrixXi &F_edges, 
        const Eigen::MatrixXd &v, 
	const Eigen::MatrixXd &w,
        const std::vector<Eigen::SparseMatrix<double> > &Ms,
        Eigen::MatrixXd &result)
{
    result.resize(v.rows(), v.cols());
    int nfaces = F.rows();
    for (int i = 0; i < nfaces; i++)
    {
        result.row(i) = w.row(i) * (Ms[i] * v);
    }
}

Eigen::MatrixXd W; // This can be thought of as 3 ``independent'' scalar fields
Eigen::MatrixXd W_test; // The derivative is tested in this direction
std::vector<Eigen::SparseMatrix<double> > Ms; // the gradient operator; Ms[i] * F gives the gradient of F on triangle i

Eigen::MatrixXd colorField;
Eigen::MatrixXd centroids_F;

Eigen::MatrixXd del_W_F;

void computeOperatorGradient( const std::vector<Eigen::SparseMatrix<double> > &Ms,
	const Eigen::MatrixXd &del_W_F,
	const Eigen::MatrixXd &v,
	Eigen::MatrixXd &Op_Grad)
{
    int nfaces = v.rows();

    Op_Grad = Eigen::MatrixXd::Zero(nfaces, 3);

    for (int i = 0; i < nfaces; i++)
    {
        Op_Grad += Ms[i].transpose() * v.row(i).transpose() * del_W_F.transpose();
//        Op_Grad.row(i) += del_W_F.transpose() * v.row(i).transpose() * Ms[i].transpose();
    }
}



double energy_OP = 0.;
void updateView(const Eigen::MatrixXd del_W_F, int step)
{
    int nFaces = F.rows(); 
    char derFName[50];
    sprintf(derFName, "der_step_%d.txt", step);
    std::ofstream myfile (derFName);


    // Set mesh colors and log operator state
    Eigen::VectorXd Z(nFaces);
    energy_OP = 0.;
//    std::cout << del_W_F.rows() << " " << del_W_F.cols() << "\n";
    for (int i = 0; i < nFaces; i++)
    {
        Z(i) = log(del_W_F.row(i).norm() + .000005);

	if (myfile.is_open())
	{
	    myfile << del_W_F.row(i) << "\n";
            energy_OP += del_W_F.row(i).squaredNorm();       
	}
	else std::cout << "Unable to open file";
    }
    myfile.close();
    std::cout << energy_OP << "\n";

    // Average edge length for sizing
    const double avg = igl::avg_edge_length(V,F);
    colorField.resize(nFaces, 3);
    
    //  igl::jet(Z,true,colorField);
    igl::colormap(igl::COLOR_MAP_TYPE_MAGMA,Z, true, colorField);


    // Plot the mesh
    viewer->data.clear();
    viewer->data.set_mesh(V, F);
    viewer->data.set_face_based(true);

    viewer->data.set_colors(colorField);

    Eigen::MatrixXd eps = Eigen::MatrixXd::Constant(nFaces,3,.001);

    Eigen::MatrixXd field;
    field.resize(nFaces, 3);
    for (int i = 0; i < nFaces; i++)
    {
	field.row(i) = W.row(i).normalized();
    }

    const Eigen::RowVector3d green(0.1,0.9,0.2),blue(0.1,0.2,0.8);
    viewer->data.add_edges(centroids_F  + del_W_F*avg/2, centroids_F, blue);
    viewer->data.add_edges(centroids_F  + field*avg/2, centroids_F, green);
}

int descentStep = 0;
void takeGradientDescentStep()
{
    Eigen::MatrixXd W_local;
    computeLocalCoordinatesForDistanceField(W, F, V, W_local);

    del_W_F.resize(F.rows(), 3);
    del_W_F.setZero();
    Eigen::MatrixXd Op_Grad;

    // Not effecient, but will make it feel more correct to update, then show
    for (int i = 0; i < 1; i++) 
    {
     //   computeCovariantOperator(W.col(i), F, F_edges, V, E, DM_local);	
//	updateWInGradientDirection(DM_local, W, i);
  
        computeCovariantOperatorNew(F, V, E, F_edges, W, W_test, Ms, del_W_F);
        computeOperatorGradient(Ms, del_W_F, W, Op_Grad);

        W -= Op_Grad * .001;

    }

    computeCovariantOperatorNew(F, V, E, F_edges, W, W_test, Ms, del_W_F);
       
    descentStep++;
    updateView(del_W_F, descentStep);
}

void showVectorField()
{
    computeCentroids(F,V,centroids_F);

    Eigen::Vector3d p(px, py,0);
    computeDistanceField(p, centroids_F, W);
//    computeWhirlpool(p, centroids_F, W);

    computeDistanceField(p, centroids_F, W_test);
//    computeTestField(p, centroids_F, W_test);

    computeCovariantOperatorNew(F, V, E, F_edges, W, W_test, Ms, del_W_F);

    int nFaces = F.rows(); 

    //  Eigen::VectorXd Z = W.col(0); // - del_W_F;// - W_recovered.col(0);
    // Eigen::VectorXd Z = del_W_F.transpose() * del_W_F;// - W_recovered.col(0);

    descentStep = 0;
    updateView(del_W_F, descentStep);   

}


int main(int argc, char *argv[])
{  
  //   assignFaceVal(F,viz);;

  igl::readOBJ("../circ.obj", V, F);
  buildEdges(F, E);
  buildEdgesPerFace(F, E, F_edges);
  computeGradientMatrices(F, V, E, F_edges, Ms);

  // Plot the mesh  
  viewer = new igl::viewer::Viewer();
  viewer->data.set_mesh(V, F);
  viewer->data.set_face_based(true);
  viewer->callback_init = [&](igl::viewer::Viewer& viewer)
  {
      // Add new group
      viewer.ngui->addGroup("Vector Field Options");

      // Expose a variable
      viewer.ngui->addVariable("Center X",px);
      viewer.ngui->addVariable("Center Y",py);

      // Add a button
      viewer.ngui->addButton("Recompute Derivative", showVectorField);
      viewer.ngui->addButton("Grad Descent Step", takeGradientDescentStep);

      // call to generate menu
      viewer.screen->performLayout();
      return false;
  };

  viewer->launch();
}
