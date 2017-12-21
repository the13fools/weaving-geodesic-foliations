#include <igl/viewer/Viewer.h>
#include <igl/avg_edge_length.h>

#include <math.h>

#include "DataLoad.h"
#include "FaceBased.h"

void computeEdgeWeights_noop(const Eigen::VectorXd &scalar_F, 
                   const Eigen::MatrixXi &E, 
		   Eigen::VectorXd &scalar_E) 
{
    int nfaces = scalar_F.rows();
    int nedges = E.rows();
    
    scalar_E.resize(nedges);
    scalar_E = Eigen::VectorXd::Constant(nedges, 0);

    for (int i = 0; i < nedges; i++)
    {
	scalar_E(i) = 0.;
    }    
}


void computeEdgeWeights(const Eigen::VectorXd &scalar_F, 
                   const Eigen::MatrixXd &V,
                   const Eigen::MatrixXi &E, 
		   Eigen::VectorXd &scalar_E) 
{
    int nfaces = scalar_F.rows();
    int nedges = E.rows();
    
    scalar_E.resize(nedges);
    scalar_E = Eigen::VectorXd::Constant(nedges, 0);

    for (int i = 0; i < nedges; i++)
    {
	for (int j = 0; j < 2; j++)
	{
	    double factor = 0.5;
	    int faceId = E(i, j + 2);
            if (faceId > -1)
	    {
		if ( E(i, ( (j + 1) % 2) + 2) == -1) 
		{
		    factor = 1.;
		}
		scalar_E(i) += factor * scalar_F(faceId);
	    }	
	}
//	scalar_E(i) = .5 * (V(E(i,0),0) + V(E(i,1),0));
    }    
}

void computeLocalCoordinatesForDistanceField(const Eigen::MatrixXd &W, 
					     const Eigen::MatrixXi &F, 
					     const Eigen::MatrixXd &V, 
					     Eigen::MatrixXd &W_local)
{
    int nfaces = F.rows();
    
    W_local.resize(nfaces, 2);
    for (int i = 0; i < nfaces; i++) 
    {
    	Eigen::Vector3d u = V.row(F(i,1)) - V.row(F(i,0));
    	Eigen::Vector3d v = V.row(F(i,2)) - V.row(F(i,0));
	double uu = u.dot(u);
	double uv = u.dot(v);
	double vv = v.dot(v);
        
	Eigen::MatrixXd uvTrans(2,3);
	uvTrans << u.transpose(), v.transpose(); 
//	uvTrans = uvTrans.transpose();

        Eigen::MatrixXd uvInv(2,2);
        uvInv << vv, -uv, 
                -uv,  uu;
        uvInv *= 1. / (vv * uu - uv * uv);
        
//	std::cout << uvInv << "\n" << vv * uu <<  uv * uv << "\n" << uvTrans << "\n\n";
        Eigen::Vector2d alpha_beta = uvInv * uvTrans * W.row(i).transpose();
        W_local.row(i) = alpha_beta.transpose();
// 	std::cout << W_local.row(i) << "\n\n";
    }
        
}

void computeCovariantDerivative(const Eigen::MatrixXd &W_local,
    const Eigen::MatrixXi &F,
    const Eigen::MatrixXi &F_edges,
    const Eigen::MatrixXd &V,
    const Eigen::VectorXd &scalar_E,
    Eigen::MatrixXd &del_W_F, int idx)
{
    int nfaces = F.rows();

    std::cout << scalar_E.size() << "\n";

    Eigen::Vector3d component(0, 0, 0);
    component(idx) = 1;


    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector3d u = V.row(F(i, 1)) - V.row(F(i, 0));
        Eigen::Vector3d v = V.row(F(i, 2)) - V.row(F(i, 0));

        double u_weight = 2 * (scalar_E(F_edges(i, 0)) - scalar_E(F_edges(i, 1)));
        double v_weight = 2 * (scalar_E(F_edges(i, 0)) - scalar_E(F_edges(i, 2)));


        double comp_weight = u_weight * W_local(i, 0) + v_weight * W_local(i, 1);

        del_W_F.row(i) += component * comp_weight;
    }
}

void computeCovariantDerivativeNew(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const Eigen::MatrixXi &F_edges, const Eigen::MatrixXd &v, const Eigen::MatrixXd &w,
    Eigen::MatrixXd &result)
{
    std::vector<Eigen::SparseMatrix<double> > Ms;
    computeGradientMatrices(F, V, E, F_edges, Ms);
    result.resize(v.rows(), v.cols());
    int nfaces = F.rows();
    for (int i = 0; i < nfaces; i++)
    {
        result.row(i) = w.row(i) * Ms[i] * v;
    }
}

void computeRecoveredDistanceField_test(const Eigen::MatrixXd &W_local, 
					const Eigen::MatrixXi &F, 
			                const Eigen::MatrixXd &V, 
				  	      Eigen::MatrixXd &W_recovered)
{
    int nfaces = F.rows();

    std::cout << W_local.size() << "\n";    
    
    W_recovered.resize(nfaces, 3);
    for (int i = 0; i < nfaces; i++) 
    {
    	Eigen::Vector3d u = V.row(F(i,1)) - V.row(F(i,0));
    	Eigen::Vector3d v = V.row(F(i,2)) - V.row(F(i,0));
//	std::cout << i << "\n";
//	std::cout << W_local.row(i) << "\n"; //  V.row(F(i,2)) - V.row(F(i,0))  << "\n";
        Eigen::Vector3d rec = u * W_local(i, 0) + v * W_local(i, 1);	    
        W_recovered.row(i) = rec;
    }
}

Eigen::MatrixXd V;
Eigen::MatrixXi F, E, F_edges;

igl::viewer::Viewer *viewer;

double px = 0;
double py = 0;

void showVectorField()
{
    Eigen::MatrixXd colorField;
    Eigen::MatrixXd centroids_F;
    computeCentroids(F,V,centroids_F);

    Eigen::Vector3d p(px, py,0);
    Eigen::MatrixXd W;
    computeDistanceField(p, centroids_F, W);
    //computeWhirlpool(p, centroids_F, W);


    Eigen::MatrixXd W_test;
    computeDistanceField(p, centroids_F, W_test);
    //computeTestField(p, centroids_F, W_test);

    Eigen::MatrixXd W_local;
    computeLocalCoordinatesForDistanceField(W_test, F, V, W_local);

    Eigen::MatrixXd del_W_F;
    del_W_F.resize(F.rows(), 3);
    del_W_F.setZero();
    Eigen::VectorXd scalar_E;
    for (int i = 0; i < 3; i++) 
    {
        computeEdgeWeights(W.col(i), V, E, scalar_E); 
        computeCovariantDerivative(W_local, F, F_edges, V, scalar_E, del_W_F, i);
    }

    Eigen::MatrixXd covresult;
    computeCovariantDerivativeNew(F, V, E, F_edges, W, W_test, covresult);
    int nfaces = F.rows();
    
    //  std::cout << del_W_F;
    //  Eigen::MatrixXd W_recovered;
    //  computeRecoveredDistanceField_test(W_local, F, V, W_recovered);

    // Average edge length for sizing
    const double avg = igl::avg_edge_length(V,F);

    int nFaces = F.rows(); 
    colorField.resize(nFaces, 3);

    Eigen::VectorXd Z(nFaces);
    double maxerror = 0;
    for (int i = 0; i < nFaces; i++)
    {
        Z(i) = log(del_W_F.row(i).norm());
        if (maxerror < Z(i))
        {
            maxerror = Z(i);
            std::cout << del_W_F.row(i) << "\n";
        }
    }
    //  Eigen::VectorXd Z = W.col(0); // - del_W_F;// - W_recovered.col(0);
    // Eigen::VectorXd Z = del_W_F.transpose() * del_W_F;// - W_recovered.col(0);

    //  igl::jet(Z,true,colorField);
    igl::colormap(igl::COLOR_MAP_TYPE_MAGMA,Z, true, colorField);


    // Plot the mesh
    viewer->data.clear();
    viewer->data.set_mesh(V, F);
    viewer->data.set_face_based(true);

    viewer->data.set_colors(colorField);

    Eigen::MatrixXd eps = Eigen::MatrixXd::Constant(nFaces,3,.001);

    const Eigen::RowVector3d red(0.8,0.2,0.2),blue(0.2,0.2,0.8);
    viewer->data.add_edges(centroids_F  + covresult*avg/2, centroids_F, blue);
}

int main(int argc, char *argv[])
{  
  //   assignFaceVal(F,viz);;

  igl::readOBJ("../circ.obj", V, F);
  buildEdges(F, E);
  buildEdgesPerFace(F, E, F_edges);

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

      // call to generate menu
      viewer.screen->performLayout();
      return false;
  };

  viewer->launch();
}
