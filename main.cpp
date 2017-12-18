#include <igl/viewer/Viewer.h>
#include <igl/avg_edge_length.h>

#include <math.h>

#include "DataLoad.h"

void computeEdgeWeights(const Eigen::VectorXd &scalar_F, 
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
			              Eigen::MatrixXd &del_W_F)
{
    int nfaces = F.rows();

    std::cout << scalar_E.size() << "\n";    
    
    del_W_F.resize(nfaces, 3);
    for (int i = 0; i < nfaces; i++) 
    {
    	Eigen::Vector3d u = V.row(F(i,1)) - V.row(F(i,0));
    	Eigen::Vector3d v = V.row(F(i,2)) - V.row(F(i,0));
        
//	std::cout << F_edges(i, 0) << " " << F_edges(i, 1) << " " << F_edges(i, 2) << "\n"; 
//	std::cout << scalar_E(F_edges(i, 0)) << " " << scalar_E(F_edges(i, 1)) << " " << scalar_E(F_edges(i, 2)) << "\n"; 

	double u_weight = 2 * ( scalar_E(F_edges(i, 0)) - scalar_E(F_edges(i, 1)) );
	double v_weight = 2 * ( scalar_E(F_edges(i, 0)) - scalar_E(F_edges(i, 2)) );
	
	Eigen::Vector3d rec = u * u_weight * W_local(i, 0) + v * v_weight * W_local(i, 1);	    
        del_W_F.row(i) = rec;
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

int main(int argc, char *argv[])
{

  Eigen::MatrixXd V;
  Eigen::MatrixXi F, E, F_edges;

  Eigen::MatrixXd colorField;
  //   assignFaceVal(F,viz);;

  igl::readOBJ("../circ.obj", V, F);
  buildEdges(F, E);
  buildEdgesPerFace(F, E, F_edges);

  Eigen::MatrixXd centroids_F;
  computeCentroids(F,V,centroids_F);
  
  Eigen::Vector3d p(.1,.8,0);
  Eigen::MatrixXd W;
  computeDistanceField(p, centroids_F, W);

  Eigen::VectorXd scalar_E;
  computeEdgeWeights(W.col(0), E, scalar_E); 

  Eigen::MatrixXd W_local;
  computeLocalCoordinatesForDistanceField(W, F, V, W_local);

  Eigen::MatrixXd del_W_F;
  computeCovariantDerivative(W_local, F, F_edges, V, scalar_E, del_W_F);
  
//  std::cout << del_W_F;
//  Eigen::MatrixXd W_recovered;
//  computeRecoveredDistanceField_test(W_local, F, V, W_recovered);

  // Average edge length for sizing
  const double avg = igl::avg_edge_length(V,F);
  
  int nFaces = F.rows(); 
  colorField.resize(nFaces, 3);

  Eigen::VectorXd Z(nFaces);
  for (int i = 0; i < nFaces; i++)
  {
      Z(i) = log(del_W_F.row(i).norm());
  }
//  Eigen::VectorXd Z = W.col(0); // - del_W_F;// - W_recovered.col(0);
 // Eigen::VectorXd Z = del_W_F.transpose() * del_W_F;// - W_recovered.col(0);
  
//  igl::jet(Z,true,colorField);
  igl::colormap(igl::COLOR_MAP_TYPE_PARULA,Z, true, colorField);


  // Plot the mesh
  igl::viewer::Viewer viewer;
  viewer.data.set_mesh(V, F);
  viewer.data.set_face_based(true);
 
  viewer.data.set_colors(colorField);

  Eigen::MatrixXd eps = Eigen::MatrixXd::Constant(nFaces,3,.001);

  const Eigen::RowVector3d red(0.8,0.2,0.2),blue(0.2,0.2,0.8);
  viewer.data.add_edges(centroids_F  - del_W_F* avg/2, centroids_F, blue);
//  viewer.data.add_edges(centroids_F + W*avg/2 - del_W_F, centroids_F + eps, blue);
//  viewer.data.add_edges(centroids_F, centroids_F *2 * avg, blue);

  viewer.launch();
}
