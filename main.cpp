#include <igl/viewer/Viewer.h>
#include <igl/avg_edge_length.h>

#include "DataLoad.h"

void computeEdgeWeights(const Eigen::VectorXd &scalar_F, 
		   const Eigen::MatrixXi &E, 
		   Eigen::VectorXd scalar_E) 
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

void computeLocalCoordinatesForDistanceField(const Eigen::MatrixXd W, 
					     const Eigen::MatrixXi &F, 
					     const Eigen::MatrixXd &V, 
					     Eigen::MatrixXd W_local)
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
        W_local.row(i) = alpha_beta;
	std::cout << alpha_beta << "\n\n";
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
//  std::cout << centroids_F;

  // Average edge length for sizing
  const double avg = igl::avg_edge_length(V,F);

  
  int nVerts = F.rows();

  int nFaces = F.rows();
  colorField.resize(nFaces, 3);
  Eigen::VectorXd Z = W.col(1) + W.col(0);
  
//  igl::jet(Z,true,colorField);
  igl::colormap(igl::COLOR_MAP_TYPE_PARULA,Z, true, colorField);


  // Plot the mesh
  igl::viewer::Viewer viewer;
  viewer.data.set_mesh(V, F);
  viewer.data.set_face_based(true);
 
  viewer.data.set_colors(colorField);

  Eigen::MatrixXd eps = Eigen::MatrixXd::Constant(nFaces,3,.001);

  const Eigen::RowVector3d red(0.8,0.2,0.2),blue(0.2,0.2,0.8);
  viewer.data.add_edges(centroids_F + W*avg/2, centroids_F, blue);
  viewer.data.add_edges(centroids_F + W*avg/2, centroids_F + eps, blue);
//  viewer.data.add_edges(centroids_F, centroids_F *2 * avg, blue);

  viewer.launch();
}
