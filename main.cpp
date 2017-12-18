#include <igl/viewer/Viewer.h>
#include <igl/avg_edge_length.h>

#include "DataLoad.h"

void computeEdgeWeights(const Eigen::VectorXd &scalar_F, 
		   const Eigen::MatrixXi &E, 
		   Eigen::VectorXd scalar_E) 
{
    int nfaces = scalar_F.rows();
    int nedges = E.rows();
    
//    scalar_E.resize(nedges);
//    scalar_E = Eigen::VectorXd::Constant(nedges, 0);

    return ;
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


int main(int argc, char *argv[])
{

  Eigen::MatrixXd V;
  Eigen::MatrixXi F, E, F_edges;

  Eigen::MatrixXd viz, colorField;
  Eigen::MatrixXd centroids_F, W;
  //   assignFaceVal(F,viz);;

  igl::readOBJ("../circ.obj", V, F);
  buildEdges(F, E);
  buildEdgesPerFace(F, E, F_edges);


  computeCentroids(F,V,centroids_F);
  
  Eigen::Vector3d p(.1,.8,0);
  computeDistanceField(p, centroids_F, W);

  Eigen::VectorXd scalar_E;
  computeEdgeWeights(W.col(0), E, scalar_E); 

//  std::cout << centroids_F;

  // Average edge length for sizing
  const double avg = igl::avg_edge_length(V,F);

  
  int nVerts = F.rows();
  viz.resize(nVerts, 3);
  
  for (int i = 0; i < nVerts; i++)
  {
//     Eigen::Vector3d v(0,1,0);
      for (int j = 0; j < 3; j++)
      { 
          viz(i,j) = j; 
      }
  }

  int nFaces = F.rows();
  colorField.resize(nFaces, 3);
  Eigen::VectorXd Z = W.col(1);
  
  igl::jet(Z,true,colorField);


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
