#include <igl/viewer/Viewer.h>
#include <igl/avg_edge_length.h>

#include "DataLoad.h"

void assignFaceVal(const Eigen::MatrixXi &F, Eigen::MatrixXd val) 
{
    int nfaces = F.rows();
 
    for (int i = 0; i < nfaces; i++)
    {
	for (int j = 0; j < 3; j++)
	{

	}
    }
}

void buildEdgesPerFace(const Eigen::MatrixXi &F, const Eigen::MatrixXi &E, Eigen::MatrixXi &F_edge)
{
    int nfaces = F.rows();
//    F_edge.resize(nfaces, 3);
    
       
}

/*

void computeCentroids(const Eigen::MatrixXi &F,const Eigen::MatrixXd &V, Eigen::MatrixXd &centroids)
{
 //   Eigen::MatrixXd interp; 
    int nfaces = F.rows();
    int nverts = V.rows();
    
    centroids.resize(nfaces, 3);
    for (int i = 0; i < nfaces; i++) 
    { 
	Eigen::Vector3d pos(0,0,0);
	for (int j = 0; j < 3; j++) 
	{
//	    std::cout << V.row(F(i,j)) << "\n" << F(i,j) << "\n\n";
	    pos += V.row(F(i,j));
	}
	centroids.row(i) = pos/3;
    }
}

void computeDistanceField(const Eigen::Vector3d p, const Eigen::MatrixXd &centroids, Eigen::MatrixXd &W)
{
    int nfaces = centroids.rows();
    
    W.resize(nfaces, 3);
    for (int i = 0; i < nfaces; i++) 
    { 
        Eigen::Vector3d blah = -centroids.row(i);
	blah += p;
	W.row(i) = blah.normalized();
    }
}

*/

int main(int argc, char *argv[])
{

  Eigen::MatrixXd V;
  Eigen::MatrixXi F, E;

  Eigen::MatrixXd viz, colorField;
  Eigen::MatrixXd centroids_F, W;
  //   assignFaceVal(F,viz);;

  igl::readOBJ("../circ.obj", V, F);
  buildEdges(F, E);
  computeCentroids(F,V,centroids_F);

  Eigen::Vector3d p(0,1,0);
  computeDistanceField(p, centroids_F, W);

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



  // Plot the mesh
  igl::viewer::Viewer viewer;
  viewer.data.set_mesh(V, F);
  viewer.data.set_face_based(true);


  const Eigen::RowVector3d red(0.8,0.2,0.2),blue(0.2,0.2,0.8);
  viewer.data.add_edges(centroids_F + W*avg/2, centroids_F, blue);
//  viewer.data.add_edges(centroids_F, centroids_F *2 * avg, blue);

  viewer.launch();
}
