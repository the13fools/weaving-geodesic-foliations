#include "InitField.h"
#include <Eigen/Geometry>
#include <cassert>
#include <iostream>

#include <igl/viewer/Viewer.h>

#include "VectorUtils.h"
#include <queue>

using namespace Eigen;
using namespace std;


void propogateField(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const Eigen::MatrixXi &F_edges, Eigen::MatrixXd &field)
{
  std::queue<int> toProp;

  for (int i = 0; i < F.rows(); i++)
  {
      if (field.row(i).norm() > .01) 
          toProp.push(i);
  }


  while (!toProp.empty())
  {
      int i = toProp.front();
      toProp.pop();
      for (int e = 0; e < 3; e++)
      {
	  int edgeIdx  = F_edges(i, e);
	  int neighbor = E( edgeIdx, 2 );
	  if (neighbor == i) { neighbor = E( edgeIdx, 3 ); }
	  if (neighbor != -1 && field.row(neighbor).norm() < .01)
	  {
	      toProp.push(neighbor);
	      field.row(neighbor) = mapVectorToAdjacentFace(F, V, E, 
		                        edgeIdx, i, neighbor, field.row(i));
	  } 
      }
  }
}


void computeTestField(const Eigen::Vector3d p, const Eigen::MatrixXd &centroids, Eigen::MatrixXd &W)
{
    int nfaces = centroids.rows();
    
    W.resize(nfaces, 3);
    for (int i = 0; i < nfaces; i++) 
    { 
        Eigen::Vector3d blah = centroids.row(i);
        W(i, 0) = 1.;
        W(i, 1) = 0.;//blah(0);
        W(i, 2) = 0.;

    }
}

void computeWhirlpool(const Eigen::Vector3d p, const Eigen::MatrixXd &centroids, Eigen::MatrixXd &W)
{
    int nfaces = centroids.rows();
    
    W.resize(nfaces, 3);
    for (int i = 0; i < nfaces; i++) 
    { 
        Eigen::Vector3d blah = centroids.row(i);
        W(i, 0) = blah(0);
        W(i, 1) = 0.;//blah(0);
        W(i, 2) = 0.;

    }
}

void computeDistanceField(const Eigen::Vector3d p, const Eigen::MatrixXd &centroids, Eigen::MatrixXd &W)
{
    int nfaces = centroids.rows();
    
    W.resize(nfaces, 3);
    for (int i = 0; i < nfaces; i++) 
    { 
        Eigen::Vector3d vec = -centroids.row(i);
	vec += p;
        vec.normalize();
//	Eigen::Vector3d eps = Eigen::Vector3d::Random();
//	W.row(i) = vec + eps * .01;
	W.row(i) = vec;
    }
}
