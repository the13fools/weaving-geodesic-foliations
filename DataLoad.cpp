#include "DataLoad.h"
#include <Eigen/Geometry>

#include <igl/viewer/Viewer.h>
using namespace Eigen;
using namespace std;



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
