#include "DataLoad.h"
#include <Eigen/Geometry>
#include <cassert>
#include <iostream>

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


void buildEdges(const Eigen::MatrixXi &F, Eigen::MatrixXi &E)
{
    map<pair<int, int>, Vector4i, std::less<pair<int, int> >,
        Eigen::aligned_allocator<std::pair<const int, Eigen::Vector4i> >>
          edgemap;

    int nfaces = F.rows();
    for (int i = 0; i < nfaces; i++)
    {
	for (int j = 0; j < 3; j++)
	{
	    int idx1 = F(i, j);
	    int idx2 = F(i, (j + 1) % 3);
	    int slot = 0;

	    if (idx1 > idx2)
	    {
		swap(idx1, idx2);
		slot = 1;
	    }
	    map<pair<int, int>, Vector4i, std::less<pair<int, int> >,
		Eigen::aligned_allocator<std::pair<const int, Eigen::Vector4i> >>::iterator it = edgemap.find(pair<int, int>(idx1, idx2));
	    if (it == edgemap.end())
	    {
		Vector4i newedge;
		newedge[0] = idx1;
		newedge[1] = idx2;
		newedge[2] = newedge[3] = -1;
		newedge[2 + slot] = i;
		edgemap[pair<int, int>(idx1, idx2)] = newedge;
	    }
	    else
	    {
		 edgemap[pair<int, int>(idx1, idx2)][2 + slot] = i;
	    }
	}
    }

    int nedges = edgemap.size();
    E.resize(nedges, 4);
    int idx = 0;
    for (map<pair<int, int>, Vector4i>::iterator it = edgemap.begin(); it != edgemap.end(); ++it)
    {
        E.row(idx) = it->second.transpose();
        idx++;
    }
}
