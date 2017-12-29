#include "DataLoad.h"
#include <Eigen/Geometry>
#include <cassert>
#include <iostream>

#include <igl/viewer/Viewer.h>
using namespace Eigen;
using namespace std;


#define MAXBUFSIZE  ((int) 1e6)
Eigen::MatrixXd readMatrix(const char *filename)
{
    int cols = 0, rows = 0;
    double buff[MAXBUFSIZE];

   //  Read numbers from file into buffer.
    ifstream infile;
    infile.open(filename);
    while (! infile.eof())
    {
        string line;
        getline(infile, line);

        int temp_cols = 0;
        stringstream stream(line);
        while(! stream.eof())
            stream >> buff[cols*rows+temp_cols++];

        if (temp_cols == 0)
            continue;

        if (cols == 0)
            cols = temp_cols;

        rows++;
    }

    infile.close();

    rows--;

    // Populate matrix with numbers.
    Eigen::MatrixXd result(rows,cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result(i,j) = buff[ cols*i+j ];

    return result;
}


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
