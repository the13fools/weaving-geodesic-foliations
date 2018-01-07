#include <Eigen/Core>
#include <fstream>
#include <iostream>

#include "StateStructs.h"

// For making directories
#include <sys/stat.h>
// #include <direct.h>

#include <Eigen/Geometry>

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


Eigen::Vector3d parallelTransport(const Eigen::Vector3d &v, const Eigen::Vector3d &e1, const Eigen::Vector3d &e2)
{
    Eigen::Vector3d t1 = e1 / e1.norm();
    Eigen::Vector3d t2 = e2 / e2.norm();
    Eigen::Vector3d n = t1.cross(t2);
    if (n.norm() < 1e-8)
        return v;
    n /= n.norm();
    Eigen::Vector3d p1 = n.cross(t1);
    Eigen::Vector3d p2 = n.cross(t2);
    return v.dot(n)*n + v.dot(t1)*t2 + v.dot(p1)*p2;
}

void logRibbonsToFile(const VisualizationState &vs, std::string foldername, std::string filename)
{
#ifndef WIN32
    char folderpath[50];
    sprintf(folderpath, "log/%s", foldername.c_str());
    mkdir("log", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(folderpath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    char logpath[50];
    sprintf(logpath, "%s/%s.rod", folderpath, filename.c_str());   
    std::ofstream myfile (logpath);
    
    if (!myfile.is_open())
    {
	std::cout << "Unable to open file";
	return;
    }


    double max_length = 0.;
    for (int i = 0; i < vs.c0.size() - 1; i++)
    {
        double seg_length = ( vs.c0[i] - vs.c0[i+1] ).norm(); 
	if (seg_length > max_length) 
	{
            max_length = seg_length;
	}
    }

    // Decimate 
    std::vector<Eigen::Vector3d> cnew;
    cnew.push_back(vs.c0[0]);
    int seg_counter = 0;
    Eigen::VectorXd desc_mapping = Eigen::VectorXd::Zero(vs.c0.size());
    Eigen::Vector3d prev_point = cnew.back();
    for (int i = 0; i < vs.c0.size(); i++)
    {
        double seg_length = ( prev_point - vs.c0[i] ).norm(); 	
	if (seg_length > max_length) 
	{
            cnew.push_back(vs.c0[i]);
	    prev_point = cnew.back();
	    seg_counter++;
	}
	desc_mapping(i) = seg_counter;
    }
    cnew.pop_back();
    cnew.push_back(vs.c0[ vs.c0.size() - 1 ]);

    // Find new normal vector - could be improved by using old normals.
    std::vector<Eigen::Vector3d> nnew (cnew.size() - 1, Eigen::Vector3d::Zero());
    Eigen::Vector3d v1 = cnew[1] - cnew[0];
    Eigen::Vector3d v2 = cnew[2] - cnew[1];
    Eigen::Vector3d perp = v1.cross(v2);

    nnew[0] = v1.cross(perp).normalized();
    Eigen::Vector3d prev_norm = nnew[0];
    for (int i = 0; i < cnew.size() - 2; i++) 
    {
        Eigen::Vector3d v1 = cnew[i+1] - cnew[i];
        Eigen::Vector3d v2 = cnew[i+2] - cnew[i+1];
	Eigen::Vector3d new_normal = parallelTransport(prev_norm, v1, v2);
   
/*	Eigen::Vector3d v1 = cnew[i + 1] - cnew[i];
	Eigen::Vector3d v2 = cnew[i + 2] - cnew[i + 1];
	Eigen::Vector3d perp = v1.cross(v2);

	nnew[i + 1] = v2.cross(perp).normalized();
        if (nnew[i+1].dot(nnew[i]) < 0) { nnew[i+1] *= -1.; }
*/	nnew[i+1] = new_normal;
	prev_norm = new_normal;
//	std::cout << nnew[i].dot(cnew[i] - cnew[i+1]) << " ";
//	std::cout << nnew[i+1].dot(cnew[i+2] - cnew[i+1]) << "\n";

    }
//    nnew[0] = vs.n0[0];

    myfile << "1\n";
    myfile << "0\n";
    myfile << "0.0001\n";
    myfile << "1e+08\n"; 
    myfile << "1\n\n\n";
    
    myfile << cnew.size() << "\n";
    myfile << "0\n";
    for(int i = 0; i < cnew.size(); i++)
    {
	    myfile << cnew[i](0) << " " << cnew[i](1) << " " << cnew[i](2) << " ";
    }
    myfile << "\n";
     
    for(int i = 0; i < nnew.size(); i++)
    {
	    myfile << nnew[i](0) << " " << nnew[i](1) << " " << nnew[i](2) << " ";
    } 
    myfile << "\n";

    for(int i = 0; i < nnew.size(); i++)
    {
	double c = .2;
	    myfile << " .1";
    //        myfile << (double) i * .1  + c<< " " << (double) i * .1 + c << " " << (double) i * .1 +c << " ";
    } 
    myfile << "\n";

    myfile.close();
#endif
}


void logToFile(const Eigen::MatrixXd W, std::string foldername, std::string filename)
{
#ifndef WIN32
    char folderpath[50];
    sprintf(folderpath, "log/%s", foldername.c_str());
    mkdir("log", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(folderpath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    char logpath[50];
    sprintf(logpath, "%s/%s.txt", folderpath, filename.c_str());   
    std::ofstream myfile (logpath);
    for(int i = 0; i < W.rows(); i++)
    {
        if (myfile.is_open())
        {
	    myfile << W.row(i) << "\n";
	}

	else
	{
	    std::cout << "Unable to open file";
	    break;
	}
    } 
    myfile.close();
#endif
}
