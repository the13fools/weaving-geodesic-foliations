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


void logRibbonsToFile(const VisualizationState &vs, std::string foldername, std::string filename)
{
#ifndef WIN32
    char folderpath[50];
    sprintf(folderpath, "log/%s", foldername.c_str());
    mkdir("log", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(folderpath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    char logpath[50];
    sprintf(logpath, "%s/%s.txt", folderpath, filename.c_str());   
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

    std::vector<Eigen::Vector3d> cnew;
    cnew.push_back(vs.c0[0]);
    int seg_counter = 0;
    Eigen::Vector3d prev_point = cnew.back();
    for (int i = 1; i < vs.c0.size(); i++)
    {
        double seg_length = ( prev_point - vs.c0[i] ).norm(); 	
	if (seg_length > max_length) 
	{
            cnew.push_back(vs.c0[i]);
	    prev_point = cnew.back();
	}
    }

    std::vector<Eigen::Vector3d> nnew (cnew.size() - 1, Eigen::Vector3d::Zero());
    for (int i = 0; i < cnew.size() - 2; i++) 
    {
        Eigen::Vector3d v1 = cnew[i+1] - cnew[i];
        Eigen::Vector3d v2 = cnew[i+2] - cnew[i+1];
	Eigen::Vector3d perp = v1.cross(v2);
	nnew[i] = v1.cross(perp);
	nnew[i+1] = perp.cross(v2);
    }


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
	    myfile << " .1 .1 .1";
    } 

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
