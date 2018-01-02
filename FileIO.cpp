#include <Eigen/Core>
#include <fstream>

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
