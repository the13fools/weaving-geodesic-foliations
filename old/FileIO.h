#ifndef FILEIO_H
#define FILEIO_H

#include <Eigen/Core>


Eigen::MatrixXd readMatrix(const char *filename);

void logToFile(const Eigen::MatrixXd W, std::string foldername, std::string filename);

#endif
