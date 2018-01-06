#ifndef FILEIO_H
#define FILEIO_H

#include <Eigen/Core>
#include "StateStructs.h"

Eigen::MatrixXd readMatrix(const char *filename);

void logRibbonsToFile(const VisualizationState &vs, std::string foldername, std::string filename);

void logToFile(const Eigen::MatrixXd W, std::string foldername, std::string filename);

#endif
