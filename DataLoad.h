#ifndef DATALOAD_H
#define DATALOAD_H

#include <Eigen/Core>



void computeCentroids(const Eigen::MatrixXi &F,const Eigen::MatrixXd &V, Eigen::MatrixXd &centroids);
void computeDistanceField(const Eigen::Vector3d p, const Eigen::MatrixXd &centroids, Eigen::MatrixXd &W);

void buildEdges(const Eigen::MatrixXi &F, Eigen::MatrixXi &E);



#endif
