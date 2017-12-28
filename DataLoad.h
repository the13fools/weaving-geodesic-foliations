#ifndef DATALOAD_H
#define DATALOAD_H

#include <Eigen/Core>


Eigen::MatrixXd readMatrix(const char *filename);

void computeCentroids(const Eigen::MatrixXi &F,const Eigen::MatrixXd &V, Eigen::MatrixXd &centroids);
void computeDistanceField(const Eigen::Vector3d p, const Eigen::MatrixXd &centroids, Eigen::MatrixXd &W);
void computeWhirlpool(const Eigen::Vector3d p, const Eigen::MatrixXd &centroids, Eigen::MatrixXd &W);
void computeTestField(const Eigen::Vector3d p, const Eigen::MatrixXd &centroids, Eigen::MatrixXd &W);

void buildEdges(const Eigen::MatrixXi &F, Eigen::MatrixXi &E);


void buildEdgesPerFace(const Eigen::MatrixXi &F, const Eigen::MatrixXi &E, Eigen::MatrixXi &F_edge);

#endif
