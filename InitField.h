#ifndef INITFIELD_H
#define INITFIELD_H

#include <Eigen/Core>
#include "StateStructs.h"


void computeDistanceField(const Eigen::Vector3d p, const Eigen::MatrixXd &centroids, Eigen::MatrixXd &W);
void computeWhirlpool(const Eigen::Vector3d p, const Eigen::MatrixXd &centroids, Eigen::MatrixXd &W);
void computeTestField(const Eigen::Vector3d p, const Eigen::MatrixXd &centroids, Eigen::MatrixXd &W);

#endif
