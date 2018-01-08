#ifndef VECTORUTILS_H
#define VECTORUTILS_H

#include <Eigen/Core>
#include <Eigen/Sparse>

Eigen::Vector3d faceNormal(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, int faceidx);


Eigen::Vector3d mapVectorToAdjacentFace(const Eigen::MatrixXi &F, 
                                        const Eigen::MatrixXd &V,
                                        const Eigen::MatrixXi &E, 
                                        int edgeIdx, int init_fidx, int neighbor_fidx,
                                        const Eigen::Vector3d &inp);

#endif
