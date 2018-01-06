#ifndef COVARIANT_H
#define COVARIANT_H

#include <Eigen/Core>

void computeEdgeWeights(const Eigen::VectorXd &scalar_F, 
                   const Eigen::MatrixXd &V,
                   const Eigen::MatrixXi &E, 
		   Eigen::VectorXd &scalar_E); 

void computeLocalCoordinatesForDistanceField(const Eigen::MatrixXd &W, 
					     const Eigen::MatrixXi &F, 
					     const Eigen::MatrixXd &V, 
					     Eigen::MatrixXd &W_local);

void computeCovariantDerivative(const Eigen::MatrixXd &W_local,
    const Eigen::MatrixXi &F,
    const Eigen::MatrixXi &F_edges,
    const Eigen::MatrixXd &V,
    const Eigen::VectorXd &scalar_E,
    Eigen::MatrixXd &del_W_F, int idx);

void computeRecoveredDistanceField_test(const Eigen::MatrixXd &W_local, 
					const Eigen::MatrixXi &F, 
			                const Eigen::MatrixXd &V, 
				  	      Eigen::MatrixXd &W_recovered);

#endif 
