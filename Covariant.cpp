
#include <math.h>
#include <iostream>
#include <fstream>

#include "Covariant.h"


void computeEdgeWeights(const Eigen::VectorXd &scalar_F, 
                   const Eigen::MatrixXd &V,
                   const Eigen::MatrixXi &E, 
		   Eigen::VectorXd &scalar_E) 
{
    int nfaces = scalar_F.rows();
    int nedges = E.rows();
    
    scalar_E.resize(nedges);
    scalar_E = Eigen::VectorXd::Constant(nedges, 0);

    for (int i = 0; i < nedges; i++)
    {
	for (int j = 0; j < 2; j++)
	{
	    double factor = 0.5;
	    int faceId = E(i, j + 2);
            if (faceId > -1)
	    {
		if ( E(i, ( (j + 1) % 2) + 2) == -1) 
		{
		    factor = 1.;
		}
		scalar_E(i) += factor * scalar_F(faceId);
	    }	
	}
//	scalar_E(i) = .5 * (V(E(i,0),0) + V(E(i,1),0));
    }    
}

void computeLocalCoordinatesForDistanceField(const Eigen::MatrixXd &W, 
					     const Eigen::MatrixXi &F, 
					     const Eigen::MatrixXd &V, 
					     Eigen::MatrixXd &W_local)
{
    int nfaces = F.rows();
    
    W_local.resize(nfaces, 2);
    for (int i = 0; i < nfaces; i++) 
    {
    	Eigen::Vector3d u = V.row(F(i,1)) - V.row(F(i,0));
    	Eigen::Vector3d v = V.row(F(i,2)) - V.row(F(i,0));
	double uu = u.dot(u);
	double uv = u.dot(v);
	double vv = v.dot(v);
        
	Eigen::MatrixXd uvTrans(2,3);
	uvTrans << u.transpose(), v.transpose(); 
//	uvTrans = uvTrans.transpose();

        Eigen::MatrixXd uvInv(2,2);
        uvInv << vv, -uv, 
                -uv,  uu;
        uvInv *= 1. / (vv * uu - uv * uv);
        
//	std::cout << uvInv << "\n" << vv * uu <<  uv * uv << "\n" << uvTrans << "\n\n";
        Eigen::Vector2d alpha_beta = uvInv * uvTrans * W.row(i).transpose();
        W_local.row(i) = alpha_beta.transpose();
// 	std::cout << W_local.row(i) << "\n\n";
    }
        
}


void computeCovariantDerivative(const Eigen::MatrixXd &W_local,
    const Eigen::MatrixXi &F,
    const Eigen::MatrixXi &F_edges,
    const Eigen::MatrixXd &V,
    const Eigen::VectorXd &scalar_E,
    Eigen::MatrixXd &del_W_F, int idx)
{
    int nfaces = F.rows();

    std::cout << scalar_E.size() << "\n";

    Eigen::Vector3d component(0, 0, 0);
    component(idx) = 1;


    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector3d u = V.row(F(i, 1)) - V.row(F(i, 0));
        Eigen::Vector3d v = V.row(F(i, 2)) - V.row(F(i, 0));

        double u_weight = 2 * (scalar_E(F_edges(i, 0)) - scalar_E(F_edges(i, 1)));
        double v_weight = 2 * (scalar_E(F_edges(i, 0)) - scalar_E(F_edges(i, 2)));


        double comp_weight = u_weight * W_local(i, 0) + v_weight * W_local(i, 1);

        del_W_F.row(i) += component * comp_weight;
    }
}

void computeRecoveredDistanceField_test(const Eigen::MatrixXd &W_local, 
					const Eigen::MatrixXi &F, 
			                const Eigen::MatrixXd &V, 
				  	      Eigen::MatrixXd &W_recovered)
{
    int nfaces = F.rows();

    std::cout << W_local.size() << "\n";    
    
    W_recovered.resize(nfaces, 3);
    for (int i = 0; i < nfaces; i++) 
    {
    	Eigen::Vector3d u = V.row(F(i,1)) - V.row(F(i,0));
    	Eigen::Vector3d v = V.row(F(i,2)) - V.row(F(i,0));
//	std::cout << i << "\n";
//	std::cout << W_local.row(i) << "\n"; //  V.row(F(i,2)) - V.row(F(i,0))  << "\n";
        Eigen::Vector3d rec = u * W_local(i, 0) + v * W_local(i, 1);	    
        W_recovered.row(i) = rec;
    }
}

