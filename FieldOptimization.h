#ifndef FIELDOPTIMIZATION_H
#define FIELDOPTIMIZATION_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

// stuff that is precomputed about the mesh and doesn't change during optimization
struct MeshData
{
    MeshData(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
    Eigen::MatrixXd V; // vertex positions
    Eigen::MatrixXi F; // face indices
    Eigen::MatrixXi E; // edge indices
    Eigen::MatrixXi F_edges; // face edge indices
    std::vector<Eigen::SparseMatrix<double> > Ms; // the gradient operator; Ms[i] * F gives the gradient of F on triangle i
    std::vector<Eigen::Matrix3d> Js; // rotations by 90 degrees

    Eigen::SparseMatrix <double> A; // precomputed matrices used in alternating minimization
    std::vector<Eigen::SparseMatrix<double> > Mbars;
    Eigen::SparseMatrix <double> Mbar;
    
};

struct OptVars
{
    Eigen::VectorXd v; //3F
    Eigen::VectorXd w; //3F
    Eigen::VectorXd D; //9F
};

void alternatingMinimization(const Eigen::MatrixXd &v0, const MeshData &mesh, double lambda, double mu);

#endif
