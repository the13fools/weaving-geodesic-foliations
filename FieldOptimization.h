#ifndef FIELDOPTIMIZATION_H
#define FIELDOPTIMIZATION_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>


struct OptVars
{
    Eigen::VectorXd v; //3F
    Eigen::VectorXd w; //3F
    Eigen::VectorXd D; //9F

    Eigen::MatrixXd V_opt;
    Eigen::MatrixXd W_opt;
};

// Keep application state in here until it gets annoying. 
// All state that varies with optimization goes into it's own structs for ease of refactoring later
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
    Eigen::SparseMatrix <double> Mbar;

    Eigen::MatrixXd centroids_F;
    Eigen::MatrixXd v0; // v at init, for visualizing change in descent and computing energy    

    OptVars optVars;
};

void initOptVars(const Eigen::MatrixXd &v0, const std::vector<Eigen::SparseMatrix<double> > Ms, OptVars &vars);


void alternatingMinimization(const MeshData &mesh, double lambda, double mu, OptVars &vars);

#endif
