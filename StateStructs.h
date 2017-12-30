#ifndef STATESTRUCTS_H
#define STATESTRUCTS_H

#include <Eigen/Core>

// stuff that is precomputed about the mesh and doesn't change during optimization
struct MeshData
{
    Eigen::MatrixXd V; // vertex positions
    Eigen::MatrixXi F; // face indices
    Eigen::MatrixXi E; // edge indices
    Eigen::MatrixXi F_edges; // face edge indices
    std::vector<Eigen::SparseMatrix<double> > Ms; // the gradient operator; Ms[i] * F gives the gradient of F on triangle i

    Eigen::MatrixXd N; // face normals TODO
    std::vector<Eigen::SparseMatrix<double> > Js; // the perp operator; Js[i] * w_i rotates w_i by pi/2 in the plane of face i
    
    Eigen::MatrixXd W_init; // Initial Vector Field, for visualization 
    Eigen::MatrixXd centroids_F; // A list of centroids per face
};

// stuff that is updated through the course of optimization.  Mostly a convenience to make visualization easier.
// Still good to unpack for most functions to maintain clarity + compile time checks
struct OptState
{
    Eigen::MatrixXd dot; // This is the "test" field
    std::vector<Eigen::MatrixXd> Ds;  // These are the operators with fields built in
    Eigen::MatrixXd W;
};

enum shading_enum {
    OP_ERROR = 0,
    INIT_DIRECTION,
    INIT_MAGNITUDE
};

#endif
