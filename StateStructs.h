#ifndef STATESTRUCTS_H
#define STATESTRUCTS_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

struct OptVars
{
    Eigen::VectorXd vbar; //2F
    Eigen::VectorXd wbar; //2F
    Eigen::VectorXd D; //9F

    Eigen::MatrixXd V_opt;
    Eigen::MatrixXd W_opt;
};

struct Weights
{
    double lambdaGeodesic; // weight of the geodesic constraints
    double lambdaVW; // weight of v = w constraint
    double lambdaVD; // weight of v, D compatibility constraint
    double lambdaDreg; // weight of regularization term on D
    double lambdaunit; // weight of v and w unit-length
    Eigen::VectorXd handleWeights; // one weight per face, 1.0 = use the input v0 on this face as a handle, 0.0 = ignore this input v0.
};


enum shading_enum {
    OP_ERROR = 0,
    INIT_DIRECTION,
    INIT_MAGNITUDE
};

struct VisualizationState
{
    bool normFaceVectors;
    std::vector<Eigen::Vector3d> curve;
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
    Eigen::MatrixXi faceWings; // the three vertices opposite the edges of each face
    std::vector<Eigen::Matrix3d> Js; // rotations by 90 degrees
    Eigen::MatrixXd Ms; // discrete gradient operator vectors
    Eigen::MatrixXd B; // 3*faces x 2 matrix of barycentric bases

    Eigen::MatrixXd centroids_F;
    Eigen::SparseMatrix<double> H; // precomputed matrices for compatibility Hessian
    Eigen::SparseMatrix<double> C;
    Eigen::MatrixXd v0; // v at init, for visualizing change in descent and computing energy
    Eigen::SparseMatrix<double> Bmat; // maps a 2xF vector of vectors in barycentric coordinates to a 3xF vector of extrinsic tangent vectors

    OptVars optVars;
    VisualizationState vs;
};


#endif
