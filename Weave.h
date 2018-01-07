#ifndef WEAVE_H
#define WEAVE_H

#include <Eigen/Core>
#include <vector>

// handles on the mesh
struct Handle
{
    int face; // which face the handle is on
    int field; // which of the m fields we're prescribing
    Eigen::Vector2d dir; // the vector itself in triangle barycentric coordinates
};

class Weave
{
public:
    Weave(const std::string &objname, int m);
    ~Weave();

    ////////////////////////
    // Precomputed Stuff
    ////////////////////////

    // Combinatorics Data Structures

    Eigen::MatrixXd V; // |V| x 3
    Eigen::MatrixXi F; // |F| x 3
    Eigen::MatrixXi E; // |E| x 2
    Eigen::MatrixXi edgeVerts; // |E| x 2
    Eigen::MatrixXi faceEdges; // |F| x 3, F(i,j) is the edge opposite vertex j in triangle i
    Eigen::MatrixXi faceNeighbors; // |F| x 3, F(i,j) is the face opposite vertex j in triangle i
    Eigen::MatrixXi faceWings; // |F| x 3, F(i,j) is vertex opposite vertex j in triangle i

    // Geometric Data Structures

    Eigen::MatrixXd cDiffs; // cDiffs.row(2*i) is the vector in barycentric coordinates of face E(i,0) from centroid on face E(i,0) to that on face E(i,1) after rotating triangle E(i,1) to E(i,0)'s tangent plane
                            // cDiffs.row(2*i+1) is the same but with the roles of the faces reversed
    std::vector<Eigen::Matrix<double, 3, 2> > Bs; // Basis vectors in ambient coordinates for the barycentric coordinates on faces
    Eigen::MatrixXd Ts;     // Transition matrices. Ts.block<2,2>(2*i,0) maps vectors from barycentric coordinates of face E(i,0) to barycentric coordinates of E(i,1). Ts.block<2,2>(2*i, 2) is opposite.
    Eigen::MatrixXd Js;     // Js.block<2,2>(2*i,0) rotates vectors on face i (in face i's barycentric coordinates) to the perpendicular vector (as measured in ambient space)

    double averageEdgeLength; // exactly what it says on the tin

    ////////////////////////
    // Design Variables
    ////////////////////////

    int nFields_;
    Eigen::VectorXd vectorFields;    // the vector fields, unrolled into a vector:
                                     // first 2m|F| entries: first vector on face 1, second vector on face 1, third vector on face 1, ..., last vector on face m
                                     // next 2m|F| entries: first beta vector on face 1, ...
                                     // next m|F| entries: first alpha on face 1, ...
    std::vector<Eigen::MatrixXi> Ps; // for each edge i, maps indices from triangle E(i,0) to indices in triangle E(i,1), with sign
    std::vector<Handle> handles; // handles on the vector fields

    bool addHandle(Handle h);  // this method will add a handle, also taking care to normalize the handle vector length

    int nVerts() const { return V.rows(); }
    int nFaces() const { return F.rows(); }
    int nEdges() const { return E.rows(); }
    int nFields() const { return nFields_; }
    int nHandles() const { return handles.size(); }

    int vidx(int face, int field) const;
    Eigen::Vector2d v(int face, int field) const;
    int betaidx(int face, int field) const;
    Eigen::Vector2d beta(int face, int field) const;
    int alphaidx(int face, int field) const;
    double alpha(int face, int field) const;

    Eigen::Vector3d faceNormal(int face);

    void createVisualizationEdges(Eigen::MatrixXd &edgePts, Eigen::MatrixXd &edgeVecs, Eigen::MatrixXi &edgeSegs, Eigen::MatrixXd &colors);

private:
    // computes E and faceedges/faceWings from V and F
    void buildConnectivityStructures();
    // compute the geometric data structs
    void buildGeometricStructures();
};

#endif