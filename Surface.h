#ifndef SURFACE_H
#define SURFACE_H

#include <Eigen/Core>
#include <vector>

struct SurfaceData
{
    // Combinatorics Data Structures

    Eigen::MatrixXd V; // |V| x 3
    Eigen::MatrixXi F; // |F| x 3
    Eigen::MatrixXi E; // |E| x 2
    Eigen::MatrixXi edgeVerts; // |E| x 2
    Eigen::MatrixXi faceEdges; // |F| x 3, F(i,j) is the edge opposite vertex j in triangle i
    Eigen::MatrixXi faceNeighbors; // |F| x 3, F(i,j) is the face opposite vertex j in triangle i
    Eigen::MatrixXi faceWings; // |F| x 3, F(i,j) is vertex opposite vertex j in triangle i
    std::vector< std::vector<int> > vertEdges; // |V|, F[i] is a list of edges neighboring vertex i

    // Geometric Data Structures

    Eigen::MatrixXd cDiffs; // cDiffs.row(2*i) is the vector in barycentric coordinates of face E(i,0) from centroid on face E(i,0) to that on face E(i,1) after rotating triangle E(i,1) to E(i,0)'s tangent plane
                            // cDiffs.row(2*i+1) is the same but with the roles of the faces reversed
    std::vector<Eigen::Matrix<double, 3, 2> > Bs; // Basis vectors in ambient coordinates for the barycentric coordinates on faces
    Eigen::MatrixXd Ts;     // Transition matrices. Ts.block<2,2>(2*i,0) maps vectors from barycentric coordinates of face E(i,0) to barycentric coordinates of E(i,1). Ts.block<2,2>(2*i, 2) is opposite.
    Eigen::MatrixXd Js;     // Js.block<2,2>(2*i,0) rotates vectors on face i (in face i's barycentric coordinates) to the perpendicular vector (as measured in ambient space)
    double averageEdgeLength; // exactly what it says on the tin
};

// Class wrapping a triangle mesh surface embedded in R^3, along with its combinatorial and geometric data structures
class Surface
{
public:
    Surface(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);

    const SurfaceData &data() { return data_; }

    int nVerts() const { return data_.V.rows(); }
    int nFaces() const { return data_.F.rows(); }
    int nEdges() const { return data_.E.rows(); }
    int numInteriorEdges() const;

    Eigen::Vector3d faceNormal(int face) const;

    // Finds shortest (combinatorial) path from start to end vertex. Each path entry is a combination of (1) the edge index along the path, and (2) the orientation: the jth path segment goes from
    // edgeVerts(path[j].first, path[j].second) to edgeVerts(path[j].first, 1 - path[j].second).
    // List will be empty if no path exists (vertices lie on disconnected components).
    void shortestPath(int startVert, int endVert, std::vector<std::pair<int, int> > &path) const;

private:
    // computes E and faceedges/faceWings from V and F
    void buildConnectivityStructures();
    // compute the geometric data structs
    void buildGeometricStructures();

    SurfaceData data_;
};

#endif
