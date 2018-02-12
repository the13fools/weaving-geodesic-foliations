#ifndef WEAVE_H
#define WEAVE_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

// handles on the mesh
struct Handle
{
    int face; // which face the handle is on
    int field; // which of the m fields we're prescribing
    Eigen::Vector2d dir; // the vector itself in triangle barycentric coordinates
};

struct Cut
{ 
    std::vector<std::pair<int, int> > path; // (edge, orientation) list        
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
    std::vector< std::vector<int> > vertEdges; // |V|, F[i] is a list of edges neighboring vertex i

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
    std::vector<Eigen::MatrixXi> Ps; // for each edge i, maps indices from triangle E(i,1) to indices in triangle E(i,0), with sign. I.e. the vector on E(i,1) corresponding to vector j on E(i,0) is \sum_k Ps(j,k) v(E(i,1),k)
    std::vector<Handle> handles; // handles on the vector fields
 
    std::vector<Cut> cuts; // list of cuts 

    bool addHandle(Handle h);  // this method will add a handle, also taking care to normalize the handle vector length

    // Finds shortest (combinatorial) path from start to end vertex. Each path entry is a combination of (1) the edge index along the path, and (2) the orientation: the jth path segment goes from
    // edgeVerts(path[j].first, path[j].second) to edgeVerts(path[j].first, 1 - path[j].second).
    // List will be empty if no path exists (vertices lie on disconnected components).
    void shortestPath(int startVert, int endVert, std::vector<std::pair<int, int> > &path);

    int nVerts() const { return V.rows(); }
    int nFaces() const { return F.rows(); }
    int nEdges() const { return E.rows(); }
    int nFields() const { return nFields_; }
    int nHandles() const { return handles.size(); }
    int numInteriorEdges() const;

    int vidx(int face, int field) const;
    Eigen::Vector2d v(int face, int field) const;
    int betaidx(int face, int field) const;
    Eigen::Vector2d beta(int face, int field) const;
    int alphaidx(int face, int field) const;
    double alpha(int face, int field) const;

    Eigen::Vector3d faceNormal(int face) const;
    void normalizeFields(); // make all vectors unit-length
    void createVisualizationEdges(Eigen::MatrixXd &edgePts, Eigen::MatrixXd &edgeVecs, Eigen::MatrixXi &edgeSegs, Eigen::MatrixXd &colors);
    void createVisualizationCuts(Eigen::MatrixXd &cutPts1, Eigen::MatrixXd &cutPts2);

    void removePointsFromMesh(std::vector<int> vIds);

    void serialize(const std::string &filename);
    void serialize_forexport(const std::string &filename);
    void deserialize(const std::string &filename);

    std::vector<long> _BFS_adj_list(std::vector<std::vector<long> > adj_list, int i);
    std::vector<Eigen::MatrixXd> _augmentPs();
    void augmentField();
    void computeFunc(double scalesInit);
    std::vector<double> theta;
    std::vector<double> simpleKron(Eigen::Vector3d vec, int augRow);
    std::vector<double> simpleKron(std::vector<double> vec, int augRow);
    Eigen::SparseMatrix<double> faceLaplacian();
    double barycentric(double val1, double val2, double target);
    bool crosses(double isoval, double val1, double val2, double minval, 
        double maxval, double &bary);
    int extractIsoline(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, 
        const Eigen::MatrixXi &faceNeighbors, const Eigen::VectorXd &func, 
        double isoval, double minval, double maxval);
    void drawISOLines(int numISOLines);
    std::vector<std::vector<Eigen::Vector3d> > isoLines;
    std::vector<std::vector<Eigen::Vector3d> > isoNormal;
    bool augmented;
    Eigen::MatrixXd VAugmented; // |V| x 3
    Eigen::MatrixXi FAugmented; // |F| x 3
    Eigen::MatrixXd VBuffer; // |V| x 3
    Eigen::MatrixXi FBuffer; // |F| x 3

private:
    // scale mesh to unit size
    void centerAndScale();
    // computes E and faceedges/faceWings from V and F
    void buildConnectivityStructures();
    // compute the geometric data structs
    void buildGeometricStructures();
};

#endif
