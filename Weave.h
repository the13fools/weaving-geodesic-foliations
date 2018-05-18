#ifndef WEAVE_H
#define WEAVE_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

class Surface;

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
   
    Surface *surf;

    ////////////////////////
    // Design Variables
    ////////////////////////

    int nFields_;
    int nFields_unaugmented;
    Eigen::VectorXd vectorFields;    // the vector fields, unrolled into a vector:
                                     // first 2m|F| entries: first vector on face 1, second vector on face 1, third vector on face 1, ..., last vector on face m
                                     // next 2m|F| entries: first beta vector on face 1, ...
                                     // next m|F| entries: first alpha on face 1, ...
    std::vector<Eigen::MatrixXi> Ps; // for each edge i, maps indices from triangle E(i,1) to indices in triangle E(i,0), with sign. I.e. the vector on E(i,1) corresponding to vector j on E(i,0) is \sum_k Ps[i](j,k) v(E(i,1),k)
    std::vector<Handle> handles; // handles on the vector fields
 
    std::vector<Cut> cuts; // list of cuts 

    bool addHandle(Handle h);  // this method will add a handle, also taking care to normalize the handle vector length

    int nFields() const { return nFields_; }
    int nHandles() const { return handles.size(); }    

    int vidx(int face, int field) const;
    Eigen::Vector2d v(int face, int field) const;
    int betaidx(int face, int field) const;
    Eigen::Vector2d beta(int face, int field) const;
    int alphaidx(int face, int field) const;
    double alpha(int face, int field) const;

    void normalizeFields(); // make all vectors unit-length
    bool fixFields;  // Do not allow vectors to change in optimization.  
    void createVisualizationEdges(Eigen::MatrixXd &edgePts, Eigen::MatrixXd &edgeVecs, Eigen::MatrixXi &edgeSegs, Eigen::MatrixXd &colors);
    void createVisualizationCuts(Eigen::MatrixXd &cutPts1, Eigen::MatrixXd &cutPts2);
    
    // compute an energy on each face that measures the failure of the vector fields on that face to parallel transport to the
    // equivalent vector on the neighboring faces, using the trivial connection between faces
    // each value will be in the range [0, 3*m*PI].
    void connectionEnergy(Eigen::VectorXd &energies);

    void removePointsFromMesh(std::vector<int> vIds);

    void serialize(const std::string &filename);
    void serialize_forexport(const std::string &filename);
    void deserialize(const std::string &filename);

    std::vector<long> _BFS_adj_list(std::vector<std::vector<long> > & relaxadj_list, int i);
    std::vector<Eigen::MatrixXd> _augmentPs();
    void augmentField();
    void computeFunc(double scalesInit);
    std::vector<double> theta;
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
    void centerAndScale(Eigen::MatrixXd &V);    
};

#endif
