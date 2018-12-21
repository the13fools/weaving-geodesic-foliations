#ifndef COVERMESH_H
#define COVERMESH_H

#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "Traces.h"

class FieldSurface;
class Weave;
class Surface;
class TraceSet;
class FieldIntegration;

struct CoverData
{
    Surface *splitMesh; // the covering mesh split into 2*m copies of the original mesh
    
    std::vector<Eigen::Vector3d> splitOffsets; // translation of each split mesh
    std::map<int, std::vector<int> > coverToSplitVerts; // map from vertex indices on the covering mesh to their "child" vertices on the split mesh
    Eigen::VectorXi splitToCoverVerts; // inverse of coverToSplitVerts

    std::vector<int> splitMeshCuts; // edges of the split mesh that are cuts
};

class CoverMesh
{
public:
    CoverMesh(const Weave &parent, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXi &oldToNewVertMap, const Eigen::MatrixXd &field, int ncovers);
    ~CoverMesh();

    FieldSurface *fs;
    Eigen::VectorXd theta;
    
    void createVisualization(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &edgePts, Eigen::MatrixXd &edgeVecs, Eigen::MatrixXi &edgeSegs, Eigen::MatrixXd &colors, 
        Eigen::MatrixXd &cutPts1, Eigen::MatrixXd &cutPts2, Eigen::MatrixXd &cutColors);

    void integrateField(FieldIntegration *method);
    double renderScale() {return renderScale_;}
    const Surface &splitMesh() const;
    
    // maps indices of vertices on the visualization mesh to corresponding "parent" vertices on the cover mesh
    int visMeshToCoverMesh(int vertid);
    void recomputeIsolines(int numISOLines, std::vector<Trace> &isotraces);
    void drawTraceOnSplitMesh(const Trace &trace, Eigen::MatrixXd &pathStarts, Eigen::MatrixXd &pathEnds) const;    
   
private:
    double inversePowerIteration(Eigen::SparseMatrix<double> &M, Eigen::VectorXd &evec, int iters);
    void initializeSplitMesh(const Eigen::VectorXi &oldToNewVertMap);

    double barycentric(double val1, double val2, double target);
    bool crosses(double isoval, double val1, double val2, double minval, 
        double maxval, double &bary);
    void extractIsoline(const Eigen::VectorXd &func, 
        double isoval, double minval, double maxval, std::vector<Trace> &isotrace);
    
    std::vector<std::vector<Eigen::Vector3d> > isoNormal;

    CoverData data_;
    int ncovers_;
    const Weave &parent_;
    double renderScale_;    
    // edges along which the multiple cover is cut to create a topological disk
    std::vector<int> slicedEdges;
};

#endif
