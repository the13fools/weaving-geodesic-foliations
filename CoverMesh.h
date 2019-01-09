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
class LocalFieldIntegration;
class GlobalFieldIntegration;

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
    CoverMesh(const Surface &originalSurf, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXi &oldToNewVertMap, const Eigen::MatrixXd &field, int ncovers);
    ~CoverMesh();

    CoverMesh(const CoverMesh &) = delete;
    CoverMesh &operator=(const CoverMesh &) = delete;

    FieldSurface *fs;
    Eigen::VectorXd theta;
    Eigen::VectorXd scales;
    
    void createVisualization(Eigen::MatrixXd &V, Eigen::MatrixXi &F, 
        Eigen::MatrixXd &edgePts, Eigen::MatrixXi &edgeSegs, Eigen::MatrixXd &colors, 
        Eigen::MatrixXd &cutPts1, Eigen::MatrixXd &cutPts2, Eigen::MatrixXd &cutColors,
        bool hideVectors,
        double vectorScale);

    void integrateField(LocalFieldIntegration *lmethod, GlobalFieldIntegration *gmethod, double globalScale);
    void roundAntipodalCovers(int numISOLines);
    double renderScale() {return renderScale_;}
    const Surface &splitMesh() const;
    void gradThetaDeviation(Eigen::VectorXd &error) const;
    
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
    Surface *originalSurf_;
    double renderScale_;    
    // edges along which the multiple cover is cut to create a topological disk
    std::vector<int> slicedEdges;
};

#endif
