#ifndef COVERMESH_H
#define COVERMESH_H

#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>

class FieldSurface;

class CoverMesh
{
public:
    CoverMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &field, int ncovers);
    ~CoverMesh();

    FieldSurface *fs;
    Eigen::VectorXd theta;

    void createVisualizationEdges(Eigen::MatrixXd &edgePts, Eigen::MatrixXd &edgeVecs, Eigen::MatrixXi &edgeSegs, Eigen::MatrixXd &colors);

    void computeFunc(double scalesInit);
    
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
    Eigen::MatrixXd VBuffer; // |V| x 3
    Eigen::MatrixXi FBuffer; // |F| x 3
};

#endif
