#ifndef COVERMESH_H
#define COVERMESH_H

#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>

class FieldSurface;

class CoverMesh
{
public:

    FieldSurface *fs;

    int nFields_unaugmented;
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
};

#endif
