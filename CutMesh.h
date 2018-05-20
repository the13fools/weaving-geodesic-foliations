#ifndef CUTMESH_H
#define CUTMESH_H

#include <Eigen/Core>
#include <vector>

void findCuts(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
    std::vector<std::vector<int> > &cuts);
    
void cutMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
    // list of cuts, each of which is a list (in order) of vertex indices of one cut.
    // Cuts can be closed loops (in which case the last vertex index should equal the
    // first) or open (in which case the two endpoint vertices should be distinct).
    // Multiple cuts can cross but there may be strange behavior if cuts share endpoint
    // vertices, or are non-edge-disjoint.
    const std::vector<std::vector<int> > &cuts,
    // new vertices and faces
    // **DO NOT ALIAS V OR F!**
    Eigen::MatrixXd &newV,
    Eigen::MatrixXi &newF
);

#endif
