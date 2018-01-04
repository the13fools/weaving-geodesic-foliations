#ifndef FIELDOPTIMIZATION_H
#define FIELDOPTIMIZATION_H

#include <Eigen/Core>
#include "StateStructs.h"

void initOptVars(const MeshData &mesh, const Eigen::MatrixXd &v0, OptVars &vars);

double energy(const OptVars &vars, const MeshData &mesh, Weights &weights, VisualizationState &vs);

void alternatingMinimization(const MeshData &mesh,VisualizationState &vs, Weights &w, OptVars &vars);

#endif
