#ifndef FIELDOPTIMIZATION_H
#define FIELDOPTIMIZATION_H

#include <Eigen/Core>
#include "StateStructs.h"

void initOptVars(const MeshData &mesh, const Eigen::MatrixXd &v0, OptVars &vars);


void alternatingMinimization(const MeshData &mesh, Weights &w, OptVars &vars);

#endif
