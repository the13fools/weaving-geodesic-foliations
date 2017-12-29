#ifndef RELAXVIEWER_H
#define RELAXVIEWER_H

#include <Eigen/Core>
#include "StateStructs.h"
#include <igl/viewer/Viewer.h>

void updateView(const MeshData &data, const OptState &state, igl::viewer::Viewer *viewer);


#endif
