#ifndef RELAXVIEWER_H
#define RELAXVIEWER_H

#include <igl/viewer/Viewer.h>

#include <Eigen/Core>
#include <vector>

#include "StateStructs.h"

void updateView(const MeshData *mdata, igl::viewer::Viewer *viewer); 

#endif
