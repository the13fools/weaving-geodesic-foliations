#ifndef RELAXVIEWER_H
#define RELAXVIEWER_H

#include <igl/viewer/Viewer.h>

#include <Eigen/Core>
#include <vector>

#include "FieldOptimization.h"

enum shading_enum {
    OP_ERROR = 0,
    INIT_DIRECTION,
    INIT_MAGNITUDE
};
void updateView(const MeshData *mdata, igl::viewer::Viewer *viewer); 

#endif
