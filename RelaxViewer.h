#ifndef RELAXVIEWER_H
#define RELAXVIEWER_H

#include <igl/viewer/Viewer.h>

#include <Eigen/Core>
#include <vector>

#include "StateStructs.h"


void computeSelfIntersections(const std::vector<Eigen::Vector3d> &curve, int idx, std::vector<Collision> &collisions);

void traceCurve(const MeshData &md, 
	        const Eigen::Vector3d dir, int faceId, 
		std::vector<Eigen::Vector3d> &curve,
		std::vector<Eigen::Vector3d> &normal);

void updateView(const MeshData *mdata, igl::viewer::Viewer *viewer); 

#endif
