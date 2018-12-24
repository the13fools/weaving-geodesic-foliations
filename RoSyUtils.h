#ifndef ROSYUTILS_H
#define ROSYUTILS_H

#include <Eigen/Core>
#include <vector>

class Surface;

// computes the angle between the barycentric u vector of the given face and the given vector, represented in the face's barycentric coordinates.
double vectorAngle(Surface &s, int face, const Eigen::Vector2d &v);

// converts a representative vector (in the barycentric coordinate of the face, and representing a RoSy with respect to the barycentric u vector of the given face)
// into a n-tuple of RoSy vectors. These vectors are in counterclockwise order, but the cyclic permutation is arbitrary.
void repVecToRoSy(Surface &s, int face, const Eigen::Vector2d &v, std::vector<Eigen::Vector2d> &rv, int rosyN);

#endif