#include "RoSyUtils.h"
#include "Surface.h"

#include <Eigen/Geometry>
#include <cmath>

double vectorAngle(Surface &s, int face, const Eigen::Vector2d &v)
{
    Eigen::Vector3d n = s.faceNormal(face);
    Eigen::Vector3d extv = s.data().Bs[face] * v;
    Eigen::Vector3d u = s.data().Bs[face].col(0);
    return 2.0 * std::atan2(u.cross(extv).dot(n), u.norm() * extv.norm() + u.dot(extv));
}

void repVecToRoSy(Surface &s, int face, const Eigen::Vector2d &v, Eigen::Vector2d &rv1, Eigen::Vector2d &rv2, Eigen::Vector2d &rv3)
{
    double theta = vectorAngle(s, face, v);
    double basetheta = theta / 3.0;
    Eigen::Vector3d n = s.faceNormal(face);
    Eigen::Matrix<double, 3, 2> B = s.data().Bs[face];
    Eigen::Vector3d u = B.col(0);
    Eigen::Matrix2d BTBinv = (B.transpose()*B).inverse();
    rv1 = BTBinv * B.transpose() * Eigen::AngleAxisd(basetheta, n).toRotationMatrix() * u;
    rv2 = BTBinv * B.transpose() * Eigen::AngleAxisd(basetheta + 2.0 * M_PI / 3.0, n).toRotationMatrix() * u;
    rv3 = BTBinv * B.transpose() * Eigen::AngleAxisd(basetheta + 4.0 * M_PI / 3.0, n).toRotationMatrix() * u;
}