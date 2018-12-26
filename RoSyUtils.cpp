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

void repVecToRoSy(Surface &s, int face, const Eigen::Vector2d &v, std::vector<Eigen::Vector2d> &rv, int rosyN)
{
    rv.clear();
    double theta = vectorAngle(s, face, v);
    double basetheta = theta / double(rosyN);
    Eigen::Vector3d n = s.faceNormal(face);
    Eigen::Matrix<double, 3, 2> B = s.data().Bs[face];
    Eigen::Vector3d u = B.col(0);
    Eigen::Matrix2d BTBinv = (B.transpose()*B).inverse();
    for (int j = 0; j < rosyN; j++)
    {
        double angle = basetheta + 2.0 * M_PI * double(j) / double(rosyN);
        Eigen::Vector2d rosyvec = BTBinv * B.transpose() * Eigen::AngleAxisd(angle, n).toRotationMatrix() * u;    
        rv.push_back(rosyvec);
    }
}