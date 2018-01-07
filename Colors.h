#ifndef COLORS_H
#define COLORS_H

#include <Eigen/Core>

double hue2rgb(double p, double q, double t);

Eigen::Vector3d hslToRgb(double h, double s, double l);

Eigen::Vector3d heatmap(double val, double valmin, double valmax);

#endif
