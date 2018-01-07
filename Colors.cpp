#include "Colors.h"

double hue2rgb(double p, double q, double t)
{
    if (t < 0) t += 1.0;
    if (t > 1.0) t -= 1.0;
    if (t < 1.0 / 6.0) return p + (q - p) * 6 * t;
    if (t < 0.5) return q;
    if (t < 2.0 / 3.0) return p + (q - p)*(2.0 / 3.0 - t) * 6;
    return p;
}

Eigen::Vector3d hslToRgb(double h, double s, double l)
{
    Eigen::Vector3d result(l, l, l);
    if (s == 0)
        return result;
    double q = (l < 0.5 ? l*(1.0 + s) : l + s - l*s);
    double p = 2.0*l - q;
    result[0] = hue2rgb(p, q, h + 1.0 / 3.0);
    result[1] = hue2rgb(p, q, h);
    result[2] = hue2rgb(p, q, h - 1.0 / 3.0);
    return result;
}

Eigen::Vector3d heatmap(double val, double valmin, double valmax)
{
	double valn = std::max(std::min((val-valmin) / (valmax-valmin), 1.0), 0.0);
	double h = (1.0 - valn) * 240.0 / 360.0;
	double s = 100.0 / 100.0;
	double l = 50.0 / 100.0;
	return hslToRgb(h, s, l);
}

