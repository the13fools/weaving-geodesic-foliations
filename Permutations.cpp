#include "Permutations.h"
#include "Weave.h"
#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>

double angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2, const Eigen::Vector3d axis)
{
    return 2.0 * atan2(v1.cross(v2).dot(axis), v1.norm() * v2.norm() + v1.dot(v2));
}

void reassignOnePermutation(Weave &weave, int edge, Eigen::MatrixXi &P)
{
    int f = weave.E(edge, 0);
    int g = weave.E(edge, 1);
    Eigen::Vector3d n = weave.faceNormal(f);
    int m = weave.nFields();
    P.resize(m, m);
    P.setZero();
    //gather the vectors
    std::vector<Eigen::Vector3d> fvecs;
    std::vector<Eigen::Vector3d> gvecs;

    Eigen::Matrix2d T = weave.Ts.block<2, 2>(2 * edge, 2);

    for (int i = 0; i < m; i++)
    {
        fvecs.push_back(weave.Bs[f] * weave.v(f, i));
        gvecs.push_back(weave.Bs[f] * T * weave.v(g, i));
    }

    double best = std::numeric_limits<double>::infinity();
    int bestsigns = -1;
    std::vector<int> bestperm;

    // try all permutations
    std::vector<int> perm;
    for (int i = 0; i < m; i++)
        perm.push_back(i);
    do
    {
        // try every sign assignment
        int signmax = 1 << m;
        for (int signs = 0; signs < signmax; signs++)
        {
            // check this permutation, signs pair
            double tottheta = 0;
            for (int i = 0; i < m; i++)
            {
                double sign = (signs & (1 << i)) ? -1.0 : 1.0;
                double theta = angle(fvecs[i], sign*gvecs[perm[i]], n);
                tottheta += theta*theta;
            }
            if (tottheta < best)
            {
                best = tottheta;
                bestsigns = signs;
                bestperm = perm;
            }
        }
    } while (next_permutation(perm.begin(), perm.end()));

    for (int i = 0; i < m; i++)
    {
        P(i, bestperm[i]) = (bestsigns & (1 << i)) ? -1 : 1;
    }
}

int reassignPermutations(Weave &weave)
{
    int nedges = weave.nEdges();
    int count = 0;
    for (int i = 0; i < nedges; i++)
    {
        Eigen::MatrixXi P;
        reassignOnePermutation(weave, i, P);
        if (P != weave.Ps[i])
            count++;
        weave.Ps[i] = P;
    }
    return count;
}