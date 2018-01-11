#include "Permutations.h"
#include "Weave.h"
#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>

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
        int sign = (bestsigns & (1 << i)) ? -1 : 1;        
        P(i, bestperm[i]) = sign;
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

void findSingularVertices(const Weave &weave, std::vector<int> &singularVerts)
{
    int nverts = weave.nVerts();
    singularVerts.clear();

    // in principle this can be done in O(|V|) using circulators which we do not currently compute
    // O(|V||F|) for now

    int nfaces = weave.nFaces();
    int m = weave.nFields();

    for (int i = 0; i < nverts; i++)
    {
        int startface = -1;
        int startspoke = -1;
        bool done = false;
        for (int j = 0; j < nfaces; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                if (weave.F(j, k) == i)
                {
                    startface = j;
                    startspoke = (k+1)%3;
                    done = true;
                    break;
                }
            }
            if (done)
                break;
        }
        assert(startface != -1);

        Eigen::MatrixXi totperm(m, m);
        totperm.setIdentity();

        int curface = startface;
        int curspoke = startspoke;
        while (true)
        {
            int edge = weave.faceEdges(curface, curspoke);
            int side = (weave.E(edge, 0) == curface) ? 0 : 1;
            int nextface = weave.E(edge, 1 - side);
	    if ( side == 0 )
	    {
		totperm *= weave.Ps[edge].transpose();
	    }
	    else 
	    { 
		totperm *= weave.Ps[edge]; 
	    }
	    
            curface = nextface;
            for (int k = 0; k < 3; k++)
            {
                if (weave.F(nextface, k) == i)
                {
                    curspoke = (k + 1) % 3;
                    break;
                }
            }

            if (curface == startface)
                break;
        }

        bool isidentity = true;
        for (int j = 0; j < m; j++)
        {
            if (totperm(j, j) != 1)
                isidentity = false;
        }
        if (!isidentity)
            singularVerts.push_back(i);
    }
}
