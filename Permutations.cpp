#include "Permutations.h"
#include "Weave.h"
#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
#include "Surface.h"

double angle(const Eigen::Vector3d &v1, const Eigen::Vector3d &v2, const Eigen::Vector3d axis)
{
    return 2.0 * atan2(v1.cross(v2).dot(axis), v1.norm() * v2.norm() + v1.dot(v2));
}

void reassignOneCutPermutation(Weave &weave, int cut, Eigen::MatrixXi &P)
{
    int m = weave.fs->nFields();
    double best = std::numeric_limits<double>::infinity();
    int bestsigns = -1;
    std::vector<int> bestperm;

    int nedges = weave.cuts[cut].path.size();

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

            for (int i = 0; i < nedges; i++)
            {
                std::vector<Eigen::Vector3d> fvecs;
                std::vector<Eigen::Vector3d> gvecs;
                int orient = weave.cuts[cut].path[i].second;
                int f = weave.fs->data().E(weave.cuts[cut].path[i].first, orient);
                int g = weave.fs->data().E(weave.cuts[cut].path[i].first, 1 - orient);
                if (f == -1 || g == -1)
                    continue;

                Eigen::Matrix2d T = weave.fs->data().Ts.block<2, 2>(2 * weave.cuts[cut].path[i].first, 2 - 2 * orient);
                Eigen::Vector3d n = weave.fs->faceNormal(f);

                for (int j = 0; j < m; j++)
                {
                    fvecs.push_back(weave.fs->data().Bs[f] * weave.fs->v(f, j));
                    gvecs.push_back(weave.fs->data().Bs[f] * T * weave.fs->v(g, j));
                }

                for (int j = 0; j < m; j++)
                {
                    double sign = (signs & (1 << j)) ? -1.0 : 1.0;
                    double theta = angle(fvecs[j], sign*gvecs[perm[j]], n);
                    tottheta += theta*theta;
                }
            }
            if (tottheta < best)
            {
                best = tottheta;
                bestsigns = signs;
                bestperm = perm;
            }
        }
    } while (next_permutation(perm.begin(), perm.end()));

    P.resize(m, m);
    P.setZero();
    for (int i = 0; i < m; i++)
    {
        int sign = (bestsigns & (1 << i)) ? -1 : 1;        
        P(i, bestperm[i]) = sign;
    }
}

void reassignOnePermutation(Weave &weave, int edge, Eigen::MatrixXi &P)
{
    int f = weave.fs->data().E(edge, 0);
    int g = weave.fs->data().E(edge, 1);
    Eigen::Vector3d n = weave.fs->faceNormal(f);
    int m = weave.fs->nFields();
    P.resize(m, m);
    P.setZero();
    //gather the vectors
    std::vector<Eigen::Vector3d> fvecs;
    std::vector<Eigen::Vector3d> gvecs;

    Eigen::Matrix2d T = weave.fs->data().Ts.block<2, 2>(2 * edge, 2);

    for (int i = 0; i < m; i++)
    {
        fvecs.push_back(weave.fs->data().Bs[f] * weave.fs->v(f, i));
        gvecs.push_back(weave.fs->data().Bs[f] * T * weave.fs->v(g, i));
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
    int nedges = weave.fs->nEdges();
    int count = 0;
    for (int i = 0; i < nedges; i++)
    {
        Eigen::MatrixXi P;
        reassignOnePermutation(weave, i, P);
        if (P != weave.fs->Ps(i))
            count++;
        weave.fs->Ps_[i] = P;
    }
    return count;
}

void findSingularVertices(const Weave &weave, std::vector<int> &topologicalSingularVerts, std::vector<std::pair<int, int> > &geometricSingularVerts)
{
    int nverts = weave.fs->nVerts();
    topologicalSingularVerts.clear();
    geometricSingularVerts.clear();

    // in principle this can be done in O(|V|) using circulators which we do not currently compute
    // O(|V||F|) for now

    int nfaces = weave.fs->nFaces();
    int m = weave.fs->nFields();

    for (int i = 0; i < nverts; i++)
    {
        int startface = -1;
        int startspoke = -1;
        bool done = false;
        for (int j = 0; j < nfaces; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                if (weave.fs->data().F(j, k) == i)
                {
                    startface = j;
                    startspoke = (k + 1) % 3;
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
        Eigen::VectorXd angles(m);
        angles.setZero();

        int curface = startface;
        int curspoke = startspoke;
        double totangle = 0;

        bool isboundary = false;
        int iter = 0;

        while (true)
        {
            int edge = weave.fs->data().faceEdges(curface, curspoke);
            int side = (weave.fs->data().E(edge, 0) == curface) ? 0 : 1;
            int nextface = weave.fs->data().E(edge, 1 - side);
            if (nextface == -1 || curface == -1)
            {
                isboundary = true;
                break;
            }
            
            Eigen::MatrixXi nextperm;
            if (side == 0)
            {
                nextperm = weave.fs->Ps(edge).transpose();
            }
            else
            {
                nextperm = weave.fs->Ps(edge);
            }
            
            Eigen::Vector3d normal = weave.fs->faceNormal(curface);

            for (int j = 0; j < m; j++)
            {
                Eigen::Vector3d curv = weave.fs->data().Bs[curface] * weave.fs->v(curface, j);
                Eigen::Vector2d nextvbary(0,0);
                for(int k=0; k<m; k++)
                {
                    nextvbary += nextperm(k,j) * weave.fs->v(nextface, k);
                }
                Eigen::Vector3d nextv = weave.fs->data().Bs[curface] * weave.fs->data().Ts.block<2, 2>(2 * edge, 2 - 2 * side) * nextvbary;
                angles[j] += angle(curv, nextv, normal);
            }
            
            Eigen::VectorXd nextangles(m);
            nextangles.setZero();
            for(int j=0; j<m; j++)
            {
                for(int k=0; k<m; k++)
                {
                    if(nextperm(j,k) != 0)
                    {
                        nextangles[j] = angles[k];
                    }
                }
            }
            
            angles = nextangles;
            
            totperm *= nextperm;

            int spokep1 = (curspoke + 1) % 3;
            int apex = (curspoke + 2) % 3;
            Eigen::Vector3d v1 = weave.fs->data().V.row(weave.fs->data().F(curface,curspoke)) - weave.fs->data().V.row(weave.fs->data().F(curface,apex));
            Eigen::Vector3d v2 = weave.fs->data().V.row(weave.fs->data().F(curface,spokep1)) - weave.fs->data().V.row(weave.fs->data().F(curface,apex));
            totangle += angle(v1, v2, normal);

            curface = nextface;
            for (int k = 0; k < 3; k++)
            {
                if (weave.fs->data().F(nextface, k) == i)
                {
                    curspoke = (k + 1) % 3;
                    break;
                }
            }

            if (curface == startface)
                break;
        }

        if (!isboundary)
        {
            bool isidentity = true;
            for (int j = 0; j < m; j++)
            {
                if (totperm(j, j) != 1)
                    isidentity = false;
            }
            if (!isidentity)
                topologicalSingularVerts.push_back(i);
            else
            {
                for (int j = 0; j < m; j++)
                {
                    const double PI = 3.1415926535898;
                    double index = angles[j] + 2 * PI - totangle;
                    if (fabs(index) > PI)
                    {
                        geometricSingularVerts.push_back(std::pair<int, int>(i, j));
                    }
                }
            }
        }
    }
}

int reassignCutPermutations(Weave &weave)
{
    int ncuts = (int)weave.cuts.size();
    int tot = 0;
    for (int i = 0; i < ncuts; i++)
    {
        Eigen::MatrixXi P;
        reassignOneCutPermutation(weave, i, P);
        for (int j = 0; j < weave.cuts[i].path.size(); j++)
        {
            Eigen::MatrixXi Pedge = P;
            if (weave.cuts[i].path[j].second == 1)
                Pedge.transposeInPlace();
            if (weave.fs->Ps(weave.cuts[i].path[j].first) != Pedge)
            {
                weave.fs->Ps_[weave.cuts[i].path[j].first] = Pedge;
                tot++;
            }
        }
    }
    return tot;
}
