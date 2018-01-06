#include "GaussNewton.h"
#include "Weave.h"
#include <iostream>
#include <fstream>

using namespace Eigen;

double energy(const Weave &weave, SolverParams params)
{
    double result = 0;
    int nhandles = weave.nHandles();
    for(int i=0; i<nhandles; i++)
    { 
        int field = weave.handles[i].field;
        int face = weave.handles[i].face;
        Vector2d dir = weave.handles[i].dir;
        result += 0.5 * (weave.v(face, field)- dir).squaredNorm();
    }

    // compatibility constraint
    int nedges = weave.nEdges();
    int m = weave.nFields();
    for (int e = 0; e < nedges; e++)
    {
        for (int i = 0; i < m; i++)
        {
            for (int side = 0; side < 2; side++)
            {
                int f = (side == 0 ? weave.E(e, 0) : weave.E(e, 1));
                int g = (side == 0 ? weave.E(e, 1) : weave.E(e, 0));
                Eigen::Matrix2d Jf = weave.Js.block<2, 2>(2 * f, 0);
                Eigen::Vector2d vif = weave.v(f, i);
                Eigen::Matrix2d Dif = weave.beta(f, i) * (Jf*vif).transpose() + weave.alpha(f,i) * vif * vif.transpose();
                Eigen::Vector2d cdiff = weave.cDiffs.row(2 * e + side);
                Eigen::Vector2d vpermut(0, 0);
                Eigen::MatrixXd permut = weave.Ps[e];
                if (side == 1)
                    permut.transposeInPlace();
                for (int field = 0; field < m; field++)
                {
                    vpermut += permut(i, field) * weave.v(g, field);
                }
                Eigen::Matrix2d Tgf = weave.Ts.block<2, 2>(2 * e, 2 - 2 * side);
                result += 0.5 * params.lambdacompat * (Dif*cdiff - (Tgf*vpermut - vif)).squaredNorm();
            }
        }
    }
    return result;
}

void GNEnergy(const Weave &weave, SolverParams params, Eigen::VectorXd &E)
{
    int nhandles = weave.nHandles();
    int nedges = weave.nEdges();
    int m = weave.nFields();

    int nterms = 2 * nhandles + 4 * nedges*m;
    E.resize(nterms);
    E.setZero();

    for(int i=0; i<nhandles; i++)
    { 
        int field = weave.handles[i].field;
        int face = weave.handles[i].face;
        Vector2d dir = weave.handles[i].dir;
        for (int k = 0; k < 2; k++)
            E[2 * i + k] = (weave.v(face, field) - dir)[k];
    }

    int term = 2*nhandles;

    // compatibility constraint
    for (int e = 0; e < nedges; e++)
    {
        for (int i = 0; i < m; i++)
        {
            for (int side = 0; side < 2; side++)
            {
                int f = (side == 0 ? weave.E(e, 0) : weave.E(e, 1));
                int g = (side == 0 ? weave.E(e, 1) : weave.E(e, 0));
                Eigen::Matrix2d Jf = weave.Js.block<2, 2>(2 * f, 0);
                Eigen::Vector2d vif = weave.v(f, i);
                Eigen::Matrix2d Dif = weave.beta(f, i) * (Jf*vif).transpose() + weave.alpha(f,i) * vif * vif.transpose();
                Eigen::Vector2d cdiff = weave.cDiffs.row(2 * e + side);
                Eigen::Vector2d vpermut(0, 0);
                Eigen::MatrixXd permut = weave.Ps[e];
                if (side == 1)
                    permut.transposeInPlace();
                for (int field = 0; field < m; field++)
                {
                    vpermut += permut(i, field) * weave.v(g, field);
                }
                Eigen::Matrix2d Tgf = weave.Ts.block<2, 2>(2 * e, 2 - 2 * side);
                for (int k = 0; k < 2; k++)
                {
                    E[term] = sqrt(params.lambdacompat) * (Dif*cdiff - (Tgf*vpermut - vif))[k];
                    term++;
                }
            }
        }
    }    
}

void trueGradient(const Weave &weave, SolverParams params, Eigen::VectorXd &dE)
{
    dE.resize(weave.vectorFields.size());
    dE.setZero();

    int nhandles = weave.nHandles();
    for(int i=0; i<nhandles; i++)
    { 
        int field = weave.handles[i].field;
        int face = weave.handles[i].face;
        Vector2d dir = weave.handles[i].dir;
        dE.segment<2>(weave.vidx(face, field)) += (weave.v(face, field) - dir);
    }

    // compatibility constraint
    int nedges = weave.nEdges();
    int m = weave.nFields();
    for (int e = 0; e < nedges; e++)
    {
        for (int i = 0; i < m; i++)
        {
            for (int side = 0; side < 2; side++)
            {
                int f = (side == 0 ? weave.E(e, 0) : weave.E(e, 1));
                int g = (side == 0 ? weave.E(e, 1) : weave.E(e, 0));
                Eigen::Matrix2d Jf = weave.Js.block<2, 2>(2 * f, 0);
                Eigen::Vector2d vif = weave.v(f, i);
                Eigen::Matrix2d Dif = weave.beta(f, i) * (Jf*vif).transpose() + weave.alpha(f,i) * vif * vif.transpose();
                Eigen::Vector2d cdiff = weave.cDiffs.row(2 * e + side);
                Eigen::Vector2d vpermut(0, 0);
                Eigen::MatrixXd permut = weave.Ps[e];
                if (side == 1)
                    permut.transposeInPlace();
                for (int field = 0; field < m; field++)
                {
                    vpermut += permut(i, field) * weave.v(g, field);
                }
                Eigen::Matrix2d Tgf = weave.Ts.block<2, 2>(2 * e, 2 - 2 * side);
                
                Eigen::Vector2d innervec = params.lambdacompat*(Dif*cdiff - (Tgf*vpermut - vif));
                dE.segment<2>(weave.vidx(f, i)) += Jf.transpose()*cdiff * weave.beta(f, i).dot(innervec);
                dE.segment<2>(weave.vidx(f, i)) += cdiff.dot(vif) * weave.alpha(f, i) * innervec;
                dE.segment<2>(weave.vidx(f, i)) += cdiff * weave.alpha(f, i) * vif.dot(innervec);
                dE.segment<2>(weave.vidx(f, i)) += innervec;
                for (int field = 0; field < m; field++)
                {
                    dE.segment<2>(weave.vidx(g, field)) -= permut(i, field) * Tgf.transpose() * innervec;
                }
                dE[weave.alphaidx(f, i)] += innervec.dot(vif) * vif.dot(cdiff);
                dE.segment<2>(weave.betaidx(f, i)) += (Jf * vif).dot(cdiff) * innervec;
            }
        }
    }
}

void GNGradient(const Weave &weave, SolverParams params, Eigen::SparseMatrix<double> &J)
{
    int nhandles = weave.nHandles();
    int nedges = weave.nEdges();
    int m = weave.nFields();

    int nterms = 2 * nhandles + 4 * nedges*m;
    J.resize(nterms, weave.vectorFields.size());

    std::vector<Eigen::Triplet<double> > Jcoeffs;

    for(int i=0; i<nhandles; i++)
    { 
        int field = weave.handles[i].field;
        int face = weave.handles[i].face;
        for(int k=0; k<2; k++)
            Jcoeffs.push_back(Triplet<double>(2*i+k, weave.vidx(face, field) + k, 1.0));        
    }

    int term = 2 * nhandles;

    // compatibility constraint
    for (int e = 0; e < nedges; e++)
    {
        for (int i = 0; i < m; i++)
        {
            for (int side = 0; side < 2; side++)
            {
                int f = (side == 0 ? weave.E(e, 0) : weave.E(e, 1));
                int g = (side == 0 ? weave.E(e, 1) : weave.E(e, 0));
                Eigen::Matrix2d Jf = weave.Js.block<2, 2>(2 * f, 0);
                Eigen::Vector2d vif = weave.v(f, i);
                Eigen::Matrix2d Dif = weave.beta(f, i) * (Jf*vif).transpose() + weave.alpha(f,i) * vif * vif.transpose();
                Eigen::Vector2d cdiff = weave.cDiffs.row(2 * e + side);
                Eigen::Vector2d vpermut(0, 0);
                Eigen::MatrixXd permut = weave.Ps[e];
                if (side == 1)
                    permut.transposeInPlace();
                for (int field = 0; field < m; field++)
                {
                    vpermut += permut(i, field) * weave.v(g, field);
                }
                Eigen::Matrix2d Tgf = weave.Ts.block<2, 2>(2 * e, 2 - 2 * side);

                for (int coeff = 0; coeff < 2; coeff++)
                {
                    Eigen::Vector2d innervec(0, 0);
                    innervec[coeff] = sqrt(params.lambdacompat);
                    Vector2d dE = Jf.transpose()*cdiff * weave.beta(f, i).dot(innervec);
                    dE += cdiff.dot(vif) * weave.alpha(f, i) * innervec;
                    dE += cdiff * weave.alpha(f, i) * vif.dot(innervec);
                    dE += innervec;
                    for (int k = 0; k < 2; k++)
                    {
                        Jcoeffs.push_back(Triplet<double>(term, weave.vidx(f, i) + k, dE[k]));
                    }
                    for (int field = 0; field < m; field++)
                    {
                        dE = -permut(i, field) * Tgf.transpose() * innervec;
                        for (int k = 0; k < 2; k++)
                        {
                            Jcoeffs.push_back(Triplet<double>(term, weave.vidx(g, field) + k, dE[k]));
                        }
                    }
                    Jcoeffs.push_back(Triplet<double>(term, weave.alphaidx(f, i), innervec.dot(vif) * vif.dot(cdiff)));
                    dE = (Jf * vif).dot(cdiff) * innervec;
                    for (int k = 0; k < 2; k++)
                    {
                        Jcoeffs.push_back(Triplet<double>(term, weave.betaidx(f, i) + k, dE[k]));
                    }
                    term++;
                }
            }
        }
    }
    J.setFromTriplets(Jcoeffs.begin(), Jcoeffs.end());
}

void testFiniteDifferences(Weave &weave, SolverParams params)
{
    Weave test = weave;
    test.vectorFields.setRandom();

    double orig = energy(test, params);
    VectorXd dE;
    trueGradient(test, params, dE);

    std::ofstream ofs("dump.txt");

    for (int i = 0; i < test.vectorFields.size(); i++)
    {
        Weave perturb = test;
        perturb.vectorFields[i] += 1e-6;
        double newe = energy(perturb, params);
        double findiff = (newe - orig) / 1e-6;
        ofs << findiff << " vs " << dE[i] << std::endl;
    }
}

void GNtestFiniteDifferences(Weave &weave, SolverParams params)
{
    Weave test = weave;
    test.vectorFields.setRandom();

    Eigen::VectorXd orig;
    GNEnergy(test, params, orig);
    Eigen::SparseMatrix<double> J;
    GNGradient(test, params, J);

    std::ofstream ofs("dump.txt");

    for (int i = 0; i < test.vectorFields.size(); i++)
    {
        Weave perturb = test;
        perturb.vectorFields[i] += 1e-6;
        Eigen::VectorXd newE;
        GNEnergy(perturb, params, newE);
        Eigen::VectorXd findiff = (newE - orig) / 1e-6;
        Eigen::VectorXd exact = J.col(i);
        ofs << (findiff - exact).norm() << " / " << exact.norm() << std::endl;
    }
}