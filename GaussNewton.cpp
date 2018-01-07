#include "GaussNewton.h"
#include "Weave.h"
#include <iostream>
#include <fstream>

using namespace Eigen;

void GNEnergy(const Weave &weave, SolverParams params, Eigen::VectorXd &E)
{
    int nhandles = weave.nHandles();
    int nedges = weave.nEdges();
    int m = weave.nFields();

    int nterms = 3 * nhandles + 6 * nedges*m;
    E.resize(nterms);
    E.setZero();

    for(int i=0; i<nhandles; i++)
    { 
        int field = weave.handles[i].field;
        int face = weave.handles[i].face;
        Vector2d dir = weave.handles[i].dir;
        for (int k = 0; k < 3; k++)
            E[3 * i + k] = (weave.Bs[face] * (weave.v(face, field) - dir))[k];
    }

    int term = 3*nhandles;

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
                Eigen::MatrixXi permut = weave.Ps[e];
                if (side == 1)
                    permut.transposeInPlace();
                for (int field = 0; field < m; field++)
                {
                    vpermut += permut(i, field) * weave.v(g, field);
                }
                Eigen::Matrix2d Tgf = weave.Ts.block<2, 2>(2 * e, 2 - 2 * side);
                for (int k = 0; k < 3; k++)
                {
                    E[term] = sqrt(params.lambdacompat) * (weave.Bs[f] * (Dif*cdiff - (Tgf*vpermut - vif)))[k];
                    term++;
                }
            }
        }
    }    
}

void GNGradient(const Weave &weave, SolverParams params, Eigen::SparseMatrix<double> &J)
{
    int nhandles = weave.nHandles();
    int nedges = weave.nEdges();
    int m = weave.nFields();

    int nterms = 3 * nhandles + 6 * nedges*m;
    J.resize(nterms, weave.vectorFields.size());

    std::vector<Eigen::Triplet<double> > Jcoeffs;

    for(int i=0; i<nhandles; i++)
    { 
        int field = weave.handles[i].field;
        int face = weave.handles[i].face;
        for (int coeff = 0; coeff < 3; coeff++)
        {
            for (int k = 0; k < 2; k++)
            {
                Jcoeffs.push_back(Triplet<double>(3 * i + coeff, weave.vidx(face, field) + k, weave.Bs[face](coeff,k)));
            }
        }
    }

    int term = 3 * nhandles;

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
                Eigen::MatrixXi permut = weave.Ps[e];
                if (side == 1)
                    permut.transposeInPlace();
                for (int field = 0; field < m; field++)
                {
                    vpermut += permut(i, field) * weave.v(g, field);
                }
                Eigen::Matrix2d Tgf = weave.Ts.block<2, 2>(2 * e, 2 - 2 * side);

                for (int coeff = 0; coeff < 3; coeff++)
                {
                    Eigen::Vector2d innervec(0, 0);
                    innervec = weave.Bs[f].row(coeff).transpose() * sqrt(params.lambdacompat);
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

void oneStep(Weave &weave, SolverParams params)
{
    Eigen::VectorXd r;
    GNEnergy(weave, params, r);


    std::cout << "original energy: " << r.squaredNorm() << std::endl;;

    Eigen::SparseMatrix<double> J;
    GNGradient(weave, params, J);
    Eigen::SparseMatrix<double> M = J.transpose() * J;
    Eigen::SparseMatrix<double> I(weave.vectorFields.size(), weave.vectorFields.size());
    I.setIdentity();
    M += params.lambdareg * I;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    solver.compute(M);
    Eigen::VectorXd rhs = J.transpose() * r;
    Eigen::VectorXd update = solver.solve(rhs);
    weave.vectorFields -= update;

    GNEnergy(weave, params, r);
    std::cout << "new energy: " << r.squaredNorm() << std::endl;;
}