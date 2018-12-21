#include "GaussNewton.h"
#include "Weave.h"
#include <iostream>
#include <fstream>

using namespace Eigen;

/*

void faceEnergies(const Weave &weave, SolverParams params, Eigen::MatrixXd &E)
{
    int nhandles = weave.nHandles();
    int nedges = weave.nEdges();
    int m = weave.nFields();
    int nfaces = weave.nFaces();
    E.resize(nfaces, m + 4);
    E.setZero();

    Eigen::VectorXd r;
    GNEnergy(weave, params, r);
    int term = 2 * nhandles;

    for (int e = 0; e < nedges; e++)
    {
        if(weave.E(e, 0) == -1 || weave.E(e, 1) == -1)
            continue;
        for (int i = 0; i < m; i++)
        {
            for (int side = 0; side < 2; side++)
            {
                int f = (side == 0 ? weave.E(e, 0) : weave.E(e, 1));

                Eigen::Vector2d termr = r.segment<2>(term);
                E(f, i) += 0.5 * termr.transpose() * weave.Bs[f].transpose() * weave.Bs[f] * termr;

                term += 2;
            }
        }
    }

    for (int i = 0; i < nfaces; i++)
    {
        E(i, m)     = weave.v(i, 0).norm();
        E(i, m + 1) = weave.beta(i, 0).transpose() * weave.Js.block<2,2>(2*i,0) * weave.v(i, 0);
        E(i, m + 2) = weave.beta(i, 0).normalized().transpose() * weave.Js.block<2,2>(2*i,0) * weave.v(i, 0).normalized();
        E(i, m + 3) = weave.beta(i, 0)[0];
        E(i, m + 1) = weave.alpha(i, 0) ;
  //      std::cout << E(i, m + 1) << '\n';
    } 
}

void GNmetric(const Weave &weave, Eigen::SparseMatrix<double> &M)
{
    int nhandles = weave.nHandles();
    int nedges = weave.nEdges();
    int nfaces = weave.nFaces();
    int m = weave.nFields();
    int intedges = weave.numInteriorEdges();
    int numterms = 2*nhandles + 5 * intedges * m;
    M.resize(numterms, numterms);

    std::vector<Eigen::Triplet<double> > Mcoeffs;

    for (int i = 0; i < nhandles; i++)
    {
        int face = weave.handles[i].face;
        Eigen::Matrix2d BTB = weave.Bs[face].transpose() * weave.Bs[face];
        for (int j = 0; j < 2; j++)
        {
            for (int k = 0; k < 2; k++)
            {
                Mcoeffs.push_back(Triplet<double>(2 * i + j, 2 * i + k, BTB(j, k)));
            }
        }
    }

    int term = 2 * nhandles;

    for (int e = 0; e < nedges; e++)
    {
        if(weave.E(e,0) == -1 || weave.E(e,1) == -1)
            continue;
        for (int i = 0; i < m; i++)
        {
            for (int side = 0; side < 2; side++)
            {
                int f = (side == 0 ? weave.E(e, 0) : weave.E(e, 1));
                
                Eigen::Matrix2d BTB = weave.Bs[f].transpose() * weave.Bs[f];
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        Mcoeffs.push_back(Triplet<double>(term + j, term + k, BTB(j, k)));
                    }
                }
                term += 2;
            }
        }
    }

    for (int e = 0; e < nedges; e++)
    {
        if(weave.E(e,0) == -1 || weave.E(e,1) == -1)
            continue;
        for (int i = 0; i < m; i++)
        {
                Mcoeffs.push_back(Triplet<double>(term, term, 1));
                term += 1;
        }
    }

    // for (int f = 0; f < nfaces; f++)
    // {
    //     for (int field = 0; field < m; field++)
    //     {
    //         Eigen::Matrix2d BTB = weave.Bs[f].transpose() * weave.Bs[f];
    //         for (int j = 0; j < 2; j++)
    //         {
    //             for (int k = 0; k < 2; k++)
    //             {
    //                 Mcoeffs.push_back(Triplet<double>(term + j, term + k, 1));
    //             }
    //         }
    //         term += 2;
    //     }
    // }

    M.setFromTriplets(Mcoeffs.begin(), Mcoeffs.end());
}

// Here we define the energy of field optimization operation. 
void GNEnergy(const Weave &weave, SolverParams params, Eigen::VectorXd &E)
{
    int nhandles = weave.nHandles();
    int nedges = weave.nEdges();
    int nfaces = weave.nFaces();
    int m = weave.nFields();
    int intedges = weave.numInteriorEdges();
    int nterms = 2 * nhandles + 5 * intedges*m;
    E.resize(nterms);
    E.setZero();

    // Energy for aligning with the user prescribed handles
    for(int i=0; i<nhandles; i++)
    { 
        int field = weave.handles[i].field;
        int face = weave.handles[i].face;
        Vector2d dir = weave.handles[i].dir;
        for (int coeff = 0; coeff < 2; coeff++)
            if (!weave.fixFields)
                E[2 * i + coeff] = (weave.v(face, field) - dir)[coeff];
    }

    int term = 2*nhandles;

    // compatibility constraint
    for (int e = 0; e < nedges; e++)
    {
        if(weave.E(e,0) == -1 || weave.E(e,1) == -1)
            continue;
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
                for (int k = 0; k < 2; k++)
                {
                    E[term] = sqrt(params.edgeWeights(e) * params.lambdacompat) * (Dif*cdiff - (Tgf*vpermut - vif))[k];
                    term++;
                }
            }
        }
    } 

    // Curl free term?
    for (int e = 0; e < nedges; e++)
    {
        if(weave.E(e,0) == -1 || weave.E(e,1) == -1)
            continue;
        for (int i = 0; i < m; i++)
        {
                int f = weave.E(e, 0);
                int g = weave.E(e, 1);
                Eigen::Vector3d edge = weave.V.row(weave.edgeVerts(e, 0)) - weave.V.row(weave.edgeVerts(e, 1));
                Eigen::Matrix2d Jf = weave.Js.block<2, 2>(2 * f, 0);
                Eigen::Matrix2d Jg = weave.Js.block<2, 2>(2 * g, 0);
                Eigen::Vector2d vif = weave.v(f, i);
                Eigen::Vector2d vpermut(0, 0);
                Eigen::MatrixXi permut = weave.Ps[e];
                for (int field = 0; field < m; field++)
                {
                    vpermut += permut(i, field) * weave.v(g, field);  // Ask about this, v confuse
                }
                E[term] = (weave.Bs[f] * (vif)).dot(edge) - (weave.Bs[g] * ( vpermut)).dot(edge);
                E[term] *= params.curlreg;
                term++;
        }
    }   

    // for (int f = 0; f < nfaces; f++)
    // {
    //     for (int i = 0; i < m; i++)
    //     {
    //         Eigen::Vector2d b = weave.beta(f, m);
    //         E[term] = b[1] * params.curlreg; // Expose handle, this is the curl free parameter.
    //         E[term + 1] = 0;//(b[1] * b[1] - 1.);
    //         term += 2;
    //     }
    // }
}

void GNGradient(const Weave &weave, SolverParams params, Eigen::SparseMatrix<double> &J)
{
    int nhandles = weave.nHandles();
    int nedges = weave.nEdges();
    int nfaces = weave.nFaces();
    int m = weave.nFields();
    int intedges = weave.numInteriorEdges();
    int nterms = 2 * nhandles + 5 * intedges*m;
    J.resize(nterms, weave.vectorFields.size());

    std::vector<Eigen::Triplet<double> > Jcoeffs;

    for(int i=0; i<nhandles; i++)
    { 
        int field = weave.handles[i].field;
        int face = weave.handles[i].face;
        for (int coeff = 0; coeff < 2; coeff++)
        {
            if (!weave.fixFields)
                Jcoeffs.push_back(Triplet<double>(2 * i + coeff, weave.vidx(face, field) + coeff, 1.0));
        }
    }

    int term = 2 * nhandles;

    // compatibility constraint
    for (int e = 0; e < nedges; e++)
    {
        if(weave.E(e,0) == -1 || weave.E(e,1) == -1)
            continue;
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

                for (int coeff = 0; coeff < 2; coeff++)
                {
                    Eigen::Vector2d innervec(0, 0);
                    innervec[coeff] = sqrt(params.edgeWeights(e) * params.lambdacompat);
                    Vector2d dE = Jf.transpose()*cdiff * weave.beta(f, i).dot(innervec);
                    dE += cdiff.dot(vif) * weave.alpha(f, i) * innervec;
                    dE += cdiff * weave.alpha(f, i) * vif.dot(innervec);
                    dE += innervec;
                    for (int k = 0; k < 2; k++)
                    {
                        if (!weave.fixFields)
                            Jcoeffs.push_back(Triplet<double>(term, weave.vidx(f, i) + k, dE[k]));
                    }
                    for (int field = 0; field < m; field++)
                    {
                        dE = -permut(i, field) * Tgf.transpose() * innervec;
                        for (int k = 0; k < 2; k++)
                        {
                            if (!weave.fixFields)
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

    for (int e = 0; e < nedges; e++)
    {
        if(weave.E(e,0) == -1 || weave.E(e,1) == -1)
            continue;
        for (int i = 0; i < m; i++)
        {
                int f = weave.E(e, 0);
                int g = weave.E(e, 1);
                Eigen::Vector3d edge = weave.V.row(weave.edgeVerts(e, 0)) - weave.V.row(weave.edgeVerts(e, 1));
                Eigen::Matrix2d Jf = weave.Js.block<2, 2>(2 * f, 0);
                Eigen::Matrix2d Jg = weave.Js.block<2, 2>(2 * g, 0);
                Eigen::Vector2d vif = weave.v(f, i);
                Eigen::Vector2d b0(1, 0);
                Eigen::Vector2d b1(0, 1);
                Eigen::MatrixXi permut = weave.Ps[e];
                for (int k = 0; k < 2; k++)
                {
                    Jcoeffs.push_back(Triplet<double>(term, weave.vidx(f, i) + k, (edge.transpose() * weave.Bs[f] )[k] * params.curlreg));
                }

                for (int field = 0; field < m; field++)
                {
                    Eigen::Vector2d dE = -permut(i, field) * (edge.transpose() * weave.Bs[g] )* params.curlreg;
                    for (int k = 0; k < 2; k++)
                    {
                        Jcoeffs.push_back(Triplet<double>(term, weave.vidx(g, field) + k, dE[k]));
                    }
                }
                
                term++;
        }
    } 
    // for (int f = 0; f < nfaces; f++)
    // {
    //     for (int i = 0; i < m; i++)
    //     {
    //   //      Eigen::Vector2d b = weave.beta(f, m);
    //         Jcoeffs.push_back(Triplet<double>(term, weave.betaidx(f, i),  params.curlreg));
    //         Jcoeffs.push_back(Triplet<double>(term + 1, weave.betaidx(f, i),  0)); //b[1] * b[1] - 1.));
    //         term += 2;
    //     }
    // }

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
    int nvars = weave.vectorFields.size();
    Eigen::VectorXd r;
    GNEnergy(weave, params, r);
    Eigen::SparseMatrix<double> M;
    GNmetric(weave, M);
  
    Eigen::VectorXd vectorFields_fix;   
    if (weave.fixFields)
    {
        vectorFields_fix = weave.vectorFields;
    }

    std::cout << "original energy: " << 0.5 * r.transpose() * M * r << std::endl;
    std::cout << "Building matrix" << std::endl;
    Eigen::SparseMatrix<double> J;
    GNGradient(weave, params, J);
    Eigen::SparseMatrix<double> optMat(nvars, nvars);
    std::vector<Eigen::Triplet<double> > coeffs;
    int nfaces = weave.nFaces();
    int m = weave.nFields();
   // for (int i = 2 * nfaces*m; i < 5 * nfaces*m; i++)
    for (int i = 0; i < 5 * nfaces*m; i++)
        coeffs.push_back(Eigen::Triplet<double>(i, i, params.lambdareg));
    optMat.setFromTriplets(coeffs.begin(), coeffs.end());
    optMat += J.transpose() * M * J;
    std::cout << "Done, " << optMat.nonZeros() << " nonzeros" << std::endl;
    optMat.makeCompressed();
    Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double> > solver;
    std::cout << "Analyzing" << std::endl;
    solver.analyzePattern(optMat);
    Eigen::VectorXd rhs = J.transpose() * M * r;
    std::cout << "Solving" << std::endl;
    solver.factorize(optMat);
    Eigen::VectorXd update = solver.solve(rhs);
    lineSearch(weave, params, update);

    if (weave.fixFields)
    { 
        int nfaces = weave.nFaces();
        int m = weave.nFields();
        for (int i = 0; i < nfaces; i++)
        {
            for (int j = 0; j < m; j++)
            {
                int vid = weave.vidx(i, j);
                weave.vectorFields.segment<2>(vid) = vectorFields_fix.segment<2>(vid);
            }
        }
    }
    
    GNEnergy(weave, params, r);
    std::cout << "Done, new energy: " << 0.5 * r.transpose()*M*r << std::endl;
}

double lineSearch(Weave &weave, SolverParams params, const Eigen::VectorXd &update)
{
    double t = 1.0;
    double c1 = 0.1;
    double c2 = 0.9;
    double alpha = 0;
    double infinity = 1e6;
    double beta = infinity;

    int nvars = weave.vectorFields.size();
    Eigen::VectorXd r;
    GNEnergy(weave, params, r);
    Eigen::SparseMatrix<double> M;
    GNmetric(weave, M);
    Eigen::SparseMatrix<double> J;
    GNGradient(weave, params, J);

    VectorXd dE;
    VectorXd newdE;
    VectorXd startVF = weave.vectorFields;
    
    double orig = 0.5 * r.transpose() * M * r;
    dE = J.transpose() * M * r;
    double deriv = -dE.dot(update);
    assert(deriv < 0);
    
    std::cout << "Starting line search, original energy " << orig << ", descent magnitude " << deriv << std::endl;

    while (true)
    {
        weave.vectorFields = startVF - t * update;
        GNEnergy(weave, params, r);
        double newenergy = 0.5 * r.transpose() * M * r;
        GNGradient(weave, params, J);
        newdE = J.transpose() * M * r;

        std::cout << "Trying t = " << t << ", energy now " << newenergy << std::endl;
        
        if (std::isnan(newenergy) || newenergy > orig + t*deriv*c1)
        {
            beta = t;
            t = 0.5*(alpha + beta);
        }
        else if (-newdE.dot(update) < c2*deriv)
        {
            alpha = t;
            if (beta == infinity)
            {
                t = 2 * alpha;
            }
            else
            {
                t = 0.5*(alpha + beta);
            }

            if (beta - alpha < 1e-8)
            {
                return t;
            }
        }
        else
        {
            return t;
        }
    }
}


*/