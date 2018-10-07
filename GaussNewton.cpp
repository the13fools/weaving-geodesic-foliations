#include "GaussNewton.h"
#include "Weave.h"
#include <iostream>
#include <fstream>
#include <Eigen/CholmodSupport>
#include "Surface.h"
#include <igl/cotmatrix.h>

#include <Eigen/Eigenvalues>
#include <vector>       // std::vector

#include <cmath>

using namespace Eigen;

void faceEnergies(const Weave &weave, SolverParams params, Eigen::MatrixXd &E)
{
    int nhandles = weave.nHandles();
    int nedges = weave.fs->nEdges();
    int m = weave.fs->nFields();
    int nfaces = weave.fs->nFaces();
    E.resize(nfaces, m);
    E.setZero();

    Eigen::VectorXd r;
    GNEnergy(weave, params, r);
    int term = 2 * nhandles;

    for (int e = 0; e < nedges; e++)
    {
        if(weave.fs->data().E(e, 0) == -1 || weave.fs->data().E(e, 1) == -1)
            continue;
        for (int i = 0; i < m; i++)
        {
            for (int side = 0; side < 2; side++)
            {
                int f = (side == 0 ? weave.fs->data().E(e, 0) : weave.fs->data().E(e, 1));

                Eigen::Vector2d termr = r.segment<2>(term);
                E(f, i) += 0.5 * termr.transpose() * weave.fs->data().Bs[f].transpose() * weave.fs->data().Bs[f] * termr;

                term += 2;
            }
        }
    }
}

void GNmetric(const Weave &weave, Eigen::SparseMatrix<double> &M)
{
    int nhandles = weave.nHandles();
    int nedges = weave.fs->nEdges();
    int m = weave.fs->nFields();
    int intedges = weave.fs->numInteriorEdges();
    int numterms = 2*nhandles + 6 * intedges * m;
    M.resize(numterms, numterms);
    M.setZero();

    std::vector<Eigen::Triplet<double> > Mcoeffs;

    for (int i = 0; i < nhandles; i++)
    {
        int face = weave.handles[i].face;
        Eigen::Matrix2d BTB = weave.fs->data().Bs[face].transpose() * weave.fs->data().Bs[face];
        for (int j = 0; j < 2; j++)
        {
            for (int k = 0; k < 2; k++)
            {
                Mcoeffs.push_back(Triplet<double>(2 * i + j, 2 * i + k, BTB(j, k)));
            }
        }
    }

    int term = 2 * nhandles;

    // compatibliity terms
    for (int e = 0; e < nedges; e++)
    {
        if(weave.fs->data().E(e,0) == -1 || weave.fs->data().E(e,1) == -1)
            continue;
        for (int i = 0; i < m; i++)
        {
            for (int side = 0; side < 2; side++)
            {
                int f = (side == 0 ? weave.fs->data().E(e, 0) : weave.fs->data().E(e, 1));
                double area = weave.fs->faceArea(f);
                Eigen::Matrix2d BTB = weave.fs->data().Bs[f].transpose() * weave.fs->data().Bs[f];
                for (int j = 0; j < 2; j++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        Mcoeffs.push_back(Triplet<double>(term + j, term + k, area*BTB(j, k)));
                    }
                }
                term += 2;
            }
        }
    }


    // build edge metric matrix and inverse (cotan weights)
    Eigen::MatrixXd C;
    igl::cotmatrix_entries(weave.fs->data().V, weave.fs->data().F, C);

    Eigen::VectorXd edgeMetric(nedges);
    edgeMetric.setZero();
    int nfaces = weave.fs->nFaces();
    for(int i=0; i<nfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int eidx = weave.fs->data().faceEdges(i, j);
            edgeMetric[eidx] += C(i,j);
        }                
    } 

    // Curl free terms
    for (int e = 0; e < nedges; e++)
    {
        if(weave.fs->data().E(e,0) == -1 || weave.fs->data().E(e,1) == -1)
            continue;
        else 
        {
            for (int i = 0; i < m; i++)
            {
                 //   Mcoeffs.push_back(Triplet<double>(term, term, edgeMetric[e]));
                    Mcoeffs.push_back(Triplet<double>(term, term, 1));
                    term += 1;
            }
        }
    }

        // Unit Smoothness terms
    for (int e = 0; e < nedges; e++)
    {
        if(weave.fs->data().E(e,0) == -1 || weave.fs->data().E(e,1) == -1)
            continue;
        else 
        {
            for (int i = 0; i < m; i++)
            {
                 //   Mcoeffs.push_back(Triplet<double>(term, term, edgeMetric[e]));
                    Mcoeffs.push_back(Triplet<double>(term, term, 1));
                    term += 1;
            }
        }
    }

    std::cout << intedges << " " << term << std::endl;

    // int s|v| = Area constraint
    // Mcoeffs.push_back(Triplet<double>(term, term, 1));
    // term++;

    M.setFromTriplets(Mcoeffs.begin(), Mcoeffs.end());
}

void GNEnergy(const Weave &weave, SolverParams params, Eigen::VectorXd &E)
{
    int nhandles = weave.nHandles();
    int nedges = weave.fs->nEdges();
    int m = weave.fs->nFields();
    int intedges = weave.fs->numInteriorEdges();
    int nterms = 2 * nhandles + 6 * intedges*m;
//    std::cout << "nhandles: " << nhandles << std::endl;
    E.resize(nterms);
    E.setZero();

    for(int i=0; i<nhandles; i++)
    { 
        int field = weave.handles[i].field;
        int face = weave.handles[i].face;
        Vector2d dir = weave.handles[i].dir;
        for (int coeff = 0; coeff < 2; coeff++)
            E[2 * i + coeff] = (weave.fs->v(face, field) - dir)[coeff] * params.handleWeight;
    }

    int term = 2*nhandles;

    // compatibility constraint
    for (int e = 0; e < nedges; e++)
    {
        if(weave.fs->data().E(e,0) == -1 || weave.fs->data().E(e,1) == -1)
            continue;
        for (int i = 0; i < m; i++)
        {
            for (int side = 0; side < 2; side++)
            {
                int f = (side == 0 ? weave.fs->data().E(e, 0) : weave.fs->data().E(e, 1));
                int g = (side == 0 ? weave.fs->data().E(e, 1) : weave.fs->data().E(e, 0));

                Eigen::Vector2d vif = weave.fs->v(f, i);

                Eigen::Vector2d vpermut(0, 0);
                Eigen::MatrixXi permut = weave.fs->Ps(e);
                if (side == 1)
                    permut.transposeInPlace();
                for (int field = 0; field < m; field++)
                {
                    vpermut += permut(i, field) * weave.fs->v(g, field);
                }
                Eigen::Matrix2d Tgf = weave.fs->data().Ts.block<2, 2>(2 * e, 2 - 2 * side);
                for (int k = 0; k < 2; k++)
                {
                    E[term] = sqrt(params.edgeWeights(e) * params.lambdacompat) * (Tgf*vpermut - vif)[k];
                    term++;
                }
            }
        }
    }

    // Curl free term
    for (int e = 0; e < nedges; e++)
    {
        if(weave.fs->data().E(e,0) == -1 || weave.fs->data().E(e,1) == -1)
            continue;
        for (int i = 0; i < m; i++)
        {
            int f = weave.fs->data().E(e, 0);
            int g = weave.fs->data().E(e, 1);
            Eigen::Vector3d edge = weave.fs->data().V.row(weave.fs->data().edgeVerts(e, 0)) - weave.fs->data().V.row(weave.fs->data().edgeVerts(e, 1));
            Eigen::Vector2d vif = weave.fs->v(f, i);
            Eigen::Vector2d vpermut(0, 0);
            Eigen::MatrixXi permut = weave.fs->Ps(e);
            for (int field = 0; field < m; field++)
            {
                vpermut += permut(i, field) * weave.fs->v(g, field);  
            }
            double n_vif = (weave.fs->data().Bs[f] * vif).norm();
            double n_vpermut = (weave.fs->data().Bs[g] * vpermut).norm();
            vif *= weave.fs->sval(f, i);
            for (int field = 0; field < m; field++)
            {
                if (permut(i, field) != 0)
                    vpermut *= weave.fs->sval(g, field);  
            }

        //    E[term] = (weave.fs->data().Bs[f] * (vif)).dot(edge) - (weave.fs->data().Bs[g] * (vpermut)).dot(edge);
            E[term] = (weave.fs->data().Bs[f] * (vif)).dot(edge) - (weave.fs->data().Bs[g] * (vpermut)).dot(edge);
            E[term] *= params.curlLambda * params.edgeWeights(e);
            term++;
        }
    }   


    // SV unit + smooth
    for (int e = 0; e < nedges; e++)
    {
        if(weave.fs->data().E(e,0) == -1 || weave.fs->data().E(e,1) == -1)
            continue;
        for (int i = 0; i < m; i++)
        {
            int f = weave.fs->data().E(e, 0);
            int g = weave.fs->data().E(e, 1);
            Eigen::Vector3d edge = weave.fs->data().V.row(weave.fs->data().edgeVerts(e, 0)) - weave.fs->data().V.row(weave.fs->data().edgeVerts(e, 1));
            Eigen::Vector2d vif = weave.fs->v(f, i);
            Eigen::Vector2d vpermut(0, 0);
            Eigen::MatrixXi permut = weave.fs->Ps(e);
            for (int field = 0; field < m; field++)
            {
                vpermut += permut(i, field) * weave.fs->v(g, field);  
            }
            double n_vif = (weave.fs->data().Bs[f] * vif).norm();
            double n_vpermut = (weave.fs->data().Bs[g] * vpermut).norm();

            double vpermut_sval = 0;
            for (int field = 0; field < m; field++)
            {
                if (permut(i, field) != 0)
                    vpermut_sval = weave.fs->sval(g, field);  
            }

        //    E[term] = (weave.fs->data().Bs[f] * (vif)).dot(edge) - (weave.fs->data().Bs[g] * (vpermut)).dot(edge);
            E[term] = weave.fs->sval(f, i) * n_vif - vpermut_sval * n_vpermut;
            E[term] *= params.unitLambda * params.edgeWeights(e);
            term++;
        }
    }


    // s|v| = area constraint.
    // int nfaces = weave.fs->nFaces();
    // double tot_product = 0.;
    // for (int f = 0; f < nfaces; f++)
    // {
    //     for (int field = 0; field < m; field++)
    //     {
    //         Eigen::Vector2d update = weave.fs->v(f, field);
    //         tot_product += weave.fs->sval(f, field) * update.norm();
    //     }
    // }

    // E[term] = (nfaces - tot_product) * params.unitAreaLambda;
    // term++;

}

void GNGradient(const Weave &weave, SolverParams params, Eigen::SparseMatrix<double> &J)
{
    int nhandles = weave.nHandles();
    int nedges = weave.fs->nEdges();
    int m = weave.fs->nFields();
    int intedges = weave.fs->numInteriorEdges();
    int nterms = 2 * nhandles + 6 * intedges*m;
    J.resize(nterms, weave.fs->vectorFields.size());
    J.setZero();

    std::vector<Eigen::Triplet<double> > Jcoeffs;

    for(int i=0; i<nhandles; i++)
    { 
        int field = weave.handles[i].field;
        int face = weave.handles[i].face;
        for (int coeff = 0; coeff < 2; coeff++)
        {
            Jcoeffs.push_back(Triplet<double>(2 * i + coeff, weave.fs->vidx(face, field) + coeff, params.handleWeight));
        }
    }

    int term = 2 * nhandles;

    // compatibility constraint
    for (int e = 0; e < nedges; e++)
    {
        if(weave.fs->data().E(e,0) == -1 || weave.fs->data().E(e,1) == -1)
            continue;
        for (int i = 0; i < m; i++)
        {
            for (int side = 0; side < 2; side++)
            {
                int f = (side == 0 ? weave.fs->data().E(e, 0) : weave.fs->data().E(e, 1));
                int g = (side == 0 ? weave.fs->data().E(e, 1) : weave.fs->data().E(e, 0));
                Eigen::Matrix2d Jf = weave.fs->data().Js.block<2, 2>(2 * f, 0);
                Eigen::Vector2d vif = weave.fs->v(f, i);
                Eigen::Matrix2d Dif = weave.fs->beta(f, i) * (Jf*vif).transpose() + weave.fs->alpha(f,i) * vif * vif.transpose();
                Eigen::Vector2d cdiff = weave.fs->data().cDiffs.row(2 * e + side);
                Eigen::Vector2d vpermut(0, 0);
                Eigen::MatrixXi permut = weave.fs->Ps(e);
                if (side == 1)
                    permut.transposeInPlace();
                for (int field = 0; field < m; field++)
                {
                    vpermut += permut(i, field) * weave.fs->v(g, field);
                }
                Eigen::Matrix2d Tgf = weave.fs->data().Ts.block<2, 2>(2 * e, 2 - 2 * side);

                for (int coeff = 0; coeff < 2; coeff++)
                {
                    Eigen::Vector2d innervec(0, 0);
                    innervec[coeff] = sqrt(params.edgeWeights(e) * params.lambdacompat);
                 //   Vector2d dE = Jf.transpose()*cdiff * weave.fs->beta(f, i).dot(innervec);
                    Vector2d dE = -innervec;
                    for (int k = 0; k < 2; k++)
                    {
                        Jcoeffs.push_back(Triplet<double>(term, weave.fs->vidx(f, i) + k, dE[k]));
                    }
                    for (int field = 0; field < m; field++)
                    {
                        dE = permut(i, field) * Tgf.transpose() * innervec;
                        for (int k = 0; k < 2; k++)
                        {
                            Jcoeffs.push_back(Triplet<double>(term, weave.fs->vidx(g, field) + k, dE[k]));
                        }
                    }

                    term++;
                }
            }
        }
    }

    // (weave.fs->data().Bs[f] * (vif)).dot(edge) - (weave.fs->data().Bs[g] * (vpermut)).dot(edge);

    // (s * (x'*x) - s2 * w'*w)^2

    // Curl correction
    for (int e = 0; e < nedges; e++)
    {
        if(weave.fs->data().E(e,0) == -1 || weave.fs->data().E(e,1) == -1)
            continue;
        for (int i = 0; i < m; i++)
        {
                int f = weave.fs->data().E(e, 0);
                int g = weave.fs->data().E(e, 1);
                Eigen::Vector3d edge = weave.fs->data().V.row(weave.fs->data().edgeVerts(e, 0)) - 
                                           weave.fs->data().V.row(weave.fs->data().edgeVerts(e, 1));
                Eigen::Vector2d vif = weave.fs->v(f, i);
                Eigen::Vector2d b0(1, 0);
                Eigen::Vector2d b1(0, 1);
                double n_v = (weave.fs->data().Bs[f] * vif).norm();
                Eigen::Matrix<double, 3, 2> B_f = weave.fs->data().Bs[f];
                Eigen::Vector2d e_f = edge.transpose() * B_f;
         //       Eigen::Vector2d dE = (e_f );
           //     std::cout << n_v << " ";
                Eigen::Vector2d dE = ( e_f ) * weave.fs->sval(f, i);
                // dE -= vif.transpose() * B_f.transpose() * edge * vif.transpose() * B_f.transpose() * B_f / (n_v * n_v * n_v);
                // dE *= weave.fs->sval(f, i);
                for (int k = 0; k < 2; k++)
                {
                    Jcoeffs.push_back(Triplet<double>(term, weave.fs->vidx(f, i) + k, dE[k] * params.curlLambda * params.edgeWeights(e)));
                }
                 
                Eigen::Vector2d vperm(0, 0);
                Eigen::MatrixXi permut = weave.fs->Ps(e);
                for (int field = 0; field < m; field++)
                {
                    vperm += permut(i, field) * weave.fs->v(g, field);  
                }
                Eigen::Matrix<double, 3, 2> B_g = weave.fs->data().Bs[g];
                double n_vperm = (B_g * (vperm)).norm();
                Eigen::Vector2d e_g = edge.transpose() * B_g; 
                for (int field = 0; field < m; field++)
                {
      //              Eigen::Vector2d dE = permut(i, field) * (e_g);
                    Eigen::Vector2d dE = (e_g); 
                //    dE -= vperm.transpose() * B_g.transpose() * edge * vperm.transpose() * B_g.transpose() * B_g  / (n_vperm * n_vperm * n_vperm);
                    dE *= permut(i, field) * weave.fs->sval(g, field);  
                    for (int k = 0; k < 2; k++)
                    {
                        Jcoeffs.push_back(Triplet<double>(term, weave.fs->vidx(g, field) + k, -dE[k] * params.curlLambda * params.edgeWeights(e)));
                    }
                }
                
                term++;
        }
    } 

    // sv unit + smooth d/dv.  
    for (int e = 0; e < nedges; e++)
    {
        if(weave.fs->data().E(e,0) == -1 || weave.fs->data().E(e,1) == -1)
            continue;
        for (int i = 0; i < m; i++)
        {
                int f = weave.fs->data().E(e, 0);
                int g = weave.fs->data().E(e, 1);
                Eigen::Vector3d edge = weave.fs->data().V.row(weave.fs->data().edgeVerts(e, 0)) - 
                                           weave.fs->data().V.row(weave.fs->data().edgeVerts(e, 1));
                Eigen::Vector2d vif = weave.fs->v(f, i);
                Eigen::Matrix<double, 3, 2> B_f = weave.fs->data().Bs[f];
                double n_v = (B_f * vif).norm();
                Eigen::Vector2d dE = (B_f * vif).transpose() * B_f * weave.fs->sval(f, i) / n_v;

                for (int k = 0; k < 2; k++)
                {
                    Jcoeffs.push_back(Triplet<double>(term, weave.fs->vidx(f, i) + k, dE[k] * params.unitLambda * params.edgeWeights(e)));
              //      Jcoeffs.push_back(Triplet<double>(nterms - 1, weave.fs->vidx(f, i) + k, dE[k] * params.unitAreaLambda * params.edgeWeights(e)));
                }
                 
                Eigen::Vector2d vperm(0, 0);
                Eigen::MatrixXi permut = weave.fs->Ps(e);
                for (int field = 0; field < m; field++)
                {
                    vperm += permut(i, field) * weave.fs->v(g, field);  
                }
                Eigen::Matrix<double, 3, 2> B_g = weave.fs->data().Bs[g];
                double n_vperm = (B_g * (vperm)).norm();
                Eigen::Vector2d dE_g = (B_g * vperm).transpose() * B_g / n_vperm;
                for (int field = 0; field < m; field++)
                {
                    Eigen::Vector2d dE = dE_g; 
                    dE *= permut(i, field) * weave.fs->sval(g, field);  
                    for (int k = 0; k < 2; k++)
                    {
                        Jcoeffs.push_back(Triplet<double>(term, weave.fs->vidx(g, field) + k, -dE[k] * params.unitLambda * params.edgeWeights(e)));
                  //      Jcoeffs.push_back(Triplet<double>(nterms - 1, weave.fs->vidx(g, field) + k, -dE[k] * params.unitAreaLambda * params.edgeWeights(e)));
                    }
                }
                
                term++;
        }
    } 

    for (auto t : Jcoeffs)
    {
        if (std::isnan(t.value()) || std::isinf(t.value()))
        {
            std::cout << "oh no!";
            exit(-1);
        }
    }
    J.setFromTriplets(Jcoeffs.begin(), Jcoeffs.end());
}

void GNtestFiniteDifferences(Weave &weave, SolverParams params)
{
    Weave test = weave;
    test.fs->vectorFields.setRandom();

    params.handleWeight = 0.;
    params.lambdacompat = 0.; // weight of compatibility term
    params.lambdareg = 0.000001;    // Tilhonov regularization
    params.curlLambda = 0.; // Weight on the curl component of v
    params.unitLambda = 0.; // Smoothness of s|v|
    params.unitAreaLambda = 1.; // reprojection of area

    Eigen::VectorXd orig;
    GNEnergy(test, params, orig);
    Eigen::SparseMatrix<double> J;
    GNGradient(test, params, J);

    std::ofstream ofs("dump.txt");


    for (int i = 0; i < test.fs->vectorFields.size(); i++)
    {
        Weave perturb = test;
        perturb.fs->vectorFields[i] += 1e-6;
        Eigen::VectorXd newE;
        GNEnergy(perturb, params, newE);
        Eigen::VectorXd findiff = (newE - orig) / 1e-6;
        Eigen::VectorXd exact = J.col(i);
        ofs << (findiff - exact).norm() << " / " << exact.norm() << std::endl;
    }
}


 void projectGradient(Weave &weave, SolverParams params, Eigen::VectorXd &gradient) 
 {
    int nfaces = weave.fs->nFaces();
    int nedges = weave.fs->nEdges();
    int nhandles = weave.nHandles();
    int m = weave.fs->nFields();

    Eigen::VectorXd constraint; 
    constraint.resize(weave.fs->vectorFields.size());
    constraint.setZero();

    for (int e = 0; e < nedges; e++)
    {
        if(weave.fs->data().E(e,0) == -1 || weave.fs->data().E(e,1) == -1)
            continue;
        for (int i = 0; i < m; i++)
        {
                int f = weave.fs->data().E(e, 0);
                int g = weave.fs->data().E(e, 1);
                Eigen::Vector3d edge = weave.fs->data().V.row(weave.fs->data().edgeVerts(e, 0)) - 
                                           weave.fs->data().V.row(weave.fs->data().edgeVerts(e, 1));
                Eigen::Vector2d vif = weave.fs->v(f, i);
                Eigen::Matrix<double, 3, 2> B_f = weave.fs->data().Bs[f];
                double n_v = (B_f * vif).norm();
                Eigen::Vector2d dE = (B_f * vif).transpose() * B_f * weave.fs->sval(f, i) / n_v;

                for (int k = 0; k < 2; k++)
                {
                    constraint( 2 * weave.fs->vidx(f, i) + k ) = dE[k] * params.edgeWeights(e);
                }
                 
                Eigen::Vector2d vperm(0, 0);
                Eigen::MatrixXi permut = weave.fs->Ps(e);
                for (int field = 0; field < m; field++)
                {
                    vperm += permut(i, field) * weave.fs->v(g, field);  
                }
                Eigen::Matrix<double, 3, 2> B_g = weave.fs->data().Bs[g];
                double n_vperm = (B_g * (vperm)).norm();
                Eigen::Vector2d dE_g = (B_g * vperm).transpose() * B_g / n_vperm;
                for (int field = 0; field < m; field++)
                {
                    Eigen::Vector2d dE = dE_g; 
                    dE *= permut(i, field) * weave.fs->sval(g, field);  
                    for (int k = 0; k < 2; k++)
                    {
                        constraint( 2 * weave.fs->vidx(g, field) + k ) = dE[k] * params.edgeWeights(e);
                    }
                }
        }
    }

    constraint.normalize();
    std::cout << "grad size: " << gradient.size() << " constraint size " << constraint.size() << std::endl;
    gradient = gradient - constraint * constraint.transpose() * gradient;

 }
// {
//     int nfaces = weave.fs->nFaces();
//     int m = weave.fs->nFields();
//     double tot_product = 0.;
//     for (int f = 0; f < nfaces; f++)
//     {
//         for (int field = 0; field < m; field++)
//         {
//             Eigen::Vector2d update = weave.fs->v(f, field) + gradient.segment<2>(2 * (f * m) + 2 * field);
//             tot_product += weave.fs->sval(f, field) * update.norm();
//         }
//     }
//     double s = gradient.size(); 

//     gradient += Eigen::VectorXd::Constant(s, (s - tot_product) / s );
// }

int counter = 0;

void computeSMatrix(Weave &weave, SolverParams params, 
                    Eigen::SparseMatrix<double> &newS, 
                    Eigen::SparseMatrix<double> &newSCurl)
{
    int nfaces = weave.fs->nFaces();
    int m = weave.fs->nFields();

    std::vector<Eigen::Triplet<double> > newScoeffs;
    std::vector<Eigen::Triplet<double> > newSCurlcoeffs;

  //  std::cout << "update S" << std::endl;


    for (int e = 0; e < weave.fs->nEdges(); e++)
    {
        if(weave.fs->data().E(e,0) == -1 || weave.fs->data().E(e,1) == -1)
            continue;
        for (int i = 0; i < m; i++)
        {
                int f = weave.fs->data().E(e, 0);
                int g = weave.fs->data().E(e, 1);
                Eigen::Vector3d edge = weave.fs->data().V.row(weave.fs->data().edgeVerts(e, 0)) - 
                                            weave.fs->data().V.row(weave.fs->data().edgeVerts(e, 1));
                Eigen::Vector2d vif = weave.fs->v(f, i);
                double n_v = (weave.fs->data().Bs[f] * vif).norm();
             //   vif *= weave.fs->sval(f, i);
                 
                Eigen::Vector2d vperm(0, 0);
                Eigen::MatrixXi permut = weave.fs->Ps(e);
                int adj_field = -1;
                for (int field = 0; field < m; field++)
                {
                    vperm += permut(i, field) * weave.fs->v(g, field); 
                    if(permut(i, field) != 0)
                        adj_field = field;
                }
                Eigen::Matrix<double, 3, 2> B_f = weave.fs->data().Bs[f];
                Eigen::Matrix<double, 3, 2> B_g = weave.fs->data().Bs[g];
                double n_vperm = (B_g * (vperm)).norm();

                // for (int field = 0; field < m; field++)
                // {
                //     if (permut(i, field) != 0)
                //         vperm *= weave.fs->sval(g, field);  
                // }

                newScoeffs.push_back(Triplet<double>(e, m*f + i, params.unitLambda * n_v));
                newScoeffs.push_back(Triplet<double>(e, m*g + adj_field, -params.unitLambda * n_vperm));
        
                newSCurlcoeffs.push_back(Triplet<double>(e, m*f + i, params.curlLambda * (B_f * vif).dot(edge) ));
                newSCurlcoeffs.push_back(Triplet<double>(e, m*g + adj_field, -params.curlLambda * (B_g * vperm).dot(edge) ));
        

        }
    }

    newS.resize(weave.fs->nEdges(), m*nfaces);
    newS.setZero();
    newS.setFromTriplets(newScoeffs.begin(), newScoeffs.end());

    newSCurl.resize(weave.fs->nEdges(), m*nfaces);
    newSCurl.setZero();
    newSCurl.setFromTriplets(newSCurlcoeffs.begin(), newSCurlcoeffs.end());
}

void oneStep(Weave &weave, SolverParams params)
{    
    int nvars = weave.fs->vectorFields.size();
    Eigen::VectorXd r;
  //  params.handleWeight = 1.;
    // params.lambdacompat = 0.;
    // params.curlLambda = 0.;
    // params.unitLambda = 0.;
    // params.unitAreaLambda = 1.;


    GNEnergy(weave, params, r);
    Eigen::SparseMatrix<double> M;
    GNmetric(weave, M);

    int nfaces = weave.fs->nFaces();
    int m = weave.fs->nFields();

    Eigen::SparseMatrix<double> oldS;
    Eigen::SparseMatrix<double> oldSCurl;
    computeSMatrix(weave, params, oldS, oldSCurl);
    Eigen::VectorXd curSVals = weave.fs->vectorFields.segment(5*nfaces*m, nfaces);
    double sEnergy = curSVals.transpose() * Eigen::MatrixXd(oldS).transpose() * Eigen::MatrixXd(oldS) * curSVals; 

    std::cout << "original energy: " << 0.5 * r.transpose() * M * r + sEnergy << std::endl;
  //  std::cout << "Building matrix" << std::endl;
    Eigen::SparseMatrix<double> J;
    GNGradient(weave, params, J);
    Eigen::SparseMatrix<double> optMat(nvars, nvars);
    std::vector<Eigen::Triplet<double> > coeffs;

    for(int i=0; i<nfaces; i++)
    {
        for(int j=0; j<m; j++)
        {
            for(int k=0; k<3; k++)
            {
                int idx = 2*nfaces*m + 3*m*i + 3*j + k;            
                double facearea = weave.fs->faceArea(i);
                coeffs.push_back(Eigen::Triplet<double>(idx, idx, params.lambdareg*facearea));
                int shift = 5*nfaces*m + i + j;
                coeffs.push_back(Eigen::Triplet<double>(shift, shift, params.lambdareg*facearea));
            }
        }
    }
    

 //   AAAAAAARHRHRRHGGGG
    for(int i=0; i<nfaces * m * 2; i++)
    {
  //      int shift = 5*nfaces*m + j;
  //      coeffs.push_back(Eigen::Triplet<double>(shift, shift, params.lambdareg*facearea));
        coeffs.push_back(Eigen::Triplet<double>(i, i, params.lambdareg * .001));
    }

    


    // for(int i=0; i<nfaces; i++)
    // {
    //     for(int j=0; j<m; j++)
    //     {
    //         int shift = 5*nfaces*m + j;
    //         double facearea = weave.fs->faceArea(i);
    //         coeffs.push_back(Eigen::Triplet<double>(i+shift, i+shift, params.lambdareg));
    //     }
    // }


    optMat.setFromTriplets(coeffs.begin(), coeffs.end());
    optMat += J.transpose() * M * J;
 //   std::cout << "Done, " << optMat.nonZeros() << " nonzeros" << std::endl;
    optMat.makeCompressed();
    Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double> > solver;
 //   std::cout << "Analyzing" << std::endl;
    solver.analyzePattern(optMat);
    Eigen::VectorXd rhs = J.transpose() * M * r;
    std::cout << "Solving" << std::endl;
    solver.factorize(optMat);
    Eigen::VectorXd update = solver.solve(rhs);
    projectGradient(weave, params, update);

    double t = lineSearch(weave, params, sEnergy, true, update);
    
    GNEnergy(weave, params, r);
    std::cout << "Done, new energy: " << 0.5 * r.transpose()*M*r + sEnergy<< std::endl;
 //   weave.fs->normalizeFields();

    GNEnergy(weave, params, r);
    std::cout << "Done, rescaled new energy: " << 0.5 * r.transpose()*M*r + sEnergy << std::endl;
    m = 1;
    
    Eigen::SparseMatrix<double> newS;
    Eigen::SparseMatrix<double> newSCurl;
    computeSMatrix(weave, params, newS, newSCurl);

    Eigen::SparseMatrix<double> iden;
    iden.resize(m*nfaces, m*nfaces);
    iden.setIdentity();
    double scale = .00000001;

  //  std::cout << "make matrix"  << std::endl;
  //  Eigen::SparseMatrix<double> STS = newS.transpose() * newS + 0.5 * newSCurl.transpose() * newSCurl + scale * iden;


    // Eigen::VectorXd largestEigenVec;
    // largestEigenVec.resize(m*nfaces);
    // largestEigenVec.setRandom();

    // for (int i = 0; i < 1000; i++)
    // {
    //     largestEigenVec = STS * largestEigenVec;
    //     largestEigenVec.normalize();
    // } 

    // double largestEigenVal = largestEigenVec.dot(STS * largestEigenVec); // this gets the norm and sign.

    // Eigen::MatrixXd STS_shift = STS - largestEigenVal * iden;

    // largestEigenVec.resize(m*nfaces);
    // largestEigenVec.setRandom();

    // for (int i = 0; i < 1000; i++)
    // {
    //     largestEigenVec = STS_shift * largestEigenVec;
    //     largestEigenVec.normalize();
    // } 

    // curSVals = largestEigenVec;
    Eigen::MatrixXd STS =    Eigen::MatrixXd(newS).transpose() * Eigen::MatrixXd(newS)
                             + 0.5 * Eigen::MatrixXd(newSCurl).transpose() * Eigen::MatrixXd(newSCurl) 
                             + scale * Eigen::MatrixXd(iden);


    SelfAdjointEigenSolver<MatrixXd> eigensolver(STS);
    curSVals = eigensolver.eigenvectors().col(params.eigenvector).real() * nfaces;

    for (int i = 0; i < m; i++) // DOUBLE CHECK FOR M > 1, probably will be buggy.
    {
        weave.fs->vectorFields.segment(5*nfaces*m + nfaces * i, nfaces) = curSVals;
    }

    //    std::cout << "The eigenvectors of A are:\n" << curSVals << std::endl;

    GNEnergy(weave, params, r);
    curSVals = weave.fs->vectorFields.segment(5*nfaces*m, nfaces);
    sEnergy = curSVals.transpose() * Eigen::MatrixXd(newS).transpose() * Eigen::MatrixXd(newS) * curSVals; 
    std::cout << "Account for S smoothness ***: " << 0.5 * r.transpose()*M*r + sEnergy << std::endl;

    // // counter++;
    // if (t <  0.0000625)
    // {
 //      GNtestFiniteDifferences(weave, params);
   //    exit(-1);
    // }
}

double lineSearch(Weave &weave, SolverParams params, double shiftEnergy, bool toProject, const Eigen::VectorXd &update)
{
    double t = 1.0;
    double c1 = 0.1;
    double c2 = 0.9;
    double alpha = 0;
    double infinity = 1e6;
    double beta = infinity;

    int nvars = weave.fs->vectorFields.size();
    Eigen::VectorXd r;
    GNEnergy(weave, params, r);
    Eigen::SparseMatrix<double> M;
    GNmetric(weave, M);
    Eigen::SparseMatrix<double> J;
    GNGradient(weave, params, J);

    VectorXd dE;
    VectorXd newdE;
    VectorXd startVF = weave.fs->vectorFields;
    
    double orig = 0.5 * r.transpose() * M * r + shiftEnergy;
    dE = J.transpose() * M * r;
    if (toProject)
        projectGradient(weave, params, dE);
    double deriv = -dE.dot(update);
    assert(deriv < 0);
    
  //  std::cout << "Starting line search, original energy " << orig << ", descent magnitude " << deriv << std::endl;

    while (true)
    {
        weave.fs->vectorFields = startVF - t * update;
        GNEnergy(weave, params, r);
        double newenergy = 0.5 * r.transpose() * M * r + shiftEnergy;
        GNGradient(weave, params, J);
        newdE = J.transpose() * M * r;
        if (toProject)
            projectGradient(weave, params, newdE);

   //     std::cout << "Trying t = " << t << ", energy now " << newenergy<< std::endl;
        
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


/*void firstStep(Weave &weave, SolverParams params)
{    
    int nvars = weave.fs->vectorFields.size();
    Eigen::VectorXd r;
    params.handleWeight = 100.;
    double lambdacompat = params.lambdacompat;
  //  params.lambdacompat = .1;
    double curlreg = params.curlreg;
    params.curlreg = 0;

    int nfaces = weave.fs->nFaces();
    int m = weave.fs->nFields();

  //  weave.fs->vectorFields.segment(2 * nfaces*m, 3 * nfaces*m).setZero();
    weave.fs->vectorFields.segment(5*nfaces*m, nfaces*m) = Eigen::VectorXd::Constant(nfaces * m, 1.);

    GNEnergy(weave, params, r);
    Eigen::SparseMatrix<double> M;
    GNmetric(weave, M);



    Eigen::SparseMatrix<double> oldS;
    computeSMatrix(weave, params, oldS);
    Eigen::VectorXd curSVals = weave.fs->vectorFields.segment(5*nfaces*m, nfaces);
    double sEnergy = curSVals.transpose() * Eigen::MatrixXd(oldS).transpose() * Eigen::MatrixXd(oldS) * curSVals; 

 //   std::cout << "original energy: " << 0.5 * r.transpose() * M * r + sEnergy << std::endl;
  //  std::cout << "Building matrix" << std::endl;
    Eigen::SparseMatrix<double> J;
    GNGradient(weave, params, J);
    Eigen::SparseMatrix<double> optMat(nvars, nvars);
    std::vector<Eigen::Triplet<double> > coeffs;

    for(int i=0; i<nfaces; i++)
    {
        for(int j=0; j<m; j++)
        {
            for(int k=0; k<3; k++)
            {
                int idx = 2*nfaces*m + 3*m*i + 3*j + k;            
                double facearea = weave.fs->faceArea(i);
                coeffs.push_back(Eigen::Triplet<double>(idx, idx, params.lambdareg*facearea));
                int shift = 5*nfaces*m + i + j;
                coeffs.push_back(Eigen::Triplet<double>(shift, shift, params.lambdareg*facearea));
            }
        }
    }

    for(int i=0; i<nfaces * m * 2; i++)
    {
        coeffs.push_back(Eigen::Triplet<double>(i, i, params.lambdareg * .001));
    }

    optMat.setFromTriplets(coeffs.begin(), coeffs.end());
    optMat += J.transpose() * M * J;
 //   std::cout << "Done, " << optMat.nonZeros() << " nonzeros" << std::endl;
    optMat.makeCompressed();
    Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double> > solver;
 //   std::cout << "Analyzing" << std::endl;
    solver.analyzePattern(optMat);
    Eigen::VectorXd rhs = J.transpose() * M * r;
    std::cout << "Solving" << std::endl;
    solver.factorize(optMat);
    Eigen::VectorXd update = solver.solve(rhs);

    double t = lineSearch(weave, params, sEnergy, false, update);
    
    GNEnergy(weave, params, r);
//    std::cout << "Done, new energy: " << 0.5 * r.transpose()*M*r + sEnergy<< std::endl;


    params.handleWeight = 0.;
    params.lambdacompat = 0;
    params.curlreg = curlreg;

    GNEnergy(weave, params, r);
    std::cout << "Done, rescaled new energy: " << 0.5 * r.transpose()*M*r + sEnergy << std::endl;
    m = 1;
    
}*/