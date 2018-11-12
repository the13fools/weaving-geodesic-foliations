#include "LinearSolver.h"
#include <map>
#include <vector>
#include <set>
#include <cassert>
#include <Eigen/Geometry>
#include <Eigen/SPQRSupport>
#include <iostream>

#include "GaussNewton.h"
#include "Weave.h"
#include "Surface.h"


void LinearSolver::addHandle(const Handle &h)
{
    handles.push_back(h);
}

void LinearSolver::clearHandles()
{
    handles.clear();
}

// int LinearSolver::numPrimalDOFs()
// {
//     return 2 * F.rows();
// }

// int LinearSolver::numDualDOFs()
// {
//     return 2 * F.rows();
// }

// void LinearSolver::generateRandomField(Eigen::VectorXd &primalVars, Eigen::VectorXd &dualVars)
// {
//     primalVars.resize(numPrimalDOFs());
//     int nfaces = F.rows();
//     Eigen::SparseMatrix<double> proj;
//     unconstrainedProjection(proj);
//     int nhandles = handles.size();
    
//     Eigen::VectorXd handlevals(2 * nfaces);
//     handlevals.setZero();
//     for (int i = 0; i < nhandles; i++)
//     {
//         handlevals.segment<2>(2 * handles[i].face) = handles[i].dir / handles[i].dir.norm() * 100;
//     }

//     Eigen::VectorXd interp(2 * nfaces - 2 * nhandles);
//     interp.setRandom();
//     primalVars = proj * interp + handlevals;
//     for (int i = 0; i < nfaces; i++)
//     {
//         primalVars.segment<2>(2 * i).normalize();
//     }

//     dualVars.resize(numDualDOFs());
//     dualVars.setZero();    
// }

// void LinearSolver::generateHarmonicField(Eigen::VectorXd &primalVars, Eigen::VectorXd &dualVars)
// {
//     primalVars.resize(numPrimalDOFs());
//     int nfaces = F.rows();
//     Eigen::SparseMatrix<double> proj;
//     unconstrainedProjection(proj);
//     int nhandles = handles.size();
//     Eigen::VectorXd handlevals(2 * nfaces);
//     handlevals.setZero();
//     for (int i = 0; i < nhandles; i++)
//     {
//         handlevals.segment<2>(2 * handles[i].face) = handles[i].dir / handles[i].dir.norm() * 100;
//     }

//     // min ||D (proj^T v + hvals) ||^2
//     // proj * D^T D * proj^T v + proj * D^T D * hvals = 0
//     Eigen::SparseMatrix<double> D;
//     differentialOperator(D);
//     Eigen::VectorXd rhs = -proj.transpose() * D.transpose() * D * handlevals;
//     Eigen::SparseMatrix<double> M = proj.transpose() * D.transpose() * D * proj;
//     Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver(M);
//     Eigen::VectorXd interp = solver.solve(rhs);
//     primalVars = proj * interp + handlevals;
//     for (int i = 0; i < nfaces; i++)
//     {
//         primalVars.segment<2>(2 * i).normalize();
//     }

//     dualVars.resize(numDualDOFs());
//     dualVars.setZero();    
// }

// TODO ********************
// *************************
// port over the right energy
// void LinearSolver::setFaceEnergies(const Weave &weave, const Eigen::VectorXd &primalVars, const Eigen::VectorXd &dualVars, Eigen::VectorXd &faceEnergies)
// {
//     faceEnergies.resize(weave.fs->data().F.rows());
//     faceEnergies.setZero();

//     int nedges = weave.fs->data().E.rows();
//     for (int i = 0; i < nedges; i++)
//     {
//         int face1 = weave.fs->data().E(i, 2);
//         int face2 = weave.fs->data().E(i, 3);
//         Eigen::Vector2d edgevec;
//         edgevec[0] = V(E(i, 1), 0) - V(E(i, 0), 0);
//         edgevec[1] = V(E(i, 1), 1) - V(E(i, 0), 1);
//         edgevec /= edgevec.norm();
        
//         double energy = (primalVars.segment<2>(2 * face2).dot(edgevec) - primalVars.segment<2>(2 * face1).dot(edgevec));
//         faceEnergies[face1] += energy * energy;
//         faceEnergies[face2] += energy * energy;
//     }
// }

static void identityMatrix(int k, int n, Eigen::SparseMatrix<double> &I)
{
    std::vector<Eigen::Triplet<double> > coeffs;
    for (int i = 0; i < k; i++)
    {
        coeffs.push_back(Eigen::Triplet<double>(i, i, 1.0));
    }
    I.resize(n, n);
    I.setFromTriplets(coeffs.begin(), coeffs.end());
}

static void massMatrix(const Weave &weave, Eigen::SparseMatrix<double> &M)
{
    std::vector<Eigen::Triplet<double> > Mcoeffs;
    int nfaces = weave.fs->data().F.rows();
    int m = weave.fs->nFields();

    M.resize(2*m*nfaces, 2*m*nfaces);
    M.setZero();

    int term = 0;
    for (int f = 0; f < nfaces; f++)
    {
        for (int i = 0; i < m; i++)
        {
            double area = weave.fs->faceArea(f);
            Eigen::Matrix2d BTB = weave.fs->data().Bs[f].transpose() * weave.fs->data().Bs[f];
            for (int j = 0; j < 2; j++)
            {
                for (int k = 0; k < 2; k++)
                {
                    Mcoeffs.push_back(Eigen::Triplet<double>(term + j, term + k, area*BTB(j, k)));
                }
            }
            term += 2;
        }
    }

    M.setFromTriplets(Mcoeffs.begin(), Mcoeffs.end());
}

void LinearSolver::computeEnergy(const Weave &weave, SolverParams params, const Eigen::VectorXd &primalVars, const Eigen::VectorXd &dualVars )
{
    int nfaces = weave.fs->data().F.rows();
    int m = weave.fs->nFields();

    Eigen::SparseMatrix<double> D;
    Eigen::SparseMatrix<double> P;
    Eigen::SparseMatrix<double> BTB;

    differentialOperator(weave, params, D);
    unconstrainedProjection(weave, P);
    massMatrix(weave, BTB);

    std::cout << " The current energy is " << dualVars.transpose() * BTB * dualVars << " + " 
              << (primalVars + dualVars).transpose() * D.transpose() * D * (primalVars + dualVars)  << " = " 
              << (primalVars + dualVars).transpose() * D.transpose() * D * (primalVars + dualVars) + dualVars.transpose() * BTB * dualVars << std::endl;
    // Eigen::SparseMatrix<double> I;
    // differentialOperator(weave, params, D);
}

void LinearSolver::updatePrimalVars(const Weave &weave, SolverParams params, Eigen::VectorXd &primalVars, Eigen::VectorXd &dualVars )
{
    int nfaces = weave.fs->data().F.rows();
    int m = weave.fs->nFields();
    std::cout << " pre-primal update ";
    computeEnergy(weave, params, primalVars, dualVars);
    
    for (int i = 0; i < nfaces * m; i++)
    {
        primalVars.segment<2>(2 * i) = primalVars.segment<2>(2 * i) + dualVars.segment<2>(2 * i);
        dualVars.segment<2>(2 * i) = primalVars.segment<2>(2 * i);
    }

    // Eigen::SparseMatrix<double> D;
    Eigen::SparseMatrix<double> P;
    // Eigen::SparseMatrix<double> I;
    // differentialOperator(weave, params, D);
    unconstrainedProjection(weave, P);
    int nhandles = handles.size();
    // identityMatrix(2 * nfaces * m, I);
    // double t = params.lambdacompat;
    // Eigen::SparseMatrix<double> op = I + t * D.transpose() * D;
    // Eigen::SparseMatrix<double> M = P.transpose() * op * P;
    // Eigen::SPQR<Eigen::SparseMatrix<double> > solver(M);

    // Eigen::VectorXd handlevals(2 * nfaces * m);
    // handlevals.setZero();
    // for (int i = 0; i < nhandles; i++)
    // {
    //     handlevals.segment<2>(2 * (handles[i].face * m + handles[i].field ) ) = handles[i].dir / handles[i].dir.norm() * 1.;
    // }
    // // Eigen::VectorXd rhs = P.transpose() * (primalVars - op * handlevals);

    // // Eigen::VectorXd newprim = solver.solve(rhs);
    // primalVars = P * primalVars + handlevals;

    for (int i = 0; i < nhandles; i++)
    {
        primalVars.segment<2>(2 * (handles[i].face * m + handles[i].field ) ) = handles[i].dir / handles[i].dir.norm();
    }

    for (int i = 0; i < nfaces; i++)
    {
        for ( int cover = 0; cover < m; cover++)
        {
            Eigen::Matrix<double, 3, 2> B_f = weave.fs->data().Bs[i];
            Eigen::Vector3d ambient = (B_f * primalVars.segment<2>(2 * i * m + 2 * cover));
            primalVars.segment<2>(2 * i * m + 2 * cover) /= ambient.norm();
            dualVars.segment<2>(2 * i * m + 2 * cover) -= primalVars.segment<2>(2 * i * m + 2 * cover);
        }
    }


    std::cout << "post-primal update ";
    computeEnergy(weave, params, primalVars, dualVars);

}

int LinearSolver::dualMatrixSize(const Weave &weave)
{
    int nfaces = weave.fs->data().F.rows();
    int m = weave.fs->nFields();
    int intedges = weave.fs->numInteriorEdges();
    int nhandles = handles.size();

    return 2*nfaces*m - 2*nhandles + intedges*m;
}

void LinearSolver::buildDualMatrix(const Weave &weave, SolverParams params, Eigen::VectorXd &primalVars, Eigen::VectorXd &dualVars)
{
    // min_delta, \lambda   0.5 delta^2 + \lambda^T L (v + delta)
    // delta + L^T \lambda = 0
    // L delta + LL^T \lambda = 0
    // -L v + LL^T \lambda = 0
    // LL^T \lambda = Lv
    // delta = - L^T \lambda

    Eigen::SparseMatrix<double> D;
    Eigen::SparseMatrix<double> P;
    Eigen::SparseMatrix<double> C;
    Eigen::SparseMatrix<double> BTB;

    differentialOperator(weave, params, D);
    unconstrainedProjection(weave, P);

    curlOperator(weave, params, C);

    int nfaces = weave.fs->data().F.rows();
    int m = weave.fs->nFields();
    int intedges = weave.fs->numInteriorEdges();
    int nhandles = handles.size();

    int matrixSize = dualMatrixSize(weave);

    massMatrix(weave, BTB);

    double t = params.lambdacompat;
    Eigen::SparseMatrix<double> op = BTB + t * D.transpose() * D;
    Eigen::SparseMatrix<double> M = P.transpose() * op * P; // this might be wrong w handles

    std::vector<Eigen::Triplet<double> > dualCoeffs;

    for (int k=0; k<M.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(M,k); it; ++it)
        {
            dualCoeffs.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
        }
    }

    Eigen::SparseMatrix<double> CP = (C*P);
    for (int k=0; k<CP.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(CP,k); it; ++it)
        {
            dualCoeffs.push_back(Eigen::Triplet<double>(it.row() + 2*m*nfaces-2*nhandles, it.col(), it.value()));
        }
    }

    Eigen::SparseMatrix<double> CPT = Eigen::SparseMatrix<double>((C*P).transpose());
    for (int k=0; k<CPT.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(CPT,k); it; ++it)
        {
            dualCoeffs.push_back(Eigen::Triplet<double>(it.row(), it.col() + 2*m*nfaces-2*nhandles, it.value()));
        }
    }
    Eigen::SparseMatrix<double> dualMat;
    dualMat.resize(matrixSize, matrixSize);
    dualMat.setFromTriplets(dualCoeffs.begin(), dualCoeffs.end());

    std::cout << matrixSize << " matrix size " << std::endl;

    Eigen::SPQR<Eigen::SparseMatrix<double> > solver(dualMat);
    dualSolver = &solver;
 //   dualSolver.analyzePattern(dualMat);   // for this step the numerical values of A are not used
 //   dualSolver.factorize(dualMat);
}

void LinearSolver::updateDualVars_new(const Weave &weave, SolverParams params, Eigen::VectorXd &primalVars, Eigen::VectorXd &dualVars)
{
    // min_delta, \lambda   0.5 delta^2 + \lambda^T L (v + delta)
    // delta + L^T \lambda = 0
    // L delta + LL^T \lambda = 0
    // -L v + LL^T \lambda = 0
    // LL^T \lambda = Lv
    // delta = - L^T \lambda

    Eigen::SparseMatrix<double> D;
    Eigen::SparseMatrix<double> P;
  //  Eigen::SparseMatrix<double> I;
    Eigen::SparseMatrix<double> curlOp;
    Eigen::SparseMatrix<double> BTB;

    differentialOperator(weave, params, D);
    unconstrainedProjection(weave, P);

    curlOperator(weave, params, curlOp);

    int nfaces = weave.fs->data().F.rows();
    int m = weave.fs->nFields();
    int nhandles = handles.size();
    int intedges = weave.fs->numInteriorEdges();

    int matrixSize = dualMatrixSize(weave);

  // //  identityMatrix(2*nfaces*m, 2*nfaces*m + intedges*m, I);
    massMatrix(weave, BTB);
    double t = params.lambdacompat;
    Eigen::SparseMatrix<double> op = BTB + t * D.transpose() * D;
    Eigen::SparseMatrix<double> M = P.transpose() * op * P; // this might be wrong w handles

  //   Eigen::SparseMatrix<double> Full;
  //   Full.resize(2  * nfaces * m + intedges * m, 2 * m * nfaces + intedges * m);
  //   Full += op;
  //   Full += curlOp;
  //   Full += Eigen::SparseMatrix<double>(curlOp.transpose());

  //   Eigen::SPQR<Eigen::SparseMatrix<double> > solver(Full);

    Eigen::VectorXd rhs(matrixSize);
    rhs.segment(0, 2*nfaces*m - 2*nhandles) = -t * P.transpose() * D.transpose() * D * primalVars;
    rhs.segment(2*nfaces*m - 2*nhandles, intedges * m) = -curlOp * primalVars;

    // Eigen::VectorXd rhs(matrixSize);
    // rhs.setZero();
    // rhs -= t * P.transpose() * D.transpose() * D * top_primal;
    // rhs -= curlOp * primalVars;

    std::cout << rhs.size() << " rhs size " << std::endl;

/////////////////////**********************************************************************************///////////////////////
    Eigen::SparseMatrix<double> C;

    curlOperator(weave, params, C);



    massMatrix(weave, BTB);

    std::vector<Eigen::Triplet<double> > dualCoeffs;

    for (int k=0; k<M.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(M,k); it; ++it)
        {
            dualCoeffs.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
        }
    }

    Eigen::SparseMatrix<double> CP = (C*P);
    for (int k=0; k<CP.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(CP,k); it; ++it)
        {
            dualCoeffs.push_back(Eigen::Triplet<double>(it.row() + 2*m*nfaces-2*nhandles, it.col(), it.value()));
        }
    }

    Eigen::SparseMatrix<double> CPT = Eigen::SparseMatrix<double>((C*P).transpose());
    for (int k=0; k<CPT.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(CPT,k); it; ++it)
        {
            dualCoeffs.push_back(Eigen::Triplet<double>(it.row(), it.col() + 2*m*nfaces-2*nhandles, it.value()));
        }
    }
    Eigen::SparseMatrix<double> dualMat;
    dualMat.resize(matrixSize, matrixSize);
    dualMat.setFromTriplets(dualCoeffs.begin(), dualCoeffs.end());

    std::cout << matrixSize << " matrix size " << std::endl;

    Eigen::SPQR<Eigen::SparseMatrix<double> > solver(dualMat);

/////////////////////**********************************************************************************///////////////////////

    Eigen::VectorXd deltalambda = solver.solve(rhs);

    dualVars = P * deltalambda.segment(0, 2  * nfaces * m - 2*nhandles);
   


    // std::vector<Eigen::Triplet<double> > regcoeffs;
    // double reg = 1e-8;
    // for (int i = 0; i < intedges * m; i++)
    // {
    //     regcoeffs.push_back(Eigen::Triplet<double>(i, i, reg));
    // }
    // Eigen::SparseMatrix<double> Reg(intedges * m, intedges * m);
    // Reg.setFromTriplets(regcoeffs.begin(), regcoeffs.end());

    // Eigen::SparseMatrix<double> M = curlOp * P * P.transpose() * curlOp.transpose() + Reg;
    // Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver(M);
    
    // Eigen::VectorXd rhs(intedges * m);
    // rhs = curlOp * primalVars;
    
    // Eigen::VectorXd lambda = solver.solve(rhs);

    // dualVars = -P * P.transpose() * curlOp.transpose() * lambda;
    
    std::cout << "  post-dual update ";
    computeEnergy(weave, params, primalVars, dualVars);
    std::cout << "Dual vars now " << dualVars.norm() << " geodesic residual " << (curlOp * (primalVars + dualVars)).norm() << std::endl;
    // if ( std::isnan( dualVars.norm() ) )
}



// TODO
// *************************
void LinearSolver::curlOperator(const Weave &weave, SolverParams params, Eigen::SparseMatrix<double> &curlOp)
{
    std::vector<Eigen::Triplet<double> > coeffs;
    int nedges = weave.fs->data().E.rows();
    int intedges = weave.fs->numInteriorEdges();
    int m = weave.fs->nFields();
    int nfaces = weave.fs->data().F.rows();
    int nhandles = handles.size();

    int term = 0;

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
                edge.normalize();
                 
                Eigen::Vector2d vperm(0, 0);
                Eigen::MatrixXi permut = weave.fs->Ps(e);

                int adj_field = -1;

                for (int field = 0; field < m; field++)
                {
                    vperm += permut(i, field) * weave.fs->v(g, field);
                    if (permut(i, field) != 0) 
                        adj_field = field;
                }
                Eigen::Matrix<double, 3, 2> B_f = weave.fs->data().Bs[f];
                Eigen::Matrix<double, 3, 2> B_g = weave.fs->data().Bs[g];

                if (params.edgeWeights(e) > 0.) 
                {
                    for (int j = 0; j < 2; j++)
                    {
                        coeffs.push_back(Eigen::Triplet<double>(term, weave.fs->vidx(f, i) + j, (B_f.transpose() * edge)[j]));
                        coeffs.push_back(Eigen::Triplet<double>(term, weave.fs->vidx(g, adj_field) + j, -(B_g.transpose() * edge)[j]));
                    }
                }
                term++;
        }
    }


    curlOp.resize(intedges*m, 2 * m * nfaces);
    curlOp.setFromTriplets(coeffs.begin(), coeffs.end());
}

// TODO: Steal this from prev version.
// *************************
void LinearSolver::differentialOperator(const Weave &weave, SolverParams params, Eigen::SparseMatrix<double> &D)
{
    std::vector<Eigen::Triplet<double> > coeffs;
    int nedges = weave.fs->data().E.rows();
    int intedges = weave.fs->numInteriorEdges();
    int m = weave.fs->nFields();

    int term = 0;

        // compatibility constraint
    for (int e = 0; e < nedges; e++)
    {
        if(weave.fs->data().E(e,0) == -1 || weave.fs->data().E(e,1) == -1)
            continue;
        for (int i = 0; i < m; i++)
        {
            for (int side = 0; side < 1; side++)
            {
                int f = (side == 0 ? weave.fs->data().E(e, 0) : weave.fs->data().E(e, 1));
                int g = (side == 0 ? weave.fs->data().E(e, 1) : weave.fs->data().E(e, 0));
                Eigen::Vector2d vif = weave.fs->v(f, i);
                Eigen::Vector2d vpermut(0, 0);
                Eigen::MatrixXi permut = weave.fs->Ps(e);
                if (side == 1)
                    permut.transposeInPlace();

                int adj_field = -1;
                for (int field = 0; field < m; field++)
                {
                    vpermut += permut(i, field) * weave.fs->v(g, field);
                    if (permut(i, field) != 0) 
                        adj_field = field;
                }

                Eigen::Matrix<double, 3, 2> B_f = weave.fs->data().Bs[f];
                Eigen::Matrix<double, 3, 2> B_g = weave.fs->data().Bs[g];

                for (int coeff = 0; coeff < 3; coeff++)
                {
                    if (params.edgeWeights(e) > 0.) 
                    {
                        for (int k = 0; k < 2; k++)
                        {
                            coeffs.push_back(Eigen::Triplet<double>(term, weave.fs->vidx(f, i) + k, B_f(coeff, k)));
                            coeffs.push_back(Eigen::Triplet<double>(term, weave.fs->vidx(g, adj_field) + k, -B_g(coeff, k)));
                        }
                    }
                    term++;
                }
            }
        }
    }

    // for (int e = 0; e < weave.fs->nEdges(); e++)
    // {
    //     if(weave.fs->data().E(e,0) == -1 || weave.fs->data().E(e,1) == -1)
    //         continue;
    //     for (int i = 0; i < m; i++)
    //     {
    //         int f = weave.fs->data().E(e, 0);
    //         int g = weave.fs->data().E(e, 1);
    //         for (int j = 0; j < 2; j++)
    //         {
    //             coeffs.push_back(Eigen::Triplet<double>(2 * e + j, 2 * f + j, 1.0));
    //             coeffs.push_back(Eigen::Triplet<double>(2 * e + j, 2 * g + j, -1.0));
    //         }
    //     }
    // }
    int nfaces = weave.fs->data().F.rows();
    D.resize(3 * intedges * m, 2*nfaces*m);
    D.setFromTriplets(coeffs.begin(), coeffs.end());
}

void LinearSolver::updateDualVars(const Weave &weave, SolverParams params, Eigen::VectorXd &primalVars, Eigen::VectorXd &dualVars)
{
    // min_delta, \lambda   0.5 delta^2 + \lambda^T L (v + delta)
    // delta + L^T \lambda = 0
    // L delta + LL^T \lambda = 0
    // -L v + LL^T \lambda = 0
    // LL^T \lambda = Lv
    // delta = - L^T \lambda
    Eigen::SparseMatrix<double> curlOp;
    curlOperator(weave, params, curlOp);
    int nfaces = weave.fs->data().F.rows();
    int m = weave.fs->nFields();
    
    Eigen::SparseMatrix<double> V;
    unconstrainedProjection(weave, V);
    int intedges = weave.fs->numInteriorEdges();

    std::vector<Eigen::Triplet<double> > regcoeffs;
    double reg = 1e-8;
    for (int i = 0; i < intedges * m; i++)
    {
        regcoeffs.push_back(Eigen::Triplet<double>(i, i, reg));
    }
    Eigen::SparseMatrix<double> Reg(intedges * m, intedges * m);
    Reg.setFromTriplets(regcoeffs.begin(), regcoeffs.end());

    Eigen::SparseMatrix<double> M = curlOp * V * V.transpose() * curlOp.transpose() + Reg;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver(M);
    
    Eigen::VectorXd rhs(intedges * m);
    rhs = curlOp * primalVars;
    
    Eigen::VectorXd lambda = solver.solve(rhs);

    dualVars = -V * V.transpose() * curlOp.transpose() * lambda;
    
    std::cout << "  post-dual update ";
    computeEnergy(weave, params, primalVars, dualVars);
    std::cout << "Dual vars now " << dualVars.norm() << " geodesic residual " << (curlOp * (primalVars + dualVars)).norm() << std::endl;
    // if ( std::isnan( dualVars.norm() ) )
    // {
    //     std::cout << "whoops" << std::endl;
    //     exit(-1);
    // }
}

void LinearSolver::unconstrainedProjection(const Weave &weave, Eigen::SparseMatrix<double> &proj)
{
    int nfaces = weave.fs->data().F.rows();
    int nhandles = handles.size();
    int m = weave.fs->nFields();
    std::set<int> handlefaces;
    for (int i = 0; i < nhandles; i++)
        handlefaces.insert( handles[i].face * m + handles[i].field );
    std::vector<Eigen::Triplet<double> > Vcoeffs;
    int col = 0;
    for (int i = 0; i < nfaces * m; i++)
    {
        if (handlefaces.count(i))
            continue;
        Vcoeffs.push_back(Eigen::Triplet<double>(2 * i, col, 1.0));
        Vcoeffs.push_back(Eigen::Triplet<double>(2 * i + 1, col + 1, 1.0));
        col += 2;
    }
    assert(col == 2*nfaces*m - 2*nhandles);
    proj.resize(2*nfaces*m, 2*nfaces * m-2*nhandles);
    proj.setFromTriplets(Vcoeffs.begin(), Vcoeffs.end());

}