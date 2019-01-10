#include "LinearSolver.h"
#include <map>
#include <vector>
#include <set>
#include <cassert>
#include <Eigen/Geometry>
#include <iostream>

#include "GaussNewton.h"
#include "Weave.h"
#include "Surface.h"



DualSolver::DualSolver(Eigen::SparseMatrix<double> &M)
{
    solver.compute(M);
}

void DualSolver::solve(const Eigen::VectorXd &rhs, Eigen::VectorXd &x)
{
    x = solver.solve(rhs);
}


void LinearSolver::addHandle(const Handle &h)
{
    handles.push_back(h);
}

void LinearSolver::clearHandles()
{
    handles.clear();
}

void LinearSolver::takeSomeSteps(const Weave &weave, SolverParams params, Eigen::VectorXd &primalVars, Eigen::VectorXd &dualVars, bool isRoSy, int numSteps)
{
    DualSolver *ds = buildDualUpdateSolver(weave, params, isRoSy);
    for(int i=0; i<numSteps; i++)
    {
        std::cout << "###############" << std::endl;
        std::cout << "Step " << i+1 << " of " << numSteps << std::endl;
        std::cout << "###############" << std::endl;
        updateDualVars_new(weave, params, primalVars, dualVars, isRoSy, ds);
        updatePrimalVars(weave, params, primalVars, dualVars, isRoSy);
    }
    
    delete ds;
}

static void identityMatrix(int n, Eigen::SparseMatrix<double> &I)
{
    std::vector<Eigen::Triplet<double> > coeffs;
    for (int i = 0; i < n; i++)
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

void LinearSolver::computeEnergy(const Weave &weave, SolverParams params, const Eigen::VectorXd &primalVars, const Eigen::VectorXd &dualVars, bool isRoSy )
{
    int nfaces = weave.fs->data().F.rows();
    int m = weave.fs->nFields();

    Eigen::SparseMatrix<double> D;
    Eigen::SparseMatrix<double> BTB;

    if(isRoSy)
        differentialOperator_rosy(weave, params, D);
    else
        differentialOperator(weave, params, D);
    massMatrix(weave, BTB);
    
    double term1 = 0.5 * dualVars.transpose() * BTB * dualVars;
    double term2 = 0.5 * params.lambdacompat * (primalVars + dualVars).transpose() * D.transpose() * D * (primalVars + dualVars);
    std::cout << " The current energy is " << term1 << " + " << term2
              << " = " 
              << term1 + term2 << std::endl;
}

void LinearSolver::updatePrimalVars(const Weave &weave, SolverParams params, Eigen::VectorXd &primalVars, Eigen::VectorXd &dualVars, bool isRoSy )
{
    int nfaces = weave.fs->data().F.rows();
    int m = weave.fs->nFields();
    std::cout << " pre-primal update ";
    computeEnergy(weave, params, primalVars, dualVars, isRoSy);
    
    for (int i = 0; i < nfaces * m; i++)
    {
        primalVars.segment<2>(2 * i) = primalVars.segment<2>(2 * i) + dualVars.segment<2>(2 * i);
        dualVars.segment<2>(2 * i) = primalVars.segment<2>(2 * i);
    }

    int nhandles = handles.size();

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
    computeEnergy(weave, params, primalVars, dualVars, isRoSy);

}

DualSolver *LinearSolver::buildDualUpdateSolver(const Weave &weave, SolverParams params, bool isRoSy)
{
    int nfaces = weave.fs->data().F.rows();
    int m = weave.fs->nFields();
    int nhandles = handles.size();
    int intedges = weave.fs->numInteriorEdges();

    Eigen::SparseMatrix<double> D;
    Eigen::SparseMatrix<double> H;
    Eigen::SparseMatrix<double> P;
    Eigen::SparseMatrix<double> curlOp;
    Eigen::SparseMatrix<double> BTB;
    
    if(isRoSy)
        differentialOperator_rosy(weave, params, D);
    else
        differentialOperator(weave, params, D);
    Eigen::VectorXd h0;
    Eigen::VectorXd dummy(2 * nfaces * m);
    dummy.setZero();
    handleConstraintOperator(weave, params, dummy, H, h0);
    int nhconstraints = h0.size();
    int matsize = 2 * nfaces * m + intedges * m + nhconstraints;
    curlOperator(weave, params, curlOp);

    if (params.disableCurlConstraint || isRoSy)
    {
        curlOp.setZero();
    }

    massMatrix(weave, BTB);
    double t = params.lambdacompat;
    Eigen::SparseMatrix<double> M = BTB + t * D.transpose() * D;
    
    std::vector<Eigen::Triplet<double> > dualCoeffs;

    for (int k=0; k<M.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(M,k); it; ++it)
        {
            dualCoeffs.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
        }
    }

    for (int k=0; k<curlOp.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(curlOp,k); it; ++it)
        {
            dualCoeffs.push_back(Eigen::Triplet<double>(it.row() + 2*m*nfaces, it.col(), it.value()));
        }
    }

    Eigen::SparseMatrix<double> CT = curlOp.transpose();
    for (int k=0; k<CT.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(CT,k); it; ++it)
        {
            dualCoeffs.push_back(Eigen::Triplet<double>(it.row(), it.col() + 2*m*nfaces, it.value()));
        }
    }

    for (int k=0; k<H.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(H,k); it; ++it)
        {
            dualCoeffs.push_back(Eigen::Triplet<double>(it.row() + 2*m*nfaces + intedges * m, it.col(), it.value()));
        }
    }

    Eigen::SparseMatrix<double> HT = H.transpose();
    for (int k=0; k<HT.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(HT,k); it; ++it)
        {
            dualCoeffs.push_back(Eigen::Triplet<double>(it.row(), it.col() + 2*m*nfaces + intedges * m, it.value()));
        }
    }

    std::cout << matsize << " matrix size " << std::endl;

    Eigen::SparseMatrix<double> dualMat;
    dualMat.resize(matsize, matsize);
    dualMat.setFromTriplets(dualCoeffs.begin(), dualCoeffs.end());
    std::cout << "Factoring matrix" << std::endl;
    DualSolver *ds = new DualSolver(dualMat);
    std::cout << "Done" << std::endl;
    return ds;
}

void LinearSolver::updateDualVars_new(const Weave &weave, SolverParams params, Eigen::VectorXd &primalVars, Eigen::VectorXd &dualVars, bool isRoSy, DualSolver *solver)
{
    // min_delta, \lambda   0.5 delta^2 + \lambda^T L (v + delta)
    // delta + L^T \lambda = 0
    // L delta + LL^T \lambda = 0
    // -L v + LL^T \lambda = 0
    // LL^T \lambda = Lv
    // delta = - L^T \lambda
    
    std::cout << "Building rhs" << std::endl;

    int nfaces = weave.fs->data().F.rows();
    int m = weave.fs->nFields();
    int nhandles = handles.size();
    int intedges = weave.fs->numInteriorEdges();

    Eigen::SparseMatrix<double> D;
    Eigen::SparseMatrix<double> H;
    Eigen::SparseMatrix<double> curlOp;
    
    if(isRoSy)
        differentialOperator_rosy(weave, params, D);
    else
        differentialOperator(weave, params, D);
    Eigen::VectorXd h0;
    handleConstraintOperator(weave, params, primalVars, H, h0);
    int nhconstraints = h0.size();
    int matsize = 2 * nfaces * m + intedges * m + nhconstraints;
    curlOperator(weave, params, curlOp);

    if (params.disableCurlConstraint || isRoSy)
    {
        curlOp.setZero();
    }

    double t = params.lambdacompat;


    Eigen::VectorXd rhs(matsize);
    rhs.setZero();
    rhs.segment(0, 2*nfaces*m) = -t * D.transpose() * D * primalVars;
    rhs.segment(2*nfaces*m, intedges * m ) = -curlOp * (primalVars);
    rhs.segment(2*nfaces*m + intedges * m, nhconstraints) = -h0;


    std::cout << "Solving" << std::endl;


    Eigen::VectorXd deltalambda;
    
    solver->solve(rhs, deltalambda);

    dualVars = deltalambda.segment(0, 2  * nfaces * m);
   
    std::cout << "  post-dual update ";
    computeEnergy(weave, params, primalVars, dualVars, isRoSy);
    std::cout << "Dual vars now " << dualVars.norm() << " geodesic residual " << (curlOp * (primalVars + dualVars)).norm() << std::endl;
}


// // TODO
// // *************************
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

      //          std::cout << permut << std::endl;

                int adj_field = -1;

                for (int field = 0; field < m; field++)
                {
                    vperm += permut(i, field) * weave.fs->v(g, field);
                    if (permut(i, field) != 0) 
                    {
                        adj_field = field;
                    }
                }
                Eigen::Matrix<double, 3, 2> B_f = weave.fs->data().Bs[f];
                Eigen::Matrix<double, 3, 2> B_g = weave.fs->data().Bs[g];

                if (params.edgeWeights(e) > 0.) 
                {
                    for (int j = 0; j < 2; j++)
                    {
                        coeffs.push_back(Eigen::Triplet<double>(term, weave.fs->vidx(f, i) + j, (B_f.transpose() * edge)[j]));
                        coeffs.push_back(Eigen::Triplet<double>(term, weave.fs->vidx(g, adj_field) + j, - permut(i, adj_field) * (B_g.transpose() * edge)[j]));
                    }
                }
                term++;
        }
    }

    curlOp.resize(intedges*m, 2 * m * nfaces);
    curlOp.setFromTriplets(coeffs.begin(), coeffs.end());
 //   curlOp.setZero();
}


// Either (# handles) x (2 * # faces * # fields) (for alignment handles only)
// or (2 * # handles) x (2 * # faces * # fields) (for hard handles)
void LinearSolver::handleConstraintOperator(const Weave &weave, SolverParams params, Eigen::VectorXd &primalVars, Eigen::SparseMatrix<double> &H, Eigen::VectorXd &h0)
{
    int m = weave.fs->nFields();
    int nfaces = weave.fs->data().F.rows();
    int nhandles = handles.size();
    std::vector<Eigen::Triplet<double> > coeffs;
    int term = 0;
    
    if(params.softHandleConstraint)
    {
        h0.resize(nhandles);
        h0.setZero();
    }
    else
    {
        h0.resize(2*nhandles);
        h0.setZero();
    }
    
    for (int i = 0; i < nhandles; i++)
    {
        int f = handles[i].face;
        if ( params.softHandleConstraint )
        {
            Eigen::Matrix2d Jf = weave.fs->data().Js.block<2, 2>(2 * f, 0);
            double area = weave.fs->faceArea(f);
            Eigen::Matrix2d BTB = weave.fs->data().Bs[f].transpose() * weave.fs->data().Bs[f];

            Eigen::Vector2d dir = (Jf * handles[i].dir).transpose() * BTB;
            int idx = 2 * (handles[i].face * m + handles[i].field);
            coeffs.push_back(Eigen::Triplet<double>(term, idx, dir(0)));
            coeffs.push_back(Eigen::Triplet<double>(term, idx + 1, dir(1)));
            term++;
        }
        else
        {
            int idx = 2 * (handles[i].face * m + handles[i].field);
            coeffs.push_back(Eigen::Triplet<double>(term, idx, 1.0));
            coeffs.push_back(Eigen::Triplet<double>(term+1, idx+1, 1.0));
            Eigen::Vector2d dir = handles[i].dir;
            Eigen::Vector3d extdir = weave.fs->data().Bs[f]*dir;
            dir /= extdir.norm();
            h0[term] = -dir[0];
            h0[term+1] = -dir[1];
            term += 2;
        }
    }
    if(params.softHandleConstraint)
    {
        H.resize(nhandles, 2 * m * nfaces);
        H.setFromTriplets(coeffs.begin(), coeffs.end());
        h0 = H * primalVars;
    }
    else
    {
        H.resize(2*nhandles, 2 * m * nfaces);
        H.setFromTriplets(coeffs.begin(), coeffs.end());
    }
}

// TODO: Steal this from prev version.
// *************************

void LinearSolver::differentialOperator_rosy(const Weave &weave, SolverParams params, Eigen::SparseMatrix<double> &D)
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
        for (int side = 0; side < 2; side++)
        {
            int f = (side == 0 ? weave.fs->data().E(e, 0) : weave.fs->data().E(e, 1));
            int g = (side == 0 ? weave.fs->data().E(e, 1) : weave.fs->data().E(e, 0));

            Eigen::Matrix2d Tgf = weave.fs->data().Ts.block<2, 2>(2 * e, 2 - 2 * side);
            Eigen::Matrix2d Tgf_rosy = weave.fs->data().Ts_rosy.block<2, 2>(2 * e, 2 - 2 * side);
            Eigen::Matrix2d Tgf_rosy_inv = Tgf_rosy.inverse();

            Eigen::Matrix2d Tgf_rosy_power = Tgf_rosy_inv;
            for (int s = 0; s < params.rosyN - 2; s++)
            {
                Tgf_rosy_power *= Tgf_rosy_inv;
            }

            Eigen::Matrix2d id;
            id.setIdentity();
            Eigen::Matrix<double, 3, 2> id_ambient = weave.fs->data().Bs[f] * id;
            Eigen::Matrix<double, 3, 2> transported = weave.fs->data().Bs[f] * Tgf_rosy_power * Tgf;
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    coeffs.push_back(Eigen::Triplet<double>(term, weave.fs->vidx(f, 0) + j, id_ambient(i,j)));
                    coeffs.push_back(Eigen::Triplet<double>(term, weave.fs->vidx(g, 0) + j, -transported(i,j) ));
                }
                term++;
            }
        }
    }


    int nfaces = weave.fs->data().F.rows();
    D.resize(6 * intedges, 2*nfaces);
    D.setFromTriplets(coeffs.begin(), coeffs.end());
}

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
            for (int side = 0; side < 2; side++)
            {
                int f = (side == 0 ? weave.fs->data().E(e, 0) : weave.fs->data().E(e, 1));
                int g = (side == 0 ? weave.fs->data().E(e, 1) : weave.fs->data().E(e, 0));
            //    Eigen::Matrix2d Jf = weave.fs->data().Js.block<2, 2>(2 * f, 0);
            //    Eigen::Vector2d vif = weave.fs->v(f, i);
            //    Eigen::Vector2d cdiff = weave.fs->data().cDiffs.row(2 * e + side);
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
                    innervec[coeff] = sqrt(params.edgeWeights(e));
                    Eigen::Vector2d dE = innervec;
                    for (int k = 0; k < 2; k++)
                    {
                        coeffs.push_back(Eigen::Triplet<double>(term, weave.fs->vidx(f, i) + k, dE[k]));
                    }
                    for (int field = 0; field < m; field++)
                    {
                        dE = -permut(i, field) * Tgf.transpose() * innervec;
                        for (int k = 0; k < 2; k++)
                        {
                            coeffs.push_back(Eigen::Triplet<double>(term, weave.fs->vidx(g, field) + k, dE[k]));
                        }
                    }

                    term++;
                }
            }
        }
    }


    int nfaces = weave.fs->data().F.rows();
    D.resize(4 * intedges * m, 2*nfaces*m);
    D.setFromTriplets(coeffs.begin(), coeffs.end());
}

