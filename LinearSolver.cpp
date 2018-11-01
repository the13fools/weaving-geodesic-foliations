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


void LinearSolver::addHandle(const Handle &h)
{
    handles.push_back(h);
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

void LinearSolver::updatePrimalVars(const Weave &weave, Eigen::VectorXd &primalVars, const Eigen::VectorXd &dualVars, double smoothingCoeff)
{
    int nfaces = weave.fs->data().F.rows();
    int m = weave.fs->nFields();
    
    for (int i = 0; i < nfaces * m; i++)
    {
        primalVars.segment<2>(2 * i) = primalVars.segment<2>(2 * i) + dualVars.segment<2>(2 * i);
    }

    Eigen::SparseMatrix<double> D;
    Eigen::SparseMatrix<double> P;
    Eigen::SparseMatrix<double> I;
    differentialOperator(weave, D);
    unconstrainedProjection(weave, P);
    int nhandles = handles.size();
    identityMatrix(2 * nfaces * m, I);
    double t = smoothingCoeff;
    Eigen::SparseMatrix<double> op = I + t * D.transpose() * D;
    Eigen::SparseMatrix<double> M = P.transpose() * op * P;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver(M);
    Eigen::VectorXd handlevals(2 * nfaces * m);
    handlevals.setZero();
    for (int i = 0; i < nhandles; i++)
    {
        handlevals.segment<2>(2 * handles[i].face) = handles[i].dir / handles[i].dir.norm() * 100;
        std::cout << " here's a handle " << std::endl;
    }
    Eigen::VectorXd rhs = P.transpose() * (primalVars - op * handlevals);

    Eigen::VectorXd newprim = solver.solve(rhs);
    primalVars = P * newprim + handlevals;
    for (int i = 0; i < nfaces; i++)
    {
        for ( int cover = 0; cover < m; cover++)
        {
            Eigen::Matrix<double, 3, 2> B_f = weave.fs->data().Bs[i];
            Eigen::Vector3d ambient = (B_f * primalVars.segment<2>(2 * i * m + 2 * cover));
            primalVars.segment<2>(2 * i * m + 2 * cover) /= ambient.norm();
        }
    }

    for (int i = 0; i < nhandles; i++)
    {
        primalVars.segment<2>(2 * handles[i].face) = handles[i].dir / handles[i].dir.norm();
    }

}



// TODO
// *************************
void LinearSolver::curlOperator(const Weave &weave, Eigen::SparseMatrix<double> &curlOp)
{
    std::vector<Eigen::Triplet<double> > coeffs;
    int nedges = weave.fs->data().E.rows();
    int intedges = weave.fs->numInteriorEdges();
    int m = weave.fs->nFields();

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

                for (int j = 0; j < 2; j++)
                {
                    coeffs.push_back(Eigen::Triplet<double>(term, weave.fs->vidx(f, i) + j, (B_f.transpose() * edge)[j]));
                    coeffs.push_back(Eigen::Triplet<double>(term, weave.fs->vidx(g, adj_field) + j, -(B_g.transpose() * edge)[j]));
                }
                term++;
        }
    }

    int nfaces = weave.fs->data().F.rows();
    curlOp.resize(intedges * m, 2 * m * nfaces);
    curlOp.setFromTriplets(coeffs.begin(), coeffs.end());
}

// TODO: Steal this from prev version.
// *************************
void LinearSolver::differentialOperator(const Weave &weave, Eigen::SparseMatrix<double> &D)
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
                    for (int k = 0; k < 2; k++)
                    {
                        coeffs.push_back(Eigen::Triplet<double>(term, weave.fs->vidx(f, i) + k, B_f(coeff, k)));
                        coeffs.push_back(Eigen::Triplet<double>(term, weave.fs->vidx(g, adj_field) + k, -B_g(coeff, k)));
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
    D.resize(3 * intedges * m, 2  * nfaces * m);
    D.setFromTriplets(coeffs.begin(), coeffs.end());
}

void LinearSolver::updateDualVars(const Weave &weave, const Eigen::VectorXd &primalVars, Eigen::VectorXd &dualVars)
{
    // min_delta, \lambda   0.5 delta^2 + \lambda^T L (v + delta)
    // delta + L^T \lambda = 0
    // L delta + LL^T \lambda = 0
    // -L v + LL^T \lambda = 0
    // LL^T \lambda = Lv
    // delta = - L^T \lambda
    Eigen::SparseMatrix<double> curlOp;
    curlOperator(weave, curlOp);
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
        handlefaces.insert(handles[i].face);
    std::vector<Eigen::Triplet<double> > Vcoeffs;
    int col = 0;
    for (int i = 0; i < nfaces * m; i++)
    {
        if (handlefaces.count(i / m))
            continue;
        Vcoeffs.push_back(Eigen::Triplet<double>(2 * i, col, 1.0));
        Vcoeffs.push_back(Eigen::Triplet<double>(2 * i + 1, col + 1, 1.0));
        col += 2;
    }
    assert(col == 2*nfaces*m - 2*nhandles);
    proj.resize(2 * nfaces * m, 2*nfaces * m - 2*nhandles);
    proj.setFromTriplets(Vcoeffs.begin(), Vcoeffs.end());

}