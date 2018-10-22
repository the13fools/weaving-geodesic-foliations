#include "LinearSolver.h"
#include <map>
#include <vector>
#include <set>
#include <cassert>
#include <Eigen/Geometry>
#include <iostream>

// LinearSolver::LinearSolver(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) : V(V), F(F)
// {
//     std::map<std::pair<int, int>, std::vector<int> > edges;
//     for (int i = 0; i < F.rows(); i++)
//     {
//         for (int j = 0; j < 3; j++)
//         {
//             int from = F(i, j);
//             int to = F(i, (j + 1) % 3);
//             int idx = 0;
//             if (from > to)
//             {
//                 idx = 1;
//                 std::swap(from, to);
//             }
//             std::pair<int, int> e(from, to);
//             auto it = edges.find(e);
//             if (it == edges.end())
//             {
//                 edges[e].resize(2);
//                 edges[e][idx] = i;
//                 edges[e][1 - idx] = -1;
//             }
//             else
//             {
//                 edges[e][idx] = i;
//             }
//         }
//     }

//     int intedges = 0;
//     for(auto &it : edges)
//     { 
//         if (it.second[0] != -1 && it.second[1] != -1)
//             intedges++;
//     }

//     E.resize(intedges, 4);
//     int row = 0;
//     for (auto &it : edges)
//     {
//         if (it.second[0] != -1 && it.second[1] != -1)
//         {
//             E(row, 0) = it.first.first;
//             E(row, 1) = it.first.second;
//             E(row, 2) = it.second[0];
//             E(row, 3) = it.second[1];
//             row++;
//         }
//     }

//     assert(row == intedges);

//     int nfaces = F.rows();
//     centroids.resize(nfaces,2);
//     for (int i = 0; i < nfaces; i++)
//     {
//         Eigen::Vector2d centroid(0,0);
//         for (int j = 0; j < 3; j++)
//         {
//             centroid[0] += V(F(i, j), 0);
//             centroid[1] += V(F(i, j), 1);
//         }
//         centroid /= 3.0;
//         centroids.row(i) = centroid;
//     }
// }

void LinearSolver::addHandle(const Handle &h)
{
    handles.push_back(h);
}

int LinearSolver::numPrimalDOFs()
{
    return 2 * F.rows();
}

int LinearSolver::numDualDOFs()
{
    return 2 * F.rows();
}

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
void LinearSolver::setFaceEnergies(const Weave &weave, const Eigen::VectorXd &primalVars, const Eigen::VectorXd &dualVars, Eigen::VectorXd &faceEnergies)
{
    faceEnergies.resize(weave.fs->data().F.rows());
    faceEnergies.setZero();

    int nedges = weave.fs->data().E.rows();
    for (int i = 0; i < nedges; i++)
    {
        int face1 = weave.fs->data().E(i, 2);
        int face2 = weave.fs->data().E(i, 3);
        Eigen::Vector2d edgevec;
        edgevec[0] = V(E(i, 1), 0) - V(E(i, 0), 0);
        edgevec[1] = V(E(i, 1), 1) - V(E(i, 0), 1);
        edgevec /= edgevec.norm();
        
        double energy = (primalVars.segment<2>(2 * face2).dot(edgevec) - primalVars.segment<2>(2 * face1).dot(edgevec));
        faceEnergies[face1] += energy * energy;
        faceEnergies[face2] += energy * energy;
    }
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

void LinearSolver::updatePrimalVars(const Weave &weave, Eigen::VectorXd &primalVars, const Eigen::VectorXd &dualVars, double smoothingCoeff)
{
    int nfaces = weave.fs->data().F.rows();
    for (int i = 0; i < nfaces; i++)
    {
        primalVars.segment<2>(2 * i) = primalVars.segment<2>(2 * i) + dualVars.segment<2>(2 * i);
    }

    Eigen::SparseMatrix<double> D;
    Eigen::SparseMatrix<double> P;
    Eigen::SparseMatrix<double> I;
    differentialOperator(D);
    unconstrainedProjection(P);
    int nhandles = handles.size();
    identityMatrix(2 * nfaces, I);
    double t = smoothingCoeff;
    Eigen::SparseMatrix<double> op = I + t * D.transpose() * D;
    Eigen::SparseMatrix<double> M = P.transpose() * op * P;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver(M);
    Eigen::VectorXd handlevals(2 * nfaces);
    handlevals.setZero();
    for (int i = 0; i < nhandles; i++)
    {
        handlevals.segment<2>(2 * handles[i].face) = handles[i].dir / handles[i].dir.norm() * 100;
    }
    Eigen::VectorXd rhs = P.transpose() * (primalVars - op * handlevals);
    Eigen::VectorXd newprim = solver.solve(rhs);
    primalVars = P * newprim + handlevals;
    for (int i = 0; i < nfaces; i++)
    {
        primalVars.segment<2>(2 * i).normalize();
    }

    for (int i = 0; i < nhandles; i++)
    {
        primalVars.segment<2>(2 * handles[i].face) = handles[i].dir / handles[i].dir.norm();
    }
}


// TODO: Steal this from prev version
// *************************
void LinearSolver::curlOperator(Eigen::SparseMatrix<double> &curlOp)
{
    std::vector<Eigen::Triplet<double> > coeffs;
    int nedges = E.rows();
    for(int i=0; i<nedges; i++)
    {
        int face1 = E(i, 2);
        int face2 = E(i, 3);
        Eigen::Vector2d edgevec;
        edgevec[0] = V(E(i, 1), 0) - V(E(i, 0), 0);
        edgevec[1] = V(E(i, 1), 1) - V(E(i, 0), 1);
        edgevec /= edgevec.norm();
        for (int j = 0; j < 2; j++)
        {
            coeffs.push_back(Eigen::Triplet<double>(i, 2 * face2 + j, edgevec[j]));
            coeffs.push_back(Eigen::Triplet<double>(i, 2 * face1 + j, -edgevec[j]));
        }
    }
    int nfaces = F.rows();
    curlOp.resize(nedges, 2 * nfaces);
    curlOp.setFromTriplets(coeffs.begin(), coeffs.end());
}

// TODO: Steal this from prev version.
// *************************
void LinearSolver::differentialOperator(Eigen::SparseMatrix<double> &D)
{
    std::vector<Eigen::Triplet<double> > coeffs;
    int nedges = E.rows();
    for (int i = 0; i < nedges; i++)
    {
        int face1 = E(i, 2);
        int face2 = E(i, 3);
        for (int j = 0; j < 2; j++)
        {
            coeffs.push_back(Eigen::Triplet<double>(2 * i + j, 2 * face2 + j, 1.0));
            coeffs.push_back(Eigen::Triplet<double>(2 * i + j, 2 * face1 + j, -1.0));
        }
    }
    int nfaces = F.rows();
    D.resize(2 * nedges, 2 * nfaces);
    D.setFromTriplets(coeffs.begin(), coeffs.end());
}
void LinearSolver::updateDualVars(const Eigen::VectorXd &primalVars, Eigen::VectorXd &dualVars)
{
    // min_delta, \lambda   0.5 delta^2 + \lambda^T L (v + delta)
    // delta + L^T \lambda = 0
    // L delta + LL^T \lambda = 0
    // -L v + LL^T \lambda = 0
    // LL^T \lambda = Lv
    // delta = - L^T \lambda
    Eigen::SparseMatrix<double> curlOp;
    curlOperator(curlOp);
    int nfaces = F.rows();
    Eigen::SparseMatrix<double> V;
    unconstrainedProjection(V);
    int nedges = E.rows();

    std::vector<Eigen::Triplet<double> > regcoeffs;
    double reg = 1e-8;
    for (int i = 0; i < nedges; i++)
    {
        regcoeffs.push_back(Eigen::Triplet<double>(i, i, reg));
    }
    Eigen::SparseMatrix<double> Reg(nedges, nedges);
    Reg.setFromTriplets(regcoeffs.begin(), regcoeffs.end());

    Eigen::SparseMatrix<double> M = curlOp * V * V.transpose() * curlOp.transpose() + Reg;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver(M);
    
    Eigen::VectorXd rhs(nedges);
    rhs = curlOp * primalVars;
    
    Eigen::VectorXd lambda = solver.solve(rhs);

    dualVars = -V * V.transpose() * curlOp.transpose() * lambda;

    std::cout << "Dual vars now " << dualVars.norm() << " geodesic residual " << (curlOp * (primalVars + dualVars)).norm() << std::endl;
}

void LinearSolver::unconstrainedProjection(Eigen::SparseMatrix<double> &proj)
{
    int nfaces = F.rows();
    int nhandles = handles.size();
    std::set<int> handlefaces;
    for (int i = 0; i < nhandles; i++)
        handlefaces.insert(handles[i].face);
    std::vector<Eigen::Triplet<double> > Vcoeffs;
    int col = 0;
    for (int i = 0; i < nfaces; i++)
    {
        if (handlefaces.count(i))
            continue;
        Vcoeffs.push_back(Eigen::Triplet<double>(2 * i, col, 1.0));
        Vcoeffs.push_back(Eigen::Triplet<double>(2 * i + 1, col + 1, 1.0));
        col += 2;
    }
    assert(col == 2*nfaces - 2*nhandles);
    proj.resize(2 * nfaces, 2*nfaces - 2*nhandles);
    proj.setFromTriplets(Vcoeffs.begin(), Vcoeffs.end());

}