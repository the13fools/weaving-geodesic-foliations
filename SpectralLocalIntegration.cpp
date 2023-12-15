#include "SpectralLocalIntegration.h"
#include <igl/cotmatrix_entries.h>
#include <vector>
#include <Eigen/Sparse>

void SpectralLocalIntegration::locallyIntegrateOneComponent(const Surface &surf, const Eigen::MatrixXd &v, Eigen::VectorXd &s)
{
    // build edge metric matrix and inverse (cotan weights)
    Eigen::MatrixXd C;
    igl::cotmatrix_entries(surf.data().V, surf.data().F, C);
    int nedges = surf.nEdges();
    Eigen::VectorXd edgeMetric(nedges);
    edgeMetric.setZero();
    int nfaces = surf.nFaces();
    for(int i=0; i<nfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int eidx = surf.data().faceEdges(i, j);
            edgeMetric[eidx] += C(i,j);
        }                
    } 

    // face mass matrix
    Eigen::VectorXd faceAreas;
    igl::doublearea(surf.data().V, surf.data().F, faceAreas);
    faceAreas *= 0.5;

    std::cout << "Built mass matrices" << std::endl;

    // face Laplacian        
    std::vector<Eigen::Triplet<double> > DfaceCoeffs;
    for (int i = 0; i < nedges; i++)
    {
        int f0 = surf.data().E(i,0);
        int f1 = surf.data().E(i,1);
        if (f0 != -1 && f1 != -1)
        {
            DfaceCoeffs.push_back(Eigen::Triplet<double>(i, f0, -1.0));
            DfaceCoeffs.push_back(Eigen::Triplet<double>(i, f1, 1.0));
        }
    }
    Eigen::SparseMatrix<double> Dface(nedges, nfaces);
    Dface.setFromTriplets(DfaceCoeffs.begin(), DfaceCoeffs.end());

    std::vector<Eigen::Triplet<double> > inverseEdgeMetricCoeffs;
    for (int i = 0; i < nedges; i++)
    {
        double denom = std::max(1e-6, edgeMetric[i]);        
        inverseEdgeMetricCoeffs.push_back(Eigen::Triplet<double>(i, i, 1.0 / denom));
    }
    Eigen::SparseMatrix<double> inverseEdgeMetric(nedges, nedges);
    inverseEdgeMetric.setFromTriplets(inverseEdgeMetricCoeffs.begin(), inverseEdgeMetricCoeffs.end());

    Eigen::SparseMatrix<double> Lface = Dface.transpose() * inverseEdgeMetric * Dface;


    // A matrix

    std::vector<Eigen::Triplet<double> > ACoeffs;

    for(int i=0; i<nfaces; i++)
    {
        for(int j=0; j<3; j++)
        {
            ACoeffs.push_back(Eigen::Triplet<double>(3*i+j, 3*i+j, faceAreas[i]));
        }            
    }

    for(int i=0; i<Lface.outerSize(); i++)
    {
        for(Eigen::SparseMatrix<double>::InnerIterator it(Lface, i); it; ++it)
        {
            ACoeffs.push_back(Eigen::Triplet<double>(3*nfaces+it.row(), 3*nfaces+it.col(), sreg_*it.value()));
        }
    }
    Eigen::SparseMatrix<double> A(4*nfaces, 4*nfaces);
    A.setFromTriplets(ACoeffs.begin(), ACoeffs.end());

    // B matrix
    std::vector<Eigen::Triplet<double> > Bcoeffs;
    for(int i=0; i<nfaces; i++)
    {
        for(int j=0; j<3; j++)
        {
            Bcoeffs.push_back(Eigen::Triplet<double>(3*i+j, 3*i+j, faceAreas[i]));
        }
        Bcoeffs.push_back(Eigen::Triplet<double>(3*nfaces+i, 3*nfaces+i, faceAreas[i]));
    }
    Eigen::SparseMatrix<double> B(4*nfaces, 4*nfaces);
    B.setFromTriplets(Bcoeffs.begin(), Bcoeffs.end());

    // B matrix
    std::vector<Eigen::Triplet<double> > BInvcoeffs;
    for(int i=0; i<nfaces; i++)
    {
        for(int j=0; j<3; j++)
        {
            BInvcoeffs.push_back(Eigen::Triplet<double>(3*i+j, 3*i+j, 1.0/faceAreas[i]));
        }
        BInvcoeffs.push_back(Eigen::Triplet<double>(3*nfaces+i, 3*nfaces+i, 1.0/faceAreas[i]));
    }
    Eigen::SparseMatrix<double> BInv(4*nfaces, 4*nfaces);
    BInv.setFromTriplets(BInvcoeffs.begin(), BInvcoeffs.end());

    // constraint matrix
    std::vector<Eigen::Triplet<double> > DCoeffs;

    int nconstraints = 0;
    for(int i=0; i<nedges; i++)
    {
        int f0 = surf.data().E(i,0);
        int f1 = surf.data().E(i,1);
        if( f0 != -1 && f1 != -1)
            nconstraints++;
    }

    int idx=0;
    for(int i=0; i<nedges; i++)
    {
        int f0 = surf.data().E(i,0);
        int f1 = surf.data().E(i,1);
        int vert0 = surf.data().edgeVerts(i, 0);
        int vert1 = surf.data().edgeVerts(i, 1);
        Eigen::Vector3d v0 = surf.data().V.row(vert0).transpose();
        Eigen::Vector3d v1 = surf.data().V.row(vert1).transpose();
        Eigen::Vector3d edgeVec = v1-v0;

        if (f0 == -1 || f1 == -1)
            continue;

        Eigen::Vector2d field0 = v.row(f0).transpose();
        Eigen::Vector2d field1 = v.row(f1).transpose();
        Eigen::Vector3d scaledvec0 = surf.data().Bs[f0] * surf.data().Js.block<2, 2>(2 * f0, 0) * field0;
        Eigen::Vector3d scaledvec1 = surf.data().Bs[f1] * surf.data().Js.block<2, 2>(2 * f1, 0) * field1;
        // w part
        for(int j=0; j<3; j++)
        {
            DCoeffs.push_back(Eigen::Triplet<double>(idx, 3*f0+j, -edgeVec(j)));
            DCoeffs.push_back(Eigen::Triplet<double>(idx, 3*f1+j, edgeVec(j)));
        }
        // s part
        DCoeffs.push_back(Eigen::Triplet<double>(idx, 3*nfaces + f0, -scaledvec0.dot(edgeVec)));
        DCoeffs.push_back(Eigen::Triplet<double>(idx, 3*nfaces + f1, scaledvec1.dot(edgeVec)));
        idx++;
    }
    Eigen::SparseMatrix<double> D(nconstraints, 4*nfaces);
    D.setFromTriplets(DCoeffs.begin(), DCoeffs.end());

    Eigen::SparseMatrix<double> Areg = A + 1e-6 * B;

    std::cout << "Factoring A" << std::endl;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solveA(Areg);
    if(solveA.info() != Eigen::Success)
        std::cout << "failed" << std::endl;

    std::vector<Eigen::Triplet<double> > DDTregCoeffs;
    for(int i=0; i<nconstraints; i++)
        DDTregCoeffs.push_back(Eigen::Triplet<double>(i, i, 1e-6));
    Eigen::SparseMatrix<double> DDTreg(nconstraints, nconstraints);
    DDTreg.setFromTriplets(DDTregCoeffs.begin(), DDTregCoeffs.end());
    Eigen::SparseMatrix<double> DT = D.transpose();
    Eigen::SparseMatrix<double> DDT = DDTreg + D * BInv * DT;
    std::cout << "Factoring DDT" << std::endl;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solveC(DDT);
    if(solveC.info() != Eigen::Success)
        std::cout << "failed" << std::endl;

    std::cout << "Starting inverse power interation" << std::endl;
    Eigen::VectorXd x(4*nfaces);
    srand(0);
    x.setRandom();
    x /= sqrt(x.transpose() * (B * x));
    for(int i=0; i<100; i++)
    {
        Eigen::VectorXd newx = solveA.solve(B*x);
        Eigen::VectorXd rhs = D*newx;
        Eigen::VectorXd lambda = solveC.solve(rhs);
        Eigen::VectorXd projx = newx - BInv * (D.transpose() * lambda);
        double xnorm = sqrt(projx.transpose() * (B*projx));
        x = projx/xnorm;            
    }

    std::cout << "Rayleigh quotient: " << x.transpose() * (A * x) / (x.transpose() * (B*x)) << std::endl;

    s.resize(nfaces);
    for (int i = 0; i < nfaces; i++)
        s[i] = x[3*nfaces + i];

    // fix s sign
    double totals = s.sum();
    if (totals < 0)
        s *= -1.0;
}
