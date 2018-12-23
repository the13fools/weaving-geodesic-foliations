#include "OurFieldIntegration.h"
#include <vector>
#include <igl/cotmatrix_entries.h>
#include "Surface.h"
#include <Eigen/Sparse>

void OurFieldIntegration::initializeS(const Surface &surf, const Eigen::MatrixXd &v, Eigen::VectorXd &s)
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
        inverseEdgeMetricCoeffs.push_back(Eigen::Triplet<double>(i, i, 1.0 / edgeMetric[i]));
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
    Eigen::SparseMatrix<double> DDT = DDTreg + D * BInv * D.transpose();
    std::cout << "Factoring DDT" << std::endl;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solveC(DDT);
    if(solveC.info() != Eigen::Success)
        std::cout << "failed" << std::endl;

    std::cout << "Starting inverse power interation" << std::endl;
    Eigen::VectorXd x(4*nfaces);
    x.setRandom();
    x /= sqrt(x.transpose() * (B * x));
    for(int i=0; i<1000; i++)
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

    double maxS = 0;
    for(int i=0; i<nfaces; i++)
    {
        if ( fabs(s[i]) > maxS ) 
        {
            maxS = fabs(s[i]);
        }
    }

    double s_scale = 3.1415 / surf.data().averageEdgeLength / maxS;
    s *= s_scale;

    // fix s sign
    double totals = s.sum();
    if (totals < 0)
        s *= -1.0;
}

void OurFieldIntegration::globalSThetaSolve(const Surface &surf, const Eigen::MatrixXd &v, Eigen::VectorXd &s, Eigen::VectorXd &theta)
{
    theta.setZero();

    int nfaces = surf.nFaces();
    int nverts = surf.nVerts();

    std::cout << "nfaces: " << nfaces << std::endl;
    std::cout << "nverts: " << nverts << std::endl;
    std::vector<int> rowsL;
    std::vector<int> colsL;
    std::vector<double> difVecUnscaled;
    for (int fId = 0; fId < nfaces; fId++)
    { // Compute rowsL, colsL, difVecUnscaled
        int vId0 = surf.data().F(fId, 0);
        int vId1 = surf.data().F(fId, 1);
        int vId2 = surf.data().F(fId, 2);
        rowsL.push_back(vId0); rowsL.push_back(vId1); rowsL.push_back(vId2);
        colsL.push_back(vId1); colsL.push_back(vId2); colsL.push_back(vId0);
        Eigen::Vector3d p0 = surf.data().V.row(vId0);
        Eigen::Vector3d p1 = surf.data().V.row(vId1);
        Eigen::Vector3d p2 = surf.data().V.row(vId2);
        Eigen::Vector3d e01 = p0 - p1;
        Eigen::Vector3d e12 = p1 - p2;
        Eigen::Vector3d e20 = p2 - p0;
        Eigen::Vector3d faceVec;
        Eigen::Vector2d v0 = v.row(fId).transpose();
        faceVec = surf.data().Bs[fId] * v0; // The original vec            

        faceVec = faceVec.cross(surf.faceNormal(fId));
        faceVec /= faceVec.norm();

        difVecUnscaled.push_back(e01.dot(faceVec));
        difVecUnscaled.push_back(e12.dot(faceVec));
        difVecUnscaled.push_back(e20.dot(faceVec));
    }
    assert((rowsL.size() == 3 * nfaces) && (colsL.size() == 3 * nfaces) && (difVecUnscaled.size() == 3 * nfaces));

    Eigen::VectorXd scales = globalScale_ * s;
    assert(scales.size() == nfaces);
    
    int totalIter = 6;
    for (int iter = 0; iter < totalIter; iter++)
    {
        std::vector<double> difVec;
        for (int i = 0; i < difVecUnscaled.size(); i++)
            difVec.push_back(difVecUnscaled[i] * scales(i / 3));
        std::vector<Eigen::Triplet<double> > sparseContent;
        for (int i = 0; i < rowsL.size(); i++)
            sparseContent.push_back(Eigen::Triplet<double>(rowsL[i], colsL[i], 1));
        Eigen::SparseMatrix<double> TP(nverts, nverts);
        Eigen::SparseMatrix<double> TPTran(nverts, nverts);
        TP.setFromTriplets(sparseContent.begin(), sparseContent.end());
        TPTran = TP.transpose();
        TP += TPTran;
        std::vector<int> degree;
        for (int i = 0; i < nverts; i++)
            degree.push_back(TP.row(i).sum());
        std::vector<Eigen::Triplet<double> > AContent;
        for (int i = 0; i < rowsL.size(); i++)
        {
            double cVal = cos(difVec[i]);
            double sVal = sin(difVec[i]);
            AContent.push_back(Eigen::Triplet<double>(2 * rowsL[i], 2 * colsL[i], cVal));
            AContent.push_back(Eigen::Triplet<double>(2 * rowsL[i], 2 * colsL[i] + 1, -sVal));
            AContent.push_back(Eigen::Triplet<double>(2 * rowsL[i] + 1, 2 * colsL[i], sVal));
            AContent.push_back(Eigen::Triplet<double>(2 * rowsL[i] + 1, 2 * colsL[i] + 1, cVal));
        }
        Eigen::SparseMatrix<double> Amat(2 * nverts, 2 * nverts);
        Eigen::SparseMatrix<double> Amat_tran(2 * nverts, 2 * nverts);
        Amat.setFromTriplets(AContent.begin(), AContent.end());
        Amat_tran = Amat.transpose();
        Amat += Amat_tran;
        //
        std::vector<Eigen::Triplet<double> > LContent;
        for (int i = 0; i < 2 * nverts; i++)
            LContent.push_back(Eigen::Triplet<double>(i, i, degree[int(i / 2)]));
        Eigen::SparseMatrix<double> Lmat(2 * nverts, 2 * nverts);
        Lmat.setFromTriplets(LContent.begin(), LContent.end());
        Lmat -= Amat;
        // Eigen Decompose
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solverL(Lmat);
        Eigen::VectorXd eigenVec(Lmat.rows());
        eigenVec.setRandom();
        eigenVec /= eigenVec.norm();
        for (int i = 0; i < 10; i++)
        {
            eigenVec = solverL.solve(eigenVec);
            eigenVec /= eigenVec.norm();
        }
        double eigenVal = eigenVec.transpose() * (Lmat * eigenVec);
        std::cout << "Current iteration = " << iter << " currents error is: " << eigenVal << std::endl;
        // Extract the function value
        theta.resize(nverts);
        for (int i = 0; i < nverts; i++)
        {
            double curCos = eigenVec(2 * i);
            double curSin = eigenVec(2 * i + 1);
            double normalizer = sqrt(curCos * curCos + curSin * curSin);
            double curFunc = acos(curCos / normalizer);
            if (curSin < 0)
                curFunc = -curFunc;
            theta[i] = curFunc;
        }
        ////
        //// Re-compute face scales
        std::vector<double> difVecPred;
        for (int i = 0; i < rowsL.size(); i++)
        {
            double curPred = theta[rowsL[i]] - theta[colsL[i]];
            if (curPred > M_PI) curPred -= 2 * M_PI;
            if (curPred < -M_PI) curPred += 2 * M_PI;
            difVecPred.push_back(curPred);
        }
        Eigen::VectorXd bScales(nfaces);
        std::vector<double> diagAScales;
        // TODO: AScalesMat is constant
        for (int i = 0; i < rowsL.size(); i = i + 3)
        {
            double bVal = 0;
            double diagAVal = 0;
            for (int j = 0; j < 3; j++)
            {
                bVal += difVecPred[i + j] * difVecUnscaled[i + j];
                diagAVal += difVecUnscaled[i + j] * difVecUnscaled[i + j];
            }
            bScales(i / 3) = bVal;
            diagAScales.push_back(diagAVal);
        }
        // Construct A
        // TODO mu and lambda
        std::vector<Eigen::Triplet<double> > AScalesContent;
        for (int i = 0; i < nfaces; i++)
            AScalesContent.push_back(Eigen::Triplet<double>(i, i, diagAScales[i]));
        Eigen::SparseMatrix<double> AScalesMat(nfaces, nfaces);
        AScalesMat.setFromTriplets(AScalesContent.begin(), AScalesContent.end());
        // Solve for scale
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solverScales(AScalesMat);
        Eigen::VectorXd curScales = solverScales.solve(bScales);

        for (int i = 0; i < nfaces; i++)
            scales(i) = curScales(i);
    }    
}

void OurFieldIntegration::integrateOneComponent(const Surface &surf, const Eigen::MatrixXd &v, Eigen::VectorXd &theta)
{
    Eigen::VectorXd s;
    initializeS(surf, v, s);
    theta.resize(surf.nVerts());
    theta.setZero();
    for (int i = 0; i < surf.nFaces(); i++)
    {
        for(int j=0; j<3; j++)
            theta[surf.data().F(i, j)] += s[i];
    }
    globalSThetaSolve(surf, v, s, theta);
}