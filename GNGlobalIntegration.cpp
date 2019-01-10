#include "GNGlobalIntegration.h"
#include <iostream>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Geometry>

static const double PI = 3.1415926535898;

void GNGlobalIntegration::globallyIntegrateOneComponent(const Surface &surf, const Eigen::MatrixXd &v, Eigen::VectorXd &scales, Eigen::VectorXd &theta)
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
    
    assert(scales.size() == nfaces);

    int totalIter = outerIters_;
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
        srand(0);
        eigenVec.setRandom();
        eigenVec /= eigenVec.norm();
        for (int i = 0; i < powerIters_; i++)
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
            if (curPred > PI) curPred -= 2 * PI;
            if (curPred < -PI) curPred += 2 * PI;
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

