#include "CurlLocalIntegration.h"
#include <igl/cotmatrix_entries.h>
#include <vector>
#include <Eigen/Sparse>

void CurlLocalIntegration::locallyIntegrateOneComponent(const Surface &surf, const Eigen::MatrixXd &v, Eigen::VectorXd &s)
{
    // build edge metric matrix and inverse (cotan weights)
    Eigen::MatrixXd C;
    igl::cotmatrix_entries(surf.data().V, surf.data().F, C);
    int nedges = surf.nEdges();
    Eigen::VectorXd edgeMetric(nedges);
    edgeMetric.setZero();
    int nfaces = surf.nFaces();
    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int eidx = surf.data().faceEdges(i, j);
            edgeMetric[eidx] += C(i, j);
        }
    }

    // face mass matrix
    Eigen::VectorXd faceAreas;
    igl::doublearea(surf.data().V, surf.data().F, faceAreas);
    faceAreas *= 0.5;

    std::vector<Eigen::Triplet<double> > Mcoeffs;
    for (int i = 0; i < nfaces; i++)
    {
        Mcoeffs.push_back(Eigen::Triplet<double>(i, i, faceAreas[i]));
    }
    Eigen::SparseMatrix<double> M(nfaces, nfaces);
    M.setFromTriplets(Mcoeffs.begin(), Mcoeffs.end());

    std::cout << "Built mass matrices" << std::endl;

    // face Laplacian        
    std::vector<Eigen::Triplet<double> > DfaceCoeffs;
    for (int i = 0; i < nedges; i++)
    {
        int f0 = surf.data().E(i, 0);
        int f1 = surf.data().E(i, 1);
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

    // curl matrix
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
        
        DCoeffs.push_back(Eigen::Triplet<double>(idx, f0, -scaledvec0.dot(edgeVec)));
        DCoeffs.push_back(Eigen::Triplet<double>(idx, f1, scaledvec1.dot(edgeVec)));
        idx++;
    }
    Eigen::SparseMatrix<double> D(nconstraints, nfaces);
    D.setFromTriplets(DCoeffs.begin(), DCoeffs.end());

    Eigen::SparseMatrix<double> op = D.transpose() * D + sreg_ * Lface;

    std::cout << "Factoring op" << std::endl;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver(op);
    if(solver.info() != Eigen::Success)
        std::cout << "failed" << std::endl;

    std::cout << "Starting inverse power interation" << std::endl;
    Eigen::VectorXd x(nfaces);
    x.setRandom();
    x /= sqrt(x.transpose() * (M * x));
    for(int i=0; i<1000; i++)
    {
        Eigen::VectorXd newx = solver.solve(M*x);
        double xnorm = sqrt(newx.transpose() * (M*newx));
        x = newx/xnorm;            
    }

    std::cout << "Rayleigh quotient: " << x.transpose() * (op * x) / (x.transpose() * (M*x)) << std::endl;

    s.resize(nfaces);
    for (int i = 0; i < nfaces; i++)
        s[i] = x[i];

    // fix s sign
    double totals = s.sum();
    if (totals < 0)
        s *= -1.0;
}