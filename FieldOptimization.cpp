#include "FieldOptimization.h"
#include <map>
#include <iostream>
#include "FaceBased.h"

#include <Eigen/Geometry>

using namespace std;
using namespace Eigen;

void buildEdges(const Eigen::MatrixXi &F, Eigen::MatrixXi &E)
{
    map<pair<int, int>, Vector4i, std::less<pair<int, int> >,
        Eigen::aligned_allocator<std::pair<const int, Eigen::Vector4i> >>
        edgemap;

    int nfaces = (int)F.rows();
    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int idx1 = F(i, j);
            int idx2 = F(i, (j + 1) % 3);
            int slot = 0;

            if (idx1 > idx2)
            {
                swap(idx1, idx2);
                slot = 1;
            }
            map<pair<int, int>, Vector4i, std::less<pair<int, int> >,
                Eigen::aligned_allocator<std::pair<const int, Eigen::Vector4i> >>::iterator it = edgemap.find(pair<int, int>(idx1, idx2));
            if (it == edgemap.end())
            {
                Vector4i newedge;
                newedge[0] = idx1;
                newedge[1] = idx2;
                newedge[2] = newedge[3] = -1;
                newedge[2 + slot] = i;
                edgemap[pair<int, int>(idx1, idx2)] = newedge;
            }
            else
            {
                edgemap[pair<int, int>(idx1, idx2)][2 + slot] = i;
            }
        }
    }

    int nedges = (int)edgemap.size();
    E.resize(nedges, 4);
    int idx = 0;
    for (map<pair<int, int>, Vector4i>::iterator it = edgemap.begin(); it != edgemap.end(); ++it)
    {
        E.row(idx) = it->second.transpose();
        idx++;
    }
}


bool consistencyCheckEdges(const Eigen::MatrixXi &F, const Eigen::MatrixXi &E, const Eigen::MatrixXi &F_edges)
{
    int nfaces = (int)F.rows();
    for (int i = 0; i < nfaces; i++)
    {
        for (int e = 0; e < 3; e++)
        {
            //each edge had better point back to the face
            int edgeidx = F_edges(i, e);
            if (E(edgeidx, 2) != i && E(edgeidx, 3) != i)
            {
                std::cerr << "Edge assigned to triangle " << i << " side " << e << " does not have triangle " << i << " as an adjacent face!" << std::endl;
                return false;
            }
            //the edge endpoints need to be opposite vertex e of triangle i
            for (int vtex = 0; vtex < 2; vtex++)
            {
                int facevertidx = (e + 1 + vtex) % 3;
                int vertidx = F(i, facevertidx);
                if (E(edgeidx, 0) != vertidx && E(edgeidx, 1) != vertidx)
                {
                    std::cerr << "Vertex " << vertidx << " is opposite vertex " << e << " on triangle " << i << " but is not part of the edge assigned opposite vertex " << e << std::endl;
                    return false;
                }
            }
        }
    }
    return true;
}

// We assume that the UV coordinates are such that the edge u := 10, v := 20 per face
// From this, we have that the 0 edge is the one not associated with u or v
// with u as the 1 edge and v as the 2 edge
void buildEdgesPerFace(const Eigen::MatrixXi &F, const Eigen::MatrixXi &E, Eigen::MatrixXi &F_edges)
{
    int nfaces = (int)F.rows();
    int nedges = (int)E.rows();

    F_edges.resize(nfaces, 3);
    F_edges = Eigen::MatrixXi::Constant(nfaces,3,-1);
    for (int i = 0; i < nedges; i++)
    {
        int f1 = E(i, 2);
        int f2 = E(i, 3);

        if (f1 > -1)
        {
            int insIdx = -1;
            for (int j = 0; j < 3; j++)
            {
                if (F(f1, j) == E(i, 1))
                {
                    insIdx = (j+1)%3;
                }
            }
            F_edges(f1, insIdx) = i;
        }

        if (f2 > -1)
        {
            int insIdx = -1;
            for (int j = 0; j < 3; j++)
            {
                if (F(f2, j) == E(i, 0))
                {
                    insIdx = (j+1)%3;
                }
            }
            F_edges(f2, insIdx) = i;
        }
    }
    if (!consistencyCheckEdges(F, E, F_edges))
    {
        assert(false);
        exit(-1);
    }
}

void computeJs(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, std::vector<Eigen::Matrix3d> &Js)
{
    int nfaces = (int)F.rows();
    Js.resize(nfaces);
    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector3d v0 = V.row(F(i, 0));
        Eigen::Vector3d v1 = V.row(F(i, 1));
        Eigen::Vector3d v2 = V.row(F(i, 2));
        Eigen::Vector3d n = (v1 - v0).cross(v2 - v0);
        n /= n.norm();
        Js[i] << 0, -n[2], n[1],
            n[2], 0, -n[0],
            -n[1], n[0], 0;
    }
}

void initOptVars(const Eigen::MatrixXd &v0,                 
	         OptVars &vars)
{
    vars.V_opt = v0;
    vars.W_opt = v0;

    // unroll optimization variables
    int nfaces = (int)v0.rows();
    vars.v.resize(3 * nfaces);
    vars.w.resize(3 * nfaces);
    for (int i = 0; i < nfaces; i++)
    {
        vars.v.segment<3>(3 * i) = v0.row(i).transpose();
        vars.w.segment<3>(3 * i) = v0.row(i).transpose();
    }
    vars.D.resize(9 * nfaces);
    vars.D.setZero();
    
}


void computeCentroids(const Eigen::MatrixXi &F,const Eigen::MatrixXd &V, Eigen::MatrixXd &centroids)
{
    //   Eigen::MatrixXd interp; 
    int nfaces = F.rows();
    int nverts = V.rows();

    centroids.resize(nfaces, 3);
    for (int i = 0; i < nfaces; i++) 
    { 
        Eigen::Vector3d pos(0,0,0);
        for (int j = 0; j < 3; j++) 
        {
            //	    std::cout << V.row(F(i,j)) << "\n" << F(i,j) << "\n\n";
            pos += V.row(F(i,j));
        }
        centroids.row(i) = pos/3;
    }
}

MeshData::MeshData(const Eigen::MatrixXd &V, 
	           const Eigen::MatrixXi &F) : V(V), F(F)
{
    buildEdges(F, E);
    buildEdgesPerFace(F, E, F_edges);
    computeJs(F, V, Js);
    computeCentroids(F,V,centroids_F);
}

double energy(const OptVars &vars, const MeshData &mesh, double lambda, double mu)
{
    double result = 0;
    int nfaces = (int)mesh.F.rows();
    for (int i = 0; i < nfaces; i++)
    {
        result += 0.5 * (mesh.v0.row(i).transpose() - vars.v.segment<3>(3 * i)).squaredNorm();
        Matrix3d Di;
        for (int j = 0; j < 3; j++)
        {
            Di.row(j) = vars.D.segment<3>(9 * i + 3 * j);
        }
        Eigen::Vector3d wi = vars.w.segment<3>(3 * i);
        Eigen::Vector3d vi = vars.v.segment<3>(3 * i);
        double term = (wi.transpose() * mesh.Js[i].transpose() * Di.transpose() * vi);
        result += 0.5 * lambda * term*term;
    }

    int nedges = (int)mesh.E.rows();
    for(int i=0; i<nedges; i++)       
    { 
        int face1 = mesh.E(i, 2);
        int face2 = mesh.E(i, 3);
        if (face1 == -1 || face2 == -1)
            continue;
        Matrix3d D1;
        for (int j = 0; j < 3; j++)
        {
            D1.row(j) = vars.D.segment<3>(9 * face1 + 3 * j);
        }
        Matrix3d D2;
        for (int j = 0; j < 3; j++)
        {
            D2.row(j) = vars.D.segment<3>(9 * face2 + 3 * j);
        }
        Eigen::Vector3d cdiff = mesh.centroids_F.row(face2) - mesh.centroids_F.row(face1);
        Eigen::Vector3d vdiff = vars.v.segment<3>(3 * face2) - vars.v.segment<3>(3 * face1);
        result += 0.5 * mu * (D1*cdiff - vdiff).squaredNorm();
        result += 0.5 * mu * (D2*cdiff - vdiff).squaredNorm();
    }

    result += 0.5 * mu * (vars.v - vars.w).squaredNorm();
    result += 0.5 * mu * (vars.D).squaredNorm();
    return result;
}

// computes energy derivative, of the form M * v + b
void dvEnergy(const OptVars &vars, const MeshData &mesh, double lambda, double mu, Eigen::SparseMatrix<double> &M, Eigen::VectorXd &b)
{
    int nfaces = (int)mesh.F.rows();
    M.resize(3 * nfaces, 3 * nfaces);
    M.setIdentity();
    M*= (1 + mu);
    
    vector<Triplet<double> > lambdatermcoeffs;
    for (int i = 0; i < nfaces; i++)
    {
        Matrix3d Di;
        for (int j = 0; j < 3; j++)
        {
            Di.row(j) = vars.D.segment<3>(9 * i + 3 * j);
        }
        Vector3d term = Di * mesh.Js[i] * vars.w.segment<3>(3 * i);
        Matrix3d termsq = term * term.transpose();
        for(int j=0; j<3; j++)
            for (int k = 0; k < 3; k++)
            {
                lambdatermcoeffs.push_back(Triplet<double>(3 * i + j, 3 * i + k, termsq(j, k)));
            }
    }
    SparseMatrix<double> lambdaterm(3 * nfaces, 3 * nfaces);
    lambdaterm.setFromTriplets(lambdatermcoeffs.begin(), lambdatermcoeffs.end());
    M += lambda * lambdaterm;

    vector<Triplet<double> > compattermcoeffs;
    int nedges = (int)mesh.E.rows();
    for (int i = 0; i < nedges; i++)
    {
        int face1 = mesh.E(i, 2);
        int face2 = mesh.E(i, 3);
        if (face1 == -1 || face2 == -1)
            continue;

        for (int j = 0; j < 3; j++)
        {
            compattermcoeffs.push_back(Triplet<double>(3 * face1 + j, 3 * face1 + j, 2.0));
            compattermcoeffs.push_back(Triplet<double>(3 * face2 + j, 3 * face2 + j, 2.0));
            compattermcoeffs.push_back(Triplet<double>(3 * face1 + j, 3 * face2 + j, -2.0));
            compattermcoeffs.push_back(Triplet<double>(3 * face2 + j, 3 * face1 + j, -2.0));
        }

    }
    SparseMatrix<double> compatmat(3 * nfaces, 3 * nfaces);
    compatmat.setFromTriplets(compattermcoeffs.begin(), compattermcoeffs.end());
    M += mu*compatmat;

    b.resize(3 * nfaces);
    b.setZero();
    for (int i = 0; i < nedges; i++)
    {
        int face1 = mesh.E(i, 2);
        int face2 = mesh.E(i, 3);
        if (face1 == -1 || face2 == -1)
            continue;
        Matrix3d D1;
        for (int j = 0; j < 3; j++)
        {
            D1.row(j) = vars.D.segment<3>(9 * face1 + 3 * j);
        }
        Matrix3d D2;
        for (int j = 0; j < 3; j++)
        {
            D2.row(j) = vars.D.segment<3>(9 * face2 + 3 * j);
        }
        Eigen::Vector3d cdiff = mesh.centroids_F.row(face2) - mesh.centroids_F.row(face1);
        b.segment<3>(3 * face2) += mu*D1*cdiff;
        b.segment<3>(3 * face1) -= mu*D1*cdiff;
        b.segment<3>(3 * face1) -= mu*D2*cdiff;
        b.segment<3>(3 * face2) += mu*D2*cdiff;
    }

    for (int i = 0; i < nfaces; i++)
    {
        b.segment<3>(3 * i) += mesh.v0.row(i).transpose();
    }
    b += mu * vars.w;
}

// computes energy derivative, of the form M * w + b
void dwEnergy(const OptVars &vars, const MeshData &mesh, double lambda, double mu, Eigen::SparseMatrix<double> &M, Eigen::VectorXd &b)
{
    int nfaces = (int)mesh.F.rows();
    M.resize(3 * nfaces, 3 * nfaces);
    M.setIdentity();
    M*= mu;
    
    vector<Triplet<double> > lambdatermcoeffs;
    for (int i = 0; i < nfaces; i++)
    {
        Matrix3d Di;
        for (int j = 0; j < 3; j++)
        {
            Di.row(j) = vars.D.segment<3>(9 * i + 3 * j);
        }
        Vector3d term = mesh.Js[i].transpose() * Di.transpose() * vars.v.segment<3>(3 * i);
        Matrix3d termsq = term * term.transpose();
        for(int j=0; j<3; j++)
            for (int k = 0; k < 3; k++)
            {
                lambdatermcoeffs.push_back(Triplet<double>(3 * i + j, 3 * i + k, termsq(j, k)));
            }
    }
    SparseMatrix<double> lambdaterm(3 * nfaces, 3 * nfaces);
    lambdaterm.setFromTriplets(lambdatermcoeffs.begin(), lambdatermcoeffs.end());
    M += lambda * lambdaterm;

    b.resize(3 * nfaces);
    b.setZero();

    b += mu * vars.v;
}

// computes energy derivative, of the form M * D + b
void dDEnergy(const OptVars &vars, const MeshData &mesh, double lambda, double mu, Eigen::SparseMatrix<double> &M, Eigen::VectorXd &b)
{
    int nfaces = (int)mesh.F.rows();
    M.resize(9 * nfaces, 9 * nfaces);
    M.setIdentity();
    M *= mu;
    
    vector<Triplet<double> > lambdatermcoeffs;
    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector3d vi = vars.v.segment<3>(3 * i);
        Eigen::Vector3d wi = vars.w.segment<3>(3 * i);
        Matrix3d term = vi * (mesh.Js[i] * wi).transpose();
        VectorXd flattened(9);
        for (int j = 0; j < 3; j++)
            flattened.segment<3>(3 * j) = term.row(j);

        MatrixXd flsq = flattened * flattened.transpose();

        for(int j=0; j<9; j++)
            for (int k = 0; k < 9; k++)
            {
                lambdatermcoeffs.push_back(Triplet<double>(9 * i + j, 9 * i + k, flsq(j, k)));
            }
    }
    SparseMatrix<double> lambdaterm(9 * nfaces, 9 * nfaces);
    lambdaterm.setFromTriplets(lambdatermcoeffs.begin(), lambdatermcoeffs.end());
    M += lambda * lambdaterm;

    b.resize(9 * nfaces);
    b.setZero();

    std::vector<Eigen::Triplet<double> > compattermcoeffs;
    int nedges = (int)mesh.E.rows();
    for (int i = 0; i < nedges; i++)
    {
        int face1 = mesh.E(i, 2);
        int face2 = mesh.E(i, 3);
        if (face1 == -1 || face2 == -1)
            continue;
        Eigen::Vector3d cdiff = mesh.centroids_F.row(face2) - mesh.centroids_F.row(face1);
        Eigen::MatrixXd cmat(9, 3);
        cmat.setZero();
        for (int j = 0; j < 3; j++)
        {
            cmat(j, 0) = cdiff[j];
            cmat(3 + j, 1) = cdiff[j];
            cmat(6 + j, 2) = cdiff[j];
        }
        Eigen::MatrixXd cTc = cmat*cmat.transpose();
        for(int j=0; j<9; j++)
            for (int k = 0; k < 9; k++)
            {
                compattermcoeffs.push_back(Eigen::Triplet<double>(9 * face1 + j, 9 * face1 + k, cTc(j, k)));
                compattermcoeffs.push_back(Eigen::Triplet<double>(9 * face2 + j, 9 * face2 + k, cTc(j, k)));
            }
        Eigen::Vector3d vdiff = vars.v.segment<3>(3 * face2) - vars.v.segment<3>(3 * face1);
        b.segment<9>(9 * face1) += mu * cmat * vdiff;
        b.segment<9>(9 * face2) += mu * cmat * vdiff;
    }

    SparseMatrix<double> compatmat(9 * nfaces, 9 * nfaces);
    compatmat.setFromTriplets(compattermcoeffs.begin(), compattermcoeffs.end());
    M += mu * compatmat;
}

void checkFiniteDifferences(const OptVars &vars, const MeshData &mesh, double lambda, double mu)
{
    int nfaces = (int)mesh.F.rows();
    OptVars test = vars;
    VectorXd noise1(3 * nfaces);
    noise1.setRandom();
    test.v += 0.01 * noise1;
    VectorXd noise2(3 * nfaces);
    noise2.setRandom();
    test.w += 0.01 * noise2;

    double orige = energy(test, mesh, lambda, mu);

    int idx = 3;
    OptVars perturbed = test;
    perturbed.v[idx] += 1e-6;
    double newe = energy(perturbed, mesh, lambda, mu);
    Eigen::SparseMatrix<double> M;
    Eigen::VectorXd b;
    dvEnergy(test, mesh, lambda, mu, M, b);
    Eigen::VectorXd deriv = M * test.v - b;
    double findiff = (newe - orige) / 1e-6;
    std::cout << "Finite difference check v " << findiff << " vs " << deriv[idx] << std::endl;

    perturbed = test;
    perturbed.w[idx] += 1e-6;
    newe = energy(perturbed, mesh, lambda, mu);
    dwEnergy(test, mesh, lambda, mu, M, b);
    deriv = M * test.w - b;
    findiff = (newe - orige) / 1e-6;
    std::cout << "Finite difference check w " << findiff << " vs " << deriv[idx] << std::endl;

    perturbed = test;
    perturbed.D[idx] += 1e-6;
    newe = energy(perturbed, mesh, lambda, mu);
    dDEnergy(test, mesh, lambda, mu, M, b);
    deriv = M * test.D - b;
    findiff = (newe - orige) / 1e-6;
    std::cout << "Finite difference check D " << findiff << " vs " << deriv[idx] << std::endl;

}


void rollupOptVars(OptVars &vars)
{
    // rollup optimization variables
    int nfaces = (int)vars.W_opt.rows();
    for (int i = 0; i < nfaces; i++)
    {
       for (int j = 0; j < 3; j++)
       {
           vars.V_opt(i, j) = vars.v(3*i + j);
           vars.W_opt(i, j) = vars.w(3*i + j);
       }
 //   	vars.V_opt.row(i) = vars.v.segment<3>(3 * i);//.transpose(); 
 //       vars.W_opt.row(i) = vars.w.segment<3>(3 * i);//.transpose(); 
    }
    
}

void alternatingMinimization(const MeshData &mesh, double lambda, double mu, OptVars &vars)
{
    checkFiniteDifferences(vars, mesh, lambda, mu);
    
    double orig = energy(vars, mesh, lambda, mu);
    std::cout << "Original energy: " << orig << std::endl;
    SparseMatrix<double> M;
    VectorXd b;
    dvEnergy(vars, mesh, lambda, mu, M, b);
    SimplicialLDLT<SparseMatrix<double> > solver;
    solver.compute(M);
    VectorXd newv = solver.solve(b);
    std::cout << "Residual: " << (M*newv - b).norm() << std::endl;
    vars.v = newv;
    double newe = energy(vars, mesh, lambda, mu);
    std::cout << "After v: " << newe << std::endl;

    dwEnergy(vars, mesh, lambda, mu, M, b);
    solver.compute(M);
    VectorXd neww = solver.solve(b);
    std::cout << "Residual: " << (M*neww - b).norm() << std::endl;
    vars.w = neww;
//    std::cout << neww << std::endl;
    newe = energy(vars, mesh, lambda, mu);
    std::cout << "After w: " << newe << std::endl;

    dDEnergy(vars, mesh, lambda, mu, M, b);
    solver.compute(M);
    VectorXd newD = solver.solve(b);
    std::cout << "Residual: " << (M*newD - b).norm() << std::endl;
    vars.D = newD;
    newe = energy(vars, mesh, lambda, mu);
    std::cout << "After D: " << newe << std::endl;

    rollupOptVars(vars);
}
