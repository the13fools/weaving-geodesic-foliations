#include "CutMesh.h"
#include "CoverMesh.h"
#include "FieldSurface.h"
#include <queue>
#include <igl/writeOBJ.h>
#include <Eigen/Dense>
#include "Weave.h"
#include <igl/is_vertex_manifold.h>
#include <igl/is_edge_manifold.h>

typedef Eigen::Triplet<double> triplet;
# define M_PI           3.14159265358979323846

using namespace std;

CoverMesh::CoverMesh(const Weave &parent, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXi &oldToNewVertMap, const Eigen::MatrixXd &field, int ncovers)
 : parent_(parent), oldToNewVertMap_(oldToNewVertMap)
{
    fs = new FieldSurface(V, F, 1);
    int nfaces = F.rows();
    ncovers_ = ncovers;
    for (int i = 0; i < nfaces; i++)
    {
        fs->vectorFields.segment<2>(fs->vidx(i, 0)) = field.row(i).transpose();
    }

    theta.resize(fs->nVerts());
    theta.setZero();
    
    s.resize(fs->nFaces());
    s.setConstant(1.0);
    
    renderScale_ = 1.0;
}

CoverMesh::~CoverMesh()
{
    delete fs;
}

double CoverMesh::barycentric(double val1, double val2, double target)
{
    return (target-val1) / (val2-val1);
}

bool CoverMesh::crosses(double isoval, double val1, double val2, double minval, double maxval, double &bary)
{
    double halfperiod = 0.5*(maxval-minval);
    if(fabs(val2-val1) <= halfperiod)
    {
        bary = barycentric(val1, val2, isoval);
        if(bary >= 0 && bary < 1)
            return true;
        return false;
    }
    if(val1 < val2)
    {
        double wrapval1 = val1 + (maxval - minval);
        bary = barycentric(wrapval1, val2, isoval);
        if(bary >= 0 && bary < 1)
            return true;
        double wrapval2 = val2 + (minval - maxval);
        bary = barycentric(val1, wrapval2, isoval);
        if(bary >= 0 && bary < 1)
            return true;
    }
    else
    {
        double wrapval1 = val1 + (minval - maxval);
        bary = barycentric(wrapval1, val2, isoval);
        if(bary >= 0 && bary < 1)
            return true;
        double wrapval2 = val2 + (maxval - minval);
        bary = barycentric(val1, wrapval2, isoval);
        if(bary >= 0 && bary < 1)
            return true;
    }
    return false;
} 


int CoverMesh::extractIsoline(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXi &faceNeighbors, const Eigen::VectorXd &func, double isoval, double minval, double maxval)
{
    int nfaces = F.rows();
    bool *visited = new bool[nfaces];
    for(int i=0; i<nfaces; i++)
        visited[i] = false;

    int ntraces = 0;

    // Iterate between faces until encountering a zero level set.  
    // Trace out the level set in both directions from this face (marking faces as visited)
    // Save isoline to centerline
    for(int i=0; i<nfaces; i++)
    {
        if(visited[i])
            continue;
        visited[i] = true;
        std::vector<std::vector<Eigen::Vector3d> > traces;
        std::vector<std::vector<Eigen::Vector3d> > normals;
        std::vector<std::vector<std::pair<int, int> > > traces_vids;
        for(int j=0; j<3; j++)
        {
            int vp1 = F(i, (j+1)%3);
            int vp2 = F(i, (j+2)%3);
            double bary;
            if(crosses(isoval, func[vp1], func[vp2], minval, maxval, bary))
            {
                std::vector<Eigen::Vector3d> trace;
                std::vector<Eigen::Vector3d> norm;
                std::vector<std::pair<int, int> > trace_vid;
                trace.push_back( (1.0 - bary) * V.row(vp1) + bary * V.row(vp2) );
                trace_vid.push_back(std::make_pair(vp1, vp2));
                int prevface = i;
                int curface = faceNeighbors(i, j);
                //         norm.push_back(faceNormal(prevface));
                while(curface != -1 && !visited[curface])
                {
                    visited[curface] = true;
                    for(int k=0; k<3; k++)
                    {
                        if(faceNeighbors(curface, k) == prevface)
                            continue;
                        int vp1 = F(curface, (k+1)%3);
                        int vp2 = F(curface, (k+2)%3);
                        double bary;
                        if(crosses(isoval, func[vp1], func[vp2], minval, maxval, bary))
                        {
                            trace.push_back( (1.0 - bary) * V.row(vp1) + bary * V.row(vp2) );
                            trace_vid.push_back(std::make_pair(vp1, vp2));
                            norm.push_back(fs->faceNormal(curface));
                            prevface = curface;
                            curface = faceNeighbors(curface, k);
                            break;
                        }                       
                    }
                }
                traces.push_back(trace);
                normals.push_back(norm);
                traces_vids.push_back(trace_vid);
            }
        }
        assert(traces.size() == traces_vids.size());
        assert(traces.size() < 3);


        // if(traces.size() == 2)
        // {
        //     ntraces++;
        //     vector<Eigen::Vector3d> curISONormal;
        //     vector<Eigen::Vector3d> curISOLine;
        //     int nterms = traces[0].size() + traces[1].size();
        //     std::cout << nterms << " 0 0 " << nterms << " 0 0 " << std::endl;
        //     for(int j=0; j < traces[0].size(); j++)
        //     {
        //         curISOLine.push_back(traces[0][j]);
        //     }
        //     isoNormal.push_back(normals[0]);
        //     isoLines.push_back(curISOLine);
        //     std::cout << "trace size is 2\n";

        // }
        if(traces.size() == 1)
        {
            ntraces++;
            vector<Eigen::Vector3d> curISONormal;
            vector<Eigen::Vector3d> curISOLine;
            std::cout << traces[0].size() << " 0 0 " << traces[0].size() << " 0 0 " << std::endl;
            int next_vid1, next_vid2;
            for(int j=0; j<traces[0].size(); j++)
            {
                curISOLine.push_back(traces[0][j]);
                std::cout << traces[0][j].transpose() << " ";
                if (j == traces[0].size()-1)
                {
                    std::cout << " 0 0 0 " << std::endl;            
                    break;
                }
                else
                {
                    int next_vid1 = std::get<0>(traces_vids[0][j+1]);
                    int next_vid2 = std::get<1>(traces_vids[0][j+1]);
                    int cur_vid1 = std::get<0>(traces_vids[0][j]);
                    int cur_vid2 = std::get<1>(traces_vids[0][j]);
                    Eigen::Vector3d e1 = V.row(next_vid1) - V.row(next_vid2);
                    Eigen::Vector3d e2 = V.row(cur_vid1) - V.row(cur_vid2);
                    Eigen::Vector3d normal = (e1.cross(e2));
                    normal = normal / normal.norm();
                    std::cout << normal.transpose() << std::endl;
                    curISONormal.push_back(normal);
                }
            }
            isoNormal.push_back(curISONormal);
            isoLines.push_back(curISOLine);
            std::cerr << "trace size is 1\n";
        }
        if(traces.size() == 2)
        {
            ntraces++;
            vector<Eigen::Vector3d> curISONormal;
            vector<Eigen::Vector3d> curISOLine;
            int nterms = traces[0].size() + traces[1].size();
            std::cout << nterms << " 0 0 " << nterms << " 0 0 " << std::endl;
            for(int j=traces[1].size()-1; j >= 0; j--)
            {
                curISOLine.push_back(traces[1][j]);
                std::cout << traces[1][j].transpose() << " ";
                if (j == 0)
                {
                    int next_vid1 = std::get<0>(traces_vids[0][0]);
                    int next_vid2 = std::get<1>(traces_vids[0][0]);
                    int cur_vid1 = std::get<0>(traces_vids[1][j]);
                    int cur_vid2 = std::get<1>(traces_vids[1][j]);
                    Eigen::Vector3d e1 = V.row(next_vid1) - V.row(next_vid2);
                    Eigen::Vector3d e2 = V.row(cur_vid1) - V.row(cur_vid2);
                    Eigen::Vector3d normal = (e1.cross(e2));
                    normal = normal / normal.norm();
                    std::cout << normal.transpose() << std::endl;
                    curISONormal.push_back(normal);
                }
                else
                {
                    int next_vid1 = std::get<0>(traces_vids[1][j-1]);
                    int next_vid2 = std::get<1>(traces_vids[1][j-1]);
                    int cur_vid1 = std::get<0>(traces_vids[1][j]);
                    int cur_vid2 = std::get<1>(traces_vids[1][j]);
                    Eigen::Vector3d e1 = V.row(next_vid1) - V.row(next_vid2);
                    Eigen::Vector3d e2 = V.row(cur_vid1) - V.row(cur_vid2);
                    Eigen::Vector3d normal = (e1.cross(e2));
                    normal = normal / normal.norm();
                    std::cout << normal.transpose() << std::endl;
                    curISONormal.push_back(normal);
                }
            }
            for(int j=0; j<traces[0].size(); j++)
            {
                curISOLine.push_back(traces[0][j]);
                std::cout << traces[0][j].transpose() << " ";
                if (j == traces[0].size()-1)
                {
                    std::cout << "0 0 0 " << std::endl;
                    break;
                }
                else
                {
                    int next_vid1 = std::get<0>(traces_vids[0][j+1]);
                    int next_vid2 = std::get<1>(traces_vids[0][j+1]);
                    int cur_vid1 = std::get<0>(traces_vids[0][j]);
                    int cur_vid2 = std::get<1>(traces_vids[0][j]);
                    Eigen::Vector3d e1 = V.row(next_vid1) - V.row(next_vid2);
                    Eigen::Vector3d e2 = V.row(cur_vid1) - V.row(cur_vid2);
                    Eigen::Vector3d normal = (e1.cross(e2));
                    normal = normal / normal.norm();
                    std::cout << normal.transpose() << std::endl;
                    curISONormal.push_back(normal);
                }
            }
            isoNormal.push_back(curISONormal);
            isoLines.push_back(curISOLine);
            std::cout << "trace size is 2\n";
        }
    }
    delete[] visited;
    return ntraces;
}

void CoverMesh::drawISOLines(int numISOLines)
{
    double minval = -M_PI;
    double maxval = M_PI;
    double numlines = numISOLines;

    int nfaces = fs->nFaces();
    int nverts = fs->nVerts();

    std::map<std::pair<int, int>, Eigen::Vector2i > edgemap;
    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int nextj = (j + 1) % 3;
            int v1 = fs->data().F(i, j);
            int v2 = fs->data().F(i, nextj);
            int idx = 0;
            if (v1 > v2)
            {
                idx = 1;
                std::swap(v1, v2);
            }
            std::pair<int, int> p(v1,v2);
            std::map<std::pair<int, int>, Eigen::Vector2i >::iterator it = edgemap.find(p);
            if(it != edgemap.end())
                edgemap[p][idx] = i;
            else
            {
                Eigen::Vector2i entry(-1,-1);
                entry[idx] = i;
                edgemap[p] = entry;
            }
        }
    }

    Eigen::MatrixXi faceNeighbors(nfaces, 3);
    faceNeighbors.setConstant(-1);
    for(int i=0; i<nfaces; i++)
    {
        for(int j=0; j<3; j++)
        {
            int vp1 = fs->data().F(i,(j+1)%3);
            int vp2 = fs->data().F(i,(j+2)%3);
            if(vp1 > vp2) std::swap(vp1, vp2);
            std::map<std::pair<int, int>, Eigen::Vector2i >::iterator it = edgemap.find(std::pair<int,int>(vp1, vp2));
            if(it == edgemap.end())
                faceNeighbors(i, j) = -1;
            else
            {
                int opp = (it->second[0] == i ? it->second[1] : it->second[0]);
                faceNeighbors(i, j) = opp;
            }
        }
    }
    int ntraces = 0;
    isoLines.clear();
    isoNormal.clear();
    for(int i=0; i<numlines; i++)
    {
        double isoval = minval + (maxval-minval) * double(i)/double(numlines);
        ntraces += extractIsoline(fs->data().V, fs->data().F, faceNeighbors, 
            Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(theta.data(), theta.size()), 
            isoval, minval, maxval);
    }
    std::cout << ntraces << " 0 0 " << ntraces <<  " 0 0 " << std::endl;
}

void CoverMesh::computeFunc(double scalesInit)
{
    std::ofstream debugOut("debug.txt");
    std::ofstream debugVectsOut("debug.field");
    int nfaces = fs->nFaces();
    int nverts = fs->nVerts();
    cout << "nfaces: " << nfaces << endl;
    cout << "nverts: " << nverts << endl;
    vector<int> rowsL;
    vector<int> colsL;
    vector<double> difVecUnscaled;
    for (int fId = 0; fId < nfaces; fId++)
    { // Compute rowsL, colsL, difVecUnscaled
        int vId0 = fs->data().F(fId, 0);
        int vId1 = fs->data().F(fId, 1);
        int vId2 = fs->data().F(fId, 2);
        rowsL.push_back(vId0); rowsL.push_back(vId1); rowsL.push_back(vId2);
        colsL.push_back(vId1); colsL.push_back(vId2); colsL.push_back(vId0);
        Eigen::Vector3d p0 = fs->data().V.row(vId0);
        Eigen::Vector3d p1 = fs->data().V.row(vId1);
        Eigen::Vector3d p2 = fs->data().V.row(vId2);
        Eigen::Vector3d e01 = p0 - p1;
        Eigen::Vector3d e12 = p1 - p2;
        Eigen::Vector3d e20 = p2 - p0;
        Eigen::Vector3d faceVec;
        if (true)
        {
            faceVec = fs->data().Bs[fId] * fs->v(fId, 0); // The original vec            
        }
        faceVec = faceVec.cross(fs->faceNormal(fId));
        faceVec /= faceVec.norm();
        debugVectsOut << faceVec.transpose() << endl;
        difVecUnscaled.push_back(e01.dot(faceVec));
        difVecUnscaled.push_back(e12.dot(faceVec));
        difVecUnscaled.push_back(e20.dot(faceVec));
    }
    assert((rowsL.size() == 3 * nfaces) && (colsL.size() == 3 * nfaces) && (difVecUnscaled.size() == 3 * nfaces));

    // Eigen::SparseMatrix<double> faceLapMat = faceLaplacian();
    Eigen::VectorXd scales(nfaces);
    scales.setConstant(scalesInit);
    int totalIter = 6;
    for (int iter = 0; iter < totalIter; iter++)
    {
        vector<double> difVec;
        for (int i = 0; i < difVecUnscaled.size(); i++)
            difVec.push_back(difVecUnscaled[i] * scales(i / 3));
        std::vector<triplet> sparseContent;
        for (int i = 0; i < rowsL.size(); i++)
            sparseContent.push_back(triplet(rowsL[i], colsL[i], 1));
        Eigen::SparseMatrix<double> TP(nverts, nverts);
        Eigen::SparseMatrix<double> TPTran(nverts, nverts);
        TP.setFromTriplets(sparseContent.begin(), sparseContent.end());
        TPTran = TP.transpose();
        TP += TPTran;
        vector<int> degree;
        for (int i = 0; i < nverts; i++)
            degree.push_back(TP.row(i).sum());
        std::vector<triplet> AContent;
        for (int i = 0; i < rowsL.size(); i++)
        {
            double cVal = cos(difVec[i]);
            double sVal = sin(difVec[i]);
            AContent.push_back(triplet(2 * rowsL[i], 2 * colsL[i], cVal));
            AContent.push_back(triplet(2 * rowsL[i], 2 * colsL[i] + 1, -sVal));
            AContent.push_back(triplet(2 * rowsL[i] + 1, 2 * colsL[i], sVal));
            AContent.push_back(triplet(2 * rowsL[i] + 1, 2 * colsL[i] + 1, cVal));
        }
        Eigen::SparseMatrix<double> Amat(2 * nverts, 2 * nverts);
        Eigen::SparseMatrix<double> Amat_tran(2 * nverts, 2 * nverts);
        Amat.setFromTriplets(AContent.begin(), AContent.end());
        Amat_tran = Amat.transpose();
        Amat += Amat_tran;
        //
        std::vector<triplet> LContent;
        for (int i = 0; i < 2 * nverts; i++)
            LContent.push_back(triplet(i, i, degree[int(i / 2)]));
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
        double eigenVal = eigenVec.transpose() * Lmat * eigenVec;
        cout << "Current iteration = " << iter << " currents error is: " << eigenVal << endl;
        // Extract the function value
        vector<double> curTheta;
        for (int i = 0; i < nverts; i++)
        {
            double curCos = eigenVec(2 * i);
            double curSin = eigenVec(2 * i + 1);
            double normalizer = sqrt(curCos * curCos + curSin * curSin);
            double curFunc = acos(curCos / normalizer);
            if (curSin < 0)
                curFunc = -curFunc;
            curTheta.push_back(curFunc);
        }
        ////
        //// Re-compute face scales
        vector<double> difVecPred;
        for (int i = 0; i < rowsL.size(); i++)
        {
            double curPred = curTheta[rowsL[i]] - curTheta[colsL[i]];
            if (curPred > M_PI) curPred -= 2 * M_PI;
            if (curPred < -M_PI) curPred += 2 * M_PI;
            difVecPred.push_back(curPred);
        }
        Eigen::VectorXd bScales(nfaces);
        vector<double> diagAScales;
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
        std::vector<triplet> AScalesContent;
        for (int i = 0; i < nfaces; i++)
            AScalesContent.push_back(triplet(i, i, diagAScales[i]));
        Eigen::SparseMatrix<double> AScalesMat(nfaces, nfaces);
        AScalesMat.setFromTriplets(AScalesContent.begin(), AScalesContent.end());
        // Solve for scale
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solverScales(AScalesMat);
        Eigen::VectorXd curScales = solverScales.solve(bScales);
        for (int i = 0; i < nfaces; i++)
            scales(i) = curScales(i);
        for(int i=0; i<nverts; i++)
            theta[i] = curTheta[i];
    }
    for (int i = 0; i < nverts; i++)
        debugOut << theta[i] << endl;
    debugOut.close();
}

Eigen::SparseMatrix<double> CoverMesh::faceLaplacian()
{ // Only augment vector
    int nfaces = fs->nFaces();
    // TODO: boundary
    // ids = find(min(adjFaces) > 0);
    // adjFaces = adjFaces(:, ids);
    std::vector<triplet> AContent;
    for (int i = 0; i < fs->data().E.rows(); i ++)
        AContent.push_back(triplet(fs->data().E(i,0), fs->data().E(i,1), 1));
    Eigen::SparseMatrix<double> AFaceMat (nfaces, nfaces);
    AFaceMat.setFromTriplets(AContent.begin(),AContent.end());
    // Ge the degree of face
    vector<int> degreeFace;
    for (int i = 0; i < nfaces; i ++)
        degreeFace.push_back(AFaceMat.row(i).sum());
    // Get LFace
    std::vector<triplet> LContent;
    for (int i = 0; i < degreeFace.size(); i ++)
        LContent.push_back(triplet(i, i, degreeFace[i]));
    Eigen::SparseMatrix<double> LFaceMat (nfaces, nfaces);
    LFaceMat.setFromTriplets(LContent.begin(), LContent.end());
    LFaceMat -= AFaceMat;
    // Eigen::SparseMatrix<double> LFaceMat;
    return LFaceMat;
}


void CoverMesh::createVisualization(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &edgePts, Eigen::MatrixXd &edgeVecs, Eigen::MatrixXi &edgeSegs, Eigen::MatrixXd &colors)
{
    int facespercover = fs->nFaces()/ncovers_;
    int rows = 2;
    int meshesperrow = ncovers_ / rows + (ncovers_ % rows == 0 ? 0 : 1);
    std::vector<Eigen::Vector3d> offsets;
    for(int i=0; i<ncovers_; i++)
    {
        int row = i/meshesperrow;
        int col = i%meshesperrow;
        double dy = (-1.0 * row + (1.0) * (rows-row-1))/double(rows);
        double dx = (1.0 * col + (-1.0) * (meshesperrow - col-1))/double(meshesperrow);
        offsets.push_back(Eigen::Vector3d(dx, dy, 0.0));
    }
    
    int origverts = parent_.fs->nVerts();
    int origfaces = parent_.fs->nFaces();
    int newverts = ncovers_*origverts;
    int newfaces = ncovers_*origfaces;
    V.resize(newverts, 3);
    F.resize(newfaces, 3);
    renderScale_ = 1.0 / std::max(rows,meshesperrow);
    for(int i=0; i<ncovers_; i++)
    {
        for(int j=0; j<origverts; j++)
        {
            V.row(i*origverts + j) = offsets[i].transpose() + renderScale_ * parent_.fs->data().V.row(j);
        }
        for(int j=0; j<origfaces; j++)
        {
            for(int k=0; k<3; k++)
            {
                F(i*origfaces + j, k) = i*origverts + parent_.fs->data().F(j,k);
            }
        }
    }
    edgePts.resize(newfaces, 3);
    edgeVecs.resize(newfaces, 3);
    edgeVecs.setZero();
    edgeSegs.resize(newfaces, 2);
    colors.resize(newfaces , 3);

    for(int c=0; c<ncovers_; c++)
    {
        for (int i = 0; i < origfaces; i++)
        {
            Eigen::Vector3d centroid;
            centroid.setZero();
            for (int j = 0; j < 3; j++)
                centroid += renderScale_ * parent_.fs->data().V.row(parent_.fs->data().F(i, j));
            centroid /= 3.0;
            centroid += offsets[c];

            edgePts.row(c*origfaces + i) = centroid;
            edgeVecs.row(c*origfaces + i) = parent_.fs->data().Bs[i] * fs->v(c*origfaces + i, 0);
            edgeSegs(c*origfaces + i, 0) = 2 * (c*origfaces + i);
            edgeSegs(c*origfaces + i, 1) = 2 * (c*origfaces + i) + 1;
            colors.row(c*origfaces + i) = Eigen::Vector3d(0,0,0).transpose();
        }
    }
}

int CoverMesh::visMeshToCoverMesh(int vertid)
{
    return oldToNewVertMap_[vertid];
}

void CoverMesh::initializeS()
{
    Eigen::VectorXi B;
    if(!igl::is_vertex_manifold(fs->data().F, B))
        std::cout << "ERROR: cover mesh not vertex-manifold!" << std::endl;
    if(!igl::is_edge_manifold(fs->data().F))
        std::cout << "ERROR: cover mesh not edge-manifold!" << std::endl;
    // cut the mesh
    std::vector<std::vector<int> > cuts;
    findCuts(fs->data().V, fs->data().F, cuts);
    std::cout << "Found " << cuts.size() << " cuts";
    if (cuts.size() == 0)
        std::cout << std::endl;
    else
    {
        std::cout << " of lengths ";
        for (int i = 0; i < cuts.size(); i++)
        {
            if (i != 0)
                std::cout << ", ";
            std::cout << cuts[i].size();
        }
        std::cout << std::endl;
    }

    Eigen::MatrixXd cutV;
    Eigen::MatrixXi cutF;    
    cutMesh(fs->data().V, fs->data().F, cuts, cutV, cutF);
    Surface surf(cutV, cutF);
    
    // build edge metric matrix
    std::vector<Eigen::Triplet<double> > edgeMetricCoeffs;
    int nedges = surf.nEdges();
    for(int i=0; i<nedges; i++)
    {
        Eigen::Vector3d v0 = surf.data().V.row(surf.data().edgeVerts(i,0)).transpose();
        Eigen::Vector3d v1 = surf.data().V.row(surf.data().edgeVerts(i,1)).transpose();
        double len = (v1-v0).norm();
        edgeMetricCoeffs.push_back(Eigen::Triplet<double>(i,i,len));
    }
    Eigen::SparseMatrix<double> edgeMetric(nedges, nedges);
    edgeMetric.setFromTriplets(edgeMetricCoeffs.begin(), edgeMetricCoeffs.end());
    
    double reg = 1e-4;
    
    // edge gradient matrix
    std::vector<Eigen::Triplet<double> > DCoeffs;
    int nfaces = surf.nFaces();
    for(int i=0; i<nedges; i++)
    {
        int f0 = surf.data().E(i,0);
        int f1 = surf.data().E(i,1);
        if(f0 == -1 || f1 == -1)
            continue;
        Eigen::Vector3d v0 = surf.data().V.row(surf.data().edgeVerts(i,0)).transpose();
        Eigen::Vector3d v1 = surf.data().V.row(surf.data().edgeVerts(i,1)).transpose();
        Eigen::Vector3d edgeVec = v1-v0;
        Eigen::Vector3d scaledvec0 = surf.data().Bs[f0] * fs->v(f0, 0);           
        Eigen::Vector3d scaledvec1 = surf.data().Bs[f1] * fs->v(f1, 0);
        DCoeffs.push_back(Eigen::Triplet<double>(i, f0, -scaledvec0.dot(edgeVec) - reg));
        DCoeffs.push_back(Eigen::Triplet<double>(i, f1, scaledvec1.dot(edgeVec) + reg));
    }
    Eigen::SparseMatrix<double> D(nedges, nfaces);
    D.setFromTriplets(DCoeffs.begin(), DCoeffs.end());
    // the integrability operator
    Eigen::SparseMatrix<double> L = D.transpose() * edgeMetric * D;
    std::cout << "Solving eigenproblem..." << std::endl;
    double eval = inversePowerIteration(L, s, 1000);
    std::cout << "Smallest eigenvalue: " << eval << std::endl;
}

double CoverMesh::inversePowerIteration(Eigen::SparseMatrix<double> &M, Eigen::VectorXd &evec, int iters)
{
    evec.resize(M.cols());
    evec.setRandom();
    evec /= evec.norm();
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver(M);
    for(int i=0; i<iters; i++)
    {
        Eigen::VectorXd newvec = solver.solve(evec);
        evec = newvec / newvec.norm();
    }
    return evec.transpose() * M * evec;
}
