#include "CutMesh.h"
#include "CoverMesh.h"
#include "FieldSurface.h"
#include <queue>
#include <igl/writeOBJ.h>
#include <Eigen/Dense>
#include "Weave.h"
#include <igl/is_vertex_manifold.h>
#include <igl/is_edge_manifold.h>
#include <set>
#include <igl/cotmatrix_entries.h>
#include <igl/facet_components.h>
#include <igl/remove_unreferenced.h>
#include "Traces.h"
#include "igl/massmatrix.h"

typedef Eigen::Triplet<double> triplet;
# define M_PI           3.14159265358979323846

using namespace std;

CoverMesh::CoverMesh(const Weave &parent, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXi &oldToNewVertMap, const Eigen::MatrixXd &field, int ncovers)
 : parent_(parent)
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

    initializeSplitMesh(oldToNewVertMap);
}

CoverMesh::~CoverMesh()
{
    delete fs;
    if (data_.splitMesh)
        delete data_.splitMesh;
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


void CoverMesh::extractIsoline(const Eigen::VectorXd &func, double isoval, double minval, double maxval, std::vector<Trace> &isotrace)
{    
    int nfaces = fs->nFaces();
    bool *visited = new bool[nfaces];
    for(int i=0; i<nfaces; i++)
        visited[i] = false;
        

    // Iterate between faces until encountering a zero level set.  
    // Trace out the level set in both directions from this face (marking faces as visited)
    // Save isoline to centerline
    for(int i=0; i<nfaces; i++)
    {
        if(visited[i])
            continue;
        visited[i] = true;
        std::vector<int> crossings;
        std::vector<double> crossingsbary;

        std::vector<std::vector<TraceSegment> > tracePieces;

        for(int j=0; j<3; j++)
        {            
            int vp1 = fs->data().F(i, (j+1)%3);
            int vp2 = fs->data().F(i, (j+2)%3);
            double bary;
            if(crosses(isoval, func[vp1], func[vp2], minval, maxval, bary))
            {
                crossings.push_back(j);
                crossingsbary.push_back(bary);
                std::vector<TraceSegment> trace;                
                
                int prevface = i;
                int curface = fs->data().faceNeighbors(i, j);
                while(curface != -1 && !visited[curface])
                {                
                    visited[curface] = true;
                    TraceSegment nextseg;
                    nextseg.face = curface;
                    for(int k=0; k<3; k++)
                    {
                        if(fs->data().faceNeighbors(curface, k) == prevface)
                        {
                            nextseg.side[0] = k;
                            nextseg.bary[0] = 1.0-bary;
                            break;
                        }
                    }

                    for(int k=0; k<3; k++)
                    {
                        if(fs->data().faceNeighbors(curface, k) == prevface)
                            continue;
                        int vp1 = fs->data().F(curface, (k+1)%3);
                        int vp2 = fs->data().F(curface, (k+2)%3);
                        if(crosses(isoval, func[vp1], func[vp2], minval, maxval, bary))
                        {
                            nextseg.side[1] = k;
                            nextseg.bary[1] = bary;
                            trace.push_back(nextseg);
                            prevface = curface;
                            curface = fs->data().faceNeighbors(curface, k);
                            break;
                        }                       
                    }
                }
                tracePieces.push_back(trace);
            }
        }
        assert(tracePieces.size() < 3);


        if(tracePieces.size() == 1)
        {
            // lucky! no stitching together needed
            Trace line(fs, Trace_Mode::FIELD);
            line.segs = tracePieces[0];
            isotrace.push_back(line);
        }
        if(tracePieces.size() == 2)
        {
            // must stitch together both traces into one isoline
            Trace line(fs, Trace_Mode::FIELD);
            // first, reverse the order and orientation of the segments in traces[0]
            for(auto it = tracePieces[0].rbegin(); it != tracePieces[0].rend(); ++it)
            {
                TraceSegment rev = *it;
                std::swap(rev.side[0], rev.side[1]);
                std::swap(rev.bary[0], rev.bary[1]);
                line.segs.push_back(rev);
            }
            // add in the connecting segment
            TraceSegment con;
            con.face = i;
            con.side[0] = crossings[0];
            con.side[1] = crossings[1];
            con.bary[0] = crossingsbary[0];
            con.bary[1] = crossingsbary[1];
            line.segs.push_back(con);
            // finally append all of traces[1]
            for(auto &it : tracePieces[1])
                line.segs.push_back(it);
                
            isotrace.push_back(line);
        }
    }
    delete[] visited;
}

void  CoverMesh::recomputeIsolines(int numISOLines, std::vector<Trace> &isotraces)
{
    double minval = -M_PI;
    double maxval = M_PI;
    double numlines = numISOLines;

    isotraces.clear();

    for(int i=0; i<numlines; i++)
    {
        double isoval = minval + (maxval-minval) * double(i)/double(numlines);
        extractIsoline(theta, isoval, minval, maxval, isotraces);
    }
    std::cout << "Extracted " << isotraces.size() << " isolines" << std::endl;
}

void CoverMesh::computeFunc(double globalScale)
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
    Eigen::VectorXd scales = globalScale * s;
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
        double eigenVal = eigenVec.transpose() * (Lmat * eigenVec);
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


void CoverMesh::createVisualization(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXd &edgePts, Eigen::MatrixXd &edgeVecs, Eigen::MatrixXi &edgeSegs, Eigen::MatrixXd &colors, 
    Eigen::MatrixXd &cutPts1, Eigen::MatrixXd &cutPts2, Eigen::MatrixXd &cutColors)
{
    int splitFace = data_.splitMesh->nFaces();
    int origverts = parent_.fs->nVerts();
    int origfaces = parent_.fs->nFaces();
    V = data_.splitMesh->data().V;
    F = data_.splitMesh->data().F;
    
    edgePts.resize(splitFace, 3);
    edgeVecs.resize(splitFace, 3);
    edgeVecs.setZero();
    edgeSegs.resize(splitFace, 2);
    colors.resize(splitFace , 3);

    for(int c=0; c<ncovers_; c++)
    {
        for (int i = 0; i < origfaces; i++)
        {
            Eigen::Vector3d centroid;
            centroid.setZero();
            for (int j = 0; j < 3; j++)
                centroid += renderScale_ * parent_.fs->data().V.row(parent_.fs->data().F(i, j));
            centroid /= 3.0;
            centroid += data_.splitOffsets[c];

            edgePts.row(c*origfaces + i) = centroid;
            edgeVecs.row(c*origfaces + i) = parent_.fs->data().Bs[i] * fs->v(c*origfaces + i, 0);
            edgeSegs(c*origfaces + i, 0) = 2 * (c*origfaces + i);
            edgeSegs(c*origfaces + i, 1) = 2 * (c*origfaces + i) + 1;
            colors.row(c*origfaces + i) = Eigen::Vector3d(0,0,0).transpose();
        }
    }

    int ncutedges = data_.splitMeshCuts.size();
    int nsliceedges = slicedEdges.size();
    cutPts1.resize(ncutedges+nsliceedges, 3);
    cutPts2.resize(ncutedges+nsliceedges, 3);
    cutColors.resize(ncutedges + nsliceedges, 3);
    for (int i = 0; i < ncutedges; i++)
    {
        int edgeid = data_.splitMeshCuts[i];
        int v0 = data_.splitMesh->data().edgeVerts(edgeid, 0);
        int v1 = data_.splitMesh->data().edgeVerts(edgeid, 1);
        Eigen::Vector3d n(0, 0, 0);
        int f0 = data_.splitMesh->data().E(edgeid, 0);
        int f1 = data_.splitMesh->data().E(edgeid, 1);
        if (f0 != -1)
            n += data_.splitMesh->faceNormal(f0);
        if (f1 != -1)
            n += data_.splitMesh->faceNormal(f1);
        Eigen::Vector3d offset = 0.0001*n / n.norm();
        cutPts1.row(i) = data_.splitMesh->data().V.row(v0) + offset.transpose();
        cutPts2.row(i) = data_.splitMesh->data().V.row(v1) + offset.transpose();
        cutColors.row(i) = Eigen::RowVector3d(0.9, .1, .9);
    }
    for (int i = 0; i < nsliceedges; i++)
    {
        int v0 = data_.splitMesh->data().edgeVerts(slicedEdges[i], 0);
        int v1 = data_.splitMesh->data().edgeVerts(slicedEdges[i], 1);
        Eigen::Vector3d n(0, 0, 0);
        int f0 = data_.splitMesh->data().E(slicedEdges[i], 0);
        int f1 = data_.splitMesh->data().E(slicedEdges[i], 1);
        if (f0 != -1)
            n += data_.splitMesh->faceNormal(f0);
        if (f1 != -1)
            n += data_.splitMesh->faceNormal(f1);
        Eigen::Vector3d offset = 0.0001*n / n.norm();
        cutPts1.row(i + ncutedges) = data_.splitMesh->data().V.row(v0) + offset.transpose();
        cutPts2.row(i + ncutedges) = data_.splitMesh->data().V.row(v1) + offset.transpose();
        cutColors.row(i + ncutedges) = Eigen::RowVector3d(0.1, .9, .9);
    }
}

int CoverMesh::visMeshToCoverMesh(int vertid)
{
    return data_.splitToCoverVerts[vertid];
}


void CoverMesh::initializeS(double reg)
{
    theta.resize(fs->nVerts());
    theta.setZero();
    Eigen::VectorXi thetacnt(fs->nVerts());
    thetacnt.setZero();

    // some sanity checks
    Eigen::VectorXi B;
    if(!igl::is_vertex_manifold(fs->data().F, B))
        std::cout << "ERROR: cover mesh not vertex-manifold!" << std::endl;
    if(!igl::is_edge_manifold(fs->data().F))
        std::cout << "ERROR: cover mesh not edge-manifold!" << std::endl;
    // compute how to cut the mesh
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

    // update data structures for visualizing the new cuts
    std::set<std::pair<int, int> > cutsegs;
    for (auto &it : cuts)
    {
        for (int i = 0; i < it.size() - 1; i++)
        {
            int v0 = it[i];
            int v1 = it[i + 1];
            if (v0 > v1)
                std::swap(v0, v1);
            cutsegs.insert(std::pair<int, int>(v0, v1));
        }
    }

    slicedEdges.clear();
    int splitedges = data_.splitMesh->nEdges();
    for (int i = 0; i < splitedges; i++)
    {
        int v0 = data_.splitMesh->data().edgeVerts(i, 0);
        int v1 = data_.splitMesh->data().edgeVerts(i, 1);
        int mappedv0 = data_.splitToCoverVerts[v0];
        int mappedv1 = data_.splitToCoverVerts[v1];
        if (mappedv0 > mappedv1)
            std::swap(mappedv0, mappedv1);
        if (cutsegs.count(std::pair<int, int>(mappedv0, mappedv1)))
        {
            slicedEdges.push_back(i);
        }
    }
    
    // cut the mesh
    Eigen::MatrixXd cutV;
    Eigen::MatrixXi cutF;    
    cutMesh(fs->data().V, fs->data().F, cuts, cutV, cutF);
    igl::writeOBJ("debug_augmented.obj",fs->data().V, fs->data().F);
    igl::writeOBJ("debug_cut.obj",cutV,cutF);


    // separate cut mesh into connected components
    Eigen::VectorXi components;
    
    igl::facet_components(cutF, components);
    int ncomponents = 0;
    for(int i=0; i<components.size(); i++)    
        ncomponents = std::max(ncomponents, components[i]);
    ncomponents++;
    std::vector<int> componentsizes;
    for(int i=0; i<ncomponents; i++)
        componentsizes.push_back(0);
    for(int i=0; i<components.size(); i++)
        componentsizes[components[i]]++;    
    std::cout << "Covering mesh has " << ncomponents << " connected components" << std::endl;
    // loop over the connected components
    for(int component = 0; component < ncomponents; component++)
    {
        std::cout << "Component " << component << ": " << componentsizes[component] << " faces" << std::endl;
        // faces for just this connected component
        Eigen::VectorXi compFacesToGlobal(componentsizes[component]);
        Eigen::MatrixXi compF(componentsizes[component], 3);
        int idx=0;
        for(int i=0; i<components.size(); i++)
        {
            if(components[i] == component)
            {
                compFacesToGlobal[idx] = i;
                compF.row(idx) = cutF.row(i);
                idx++;
            }
        }
        
        Eigen::MatrixXd prunedV;
        Eigen::MatrixXi prunedF;
        Eigen::VectorXi I;
        igl::remove_unreferenced(cutV, compF, prunedV, prunedF, I);
        // connected component surface
        Surface surf(prunedV, prunedF);
        
        std::cout << "Built connected component surface" << std::endl;

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

        // vertex mass matrix
        Eigen::SparseMatrix<double> Mvert;
        igl::massmatrix(surf.data().V, surf.data().F, igl::MassMatrixType::MASSMATRIX_TYPE_BARYCENTRIC, Mvert);

        // face mass matrix
        Eigen::VectorXd faceAreas;
        igl::doublearea(surf.data().V, surf.data().F, faceAreas);
        faceAreas *= 0.5;
        
        std::cout << "Built mass matrices" << std::endl;

        std::vector<Eigen::Triplet<double> > AhalfCoeffs;

        int nverts = surf.nVerts();

        for(int i=0; i<nedges; i++)
        {
            int f0 = surf.data().E(i,0);
            int f1 = surf.data().E(i,1);
            int vert0 = surf.data().edgeVerts(i, 0);
            int vert1 = surf.data().edgeVerts(i, 1);
            AhalfCoeffs.push_back(Eigen::Triplet<double>(i, vert0, -1));
            AhalfCoeffs.push_back(Eigen::Triplet<double>(i, vert1, 1));
            double denom = 0.0;
            if (f0 != -1)
                denom += 1.0;
            if (f1 != -1)
                denom += 1.0;

            Eigen::Vector3d v0 = surf.data().V.row(vert0).transpose();
            Eigen::Vector3d v1 = surf.data().V.row(vert1).transpose();
            Eigen::Vector3d edgeVec = v1-v0;

            if (f0 != -1)
            {
                Eigen::Vector3d scaledvec0 = surf.data().Bs[f0] * surf.data().Js.block<2, 2>(2 * f0, 0) * fs->v(compFacesToGlobal[f0], 0);
                AhalfCoeffs.push_back(Eigen::Triplet<double>(i, nverts + f0, scaledvec0.dot(edgeVec) / denom));
            }
            if (f1 != -1)
            {
                Eigen::Vector3d scaledvec1 = surf.data().Bs[f1] * surf.data().Js.block<2, 2>(2 * f1, 0) * fs->v(compFacesToGlobal[f1], 0);
                AhalfCoeffs.push_back(Eigen::Triplet<double>(i, nverts + f1, scaledvec1.dot(edgeVec) / denom));
            }
        }

        Eigen::SparseMatrix<double> Ahalf(nedges, nverts + nfaces);
        Ahalf.setFromTriplets(AhalfCoeffs.begin(), AhalfCoeffs.end());
        
        std::vector<Eigen::Triplet<double> > edgeMetricMatrixCoeffs;
        for(int i = 0; i < nedges; i++)
        {
            edgeMetricMatrixCoeffs.push_back(Eigen::Triplet<double>(i, i, edgeMetric[i]));
        }
        
        Eigen::SparseMatrix<double> edgeMetricMatrix(nedges, nedges);
        edgeMetricMatrix.setFromTriplets(edgeMetricMatrixCoeffs.begin(), edgeMetricMatrixCoeffs.end());
        
        Eigen::SparseMatrix<double> A = Ahalf.transpose() * edgeMetricMatrix * Ahalf;

        std::vector<Eigen::Triplet<double> > DfaceCoeffs;
        for (int i = 0; i < nedges; i++)
        {
            int f0 = surf.data().E(i,0);
            int f1 = surf.data().E(i,1);
            if (f0 != -1 && f1 != -1)
            {
                DfaceCoeffs.push_back(Eigen::Triplet<double>(i, nverts + f0, -1.0));
                DfaceCoeffs.push_back(Eigen::Triplet<double>(i, nverts + f1, 1.0));
            }
        }
        Eigen::SparseMatrix<double> Dface(nedges, nverts + nfaces);
        Dface.setFromTriplets(DfaceCoeffs.begin(), DfaceCoeffs.end());
        
        std::vector<Eigen::Triplet<double> > inverseEdgeMetricCoeffs;
        for (int i = 0; i < nedges; i++)
        {
            inverseEdgeMetricCoeffs.push_back(Eigen::Triplet<double>(i, i, 1.0 / edgeMetric[i]));
        }
        Eigen::SparseMatrix<double> inverseEdgeMetric(nedges, nedges);
        inverseEdgeMetric.setFromTriplets(inverseEdgeMetricCoeffs.begin(), inverseEdgeMetricCoeffs.end());

        Eigen::SparseMatrix<double> Lface = Dface.transpose() * inverseEdgeMetric * Dface;

        Eigen::SparseMatrix<double> Areg = A + reg * Lface;

        std::vector<Eigen::Triplet<double> > Bcoeffs;
        for (int i = 0; i < nverts; i++)
        {
            Bcoeffs.push_back(Eigen::Triplet<double>(i, i, Mvert.coeff(i, i)));
        }
        for (int i = 0; i < nfaces; i++)
        {
            Bcoeffs.push_back(Eigen::Triplet<double>(nverts + i, nverts + i, faceAreas[i]));
        }

        Eigen::SparseMatrix<double> B(nverts+nfaces, nverts+nfaces);
        B.setFromTriplets(Bcoeffs.begin(), Bcoeffs.end());
        
        std::cout << "Factoring A" << std::endl;

        // inverse power iteration
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver(Areg);
        Eigen::VectorXd x(nverts + nfaces);
        x.setRandom();
        double xnorm = x.transpose() * (B * x);
        x /= sqrt(xnorm);
        Eigen::VectorXd ones(nverts + nfaces);
        for (int i = 0; i < nverts; i++)
            ones[i] = 1.0;
        for (int i = 0; i < nfaces; i++)
            ones[nverts + i] = 0.0;
        ones.normalize();

        std::cout << "Starting power iteration" << std::endl;

        for (int i = 0; i < 1000; i++)
        {
            Eigen::VectorXd newx = solver.solve(B*x);
            // project out ones component
            newx = newx - ones.dot(newx)*ones;
            x = newx / sqrt(newx.transpose() * A * newx);
        }

        std::cout << "Rayleigh quotient: " << (x.transpose() * (A * x)) / (x.transpose() * (B * x)) << std::endl;

        Eigen::VectorXd componentS(nfaces);
        for (int i = 0; i < nfaces; i++)
            componentS[i] = x[nverts + i];

        double maxS = 0;
        for(int i=0; i<nfaces; i++)
        {
            if ( fabs(componentS[i]) > maxS ) 
            {
                maxS = fabs(componentS[i]);
            }
        }

        double s_scale = 3.1415 / maxS;

        
        // map component s to the global s vector
        for(int i=0; i<nfaces; i++)
        {
            componentS[i] *= s_scale;
            s[compFacesToGlobal[i]] = componentS[i] ;
        }

        // compute theta
        Eigen::VectorXd componentTheta(nverts);
        for (int i = 0; i < nverts; i++)
            componentTheta[i] = x[i];
        
        for (int i = 0; i < nfaces; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                int cutv = surf.data().F(i, j);
                int coverv = fs->data().F(compFacesToGlobal[i], j);
                theta[coverv] += componentTheta[cutv];
                thetacnt[coverv]++;
            }
        }
    }
    
    // finally, we have the global theta initialization
    for (int i = 0; i < fs->nVerts(); i++)
    {
        if (thetacnt[i])
            theta[i] /= thetacnt[i];
        else
            theta[i] = 0;
    }        
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

void CoverMesh::initializeSplitMesh(const Eigen::VectorXi &oldToNewVertMap)
{
    data_.splitToCoverVerts = oldToNewVertMap;
    int facespercover = fs->nFaces() / ncovers_;
    int rows = 2;
    int meshesperrow = ncovers_ / rows + (ncovers_ % rows == 0 ? 0 : 1);
    data_.splitOffsets.clear();
    for (int i = 0; i < ncovers_; i++)
    {
        int row = i / meshesperrow;
        int col = i%meshesperrow;
        double dy = (-1.1 * row + (1.1) * (rows - row - 1)) / double(rows);
        double dx = (1.1 * col + (-1.1) * (meshesperrow - col - 1)) / double(meshesperrow);
        data_.splitOffsets.push_back(Eigen::Vector3d(dx, dy, 0.0));
    }

    int origverts = parent_.fs->nVerts();
    int origfaces = parent_.fs->nFaces();
    int newverts = ncovers_*origverts;
    int newfaces = ncovers_*origfaces;
    Eigen::MatrixXd V(newverts, 3);
    Eigen::MatrixXi F(newfaces, 3);
    renderScale_ = 1.0 / std::max(rows, meshesperrow);
    for (int i = 0; i < ncovers_; i++)
    {
        for (int j = 0; j < origverts; j++)
        {
            V.row(i*origverts + j) = data_.splitOffsets[i].transpose() + renderScale_ * parent_.fs->data().V.row(j);
        }
        for (int j = 0; j < origfaces; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                F(i*origfaces + j, k) = i*origverts + parent_.fs->data().F(j, k);
            }
        }
    }
    data_.splitMesh = new Surface(V, F);

    data_.coverToSplitVerts.clear();
    for (int i = 0; i < oldToNewVertMap.rows(); i++)
    {
        data_.coverToSplitVerts[oldToNewVertMap[i]].push_back(i);
    }

    data_.splitMeshCuts.clear();
    for (int i = 0; i < newfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int edge = fs->data().faceEdges(i, j);
            int f0 = fs->data().E(edge, 0);
            int f1 = fs->data().E(edge, 1);
            if (f0 == -1 || f1 == -1)
                continue;
            int face0copy = f0 / origfaces;
            int face1copy = f1 / origfaces;
            if (face0copy != face1copy)
            {
                data_.splitMeshCuts.push_back(data_.splitMesh->data().faceEdges(i, j));
            }
        }
    }    
}

const Surface &CoverMesh::splitMesh() const
{
    return *data_.splitMesh;
}

void CoverMesh::drawTraceOnSplitMesh(const Trace &trace, Eigen::MatrixXd &pathStarts, Eigen::MatrixXd &pathEnds) const
{
    int nsegs = trace.segs.size();
    pathStarts.resize(nsegs, 3);
    pathEnds.resize(nsegs, 3);
    for(int i=0; i<nsegs; i++)
    {
        Eigen::Vector3d offset = 0.0001 * data_.splitMesh->faceNormal(trace.segs[i].face);        
        int v0 = data_.splitMesh->data().F(trace.segs[i].face, (trace.segs[i].side[0]+1)%3);
        int v1 = data_.splitMesh->data().F(trace.segs[i].face, (trace.segs[i].side[0]+2)%3);
        Eigen::Vector3d pos = (1.0 - trace.segs[i].bary[0])*data_.splitMesh->data().V.row(v0).transpose() + trace.segs[i].bary[0] * data_.splitMesh->data().V.row(v1).transpose();
        pathStarts.row(i) = pos.transpose() + offset.transpose();
                
        v0 = data_.splitMesh->data().F(trace.segs[i].face, (trace.segs[i].side[1]+1)%3);
        v1 = data_.splitMesh->data().F(trace.segs[i].face, (trace.segs[i].side[1]+2)%3);
        pos = (1.0 - trace.segs[i].bary[1])*data_.splitMesh->data().V.row(v0).transpose() + trace.segs[i].bary[1] * data_.splitMesh->data().V.row(v1).transpose();
        pathEnds.row(i) = pos.transpose() + offset.transpose();
    }
}



