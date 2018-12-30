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
#include "FieldIntegration.h"
#include "igl/cotmatrix.h"
#include "CoMISoWrapper.h"

# define M_PI           3.14159265358979323846

using namespace std;

CoverMesh::CoverMesh(const Surface &originalSurf, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXi &oldToNewVertMap, const Eigen::MatrixXd &field, int ncovers)
{
    originalSurf_ = new Surface(originalSurf);
    fs = new FieldSurface(V, F, 1);
    int nfaces = F.rows();
    ncovers_ = ncovers;
    for (int i = 0; i < nfaces; i++)
    {
        fs->vectorFields.segment<2>(fs->vidx(i, 0)) = field.row(i).transpose();
    }

    theta.resize(fs->nVerts());
    theta.setZero();

    scales.resize(fs->nFaces());
    scales.setZero();
    
    renderScale_ = 1.0;

    initializeSplitMesh(oldToNewVertMap);
}

CoverMesh::~CoverMesh()
{
    delete fs;
    if (data_.splitMesh)
        delete data_.splitMesh;
    delete originalSurf_;
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
        visited[i] = fs->isFaceDeleted(i);
        

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

void CoverMesh::createVisualization(Eigen::MatrixXd &V, Eigen::MatrixXi &F, 
    Eigen::MatrixXd &edgePts, Eigen::MatrixXi &edgeSegs, Eigen::MatrixXd &colors, 
    Eigen::MatrixXd &cutPts1, Eigen::MatrixXd &cutPts2, Eigen::MatrixXd &cutColors,
    bool hideVectors, double vectorLength)
{
    int splitFace = data_.splitMesh->nFaces();
    int origverts = originalSurf_->nVerts();
    int origfaces = originalSurf_->nFaces();
    V = data_.splitMesh->data().V;
    F = data_.splitMesh->data().F;
    
    if (hideVectors)
    {
        edgePts.resize(0, 3);
        edgeSegs.resize(0, 3);
        colors.resize(0, 3);
    }
    else
    {
        edgePts.resize(2 * splitFace, 3);
        edgePts.setZero();
        edgeSegs.resize(splitFace, 2);
        colors.resize(splitFace, 3);

        for (int c = 0; c < ncovers_; c++)
        {
            for (int i = 0; i < origfaces; i++)
            {
                Eigen::Vector3d centroid;
                centroid.setZero();
                for (int j = 0; j < 3; j++)
                    centroid += renderScale_ * originalSurf_->data().V.row(originalSurf_->data().F(i, j));
                centroid /= 3.0;
                centroid += data_.splitOffsets[c];

                edgePts.row(2 * c*origfaces + 2 * i) = centroid.transpose();
                Eigen::Vector3d vec = originalSurf_->data().Bs[i] * fs->v(c*origfaces + i, 0);
                vec *= vectorLength * renderScale_ * fs->data().averageEdgeLength / vec.norm() * sqrt(3.0) / 6.0 * 0.75;
                edgePts.row(2 * c * origfaces + 2 * i + 1) = (centroid + vec).transpose();
                edgeSegs(c*origfaces + i, 0) = 2 * (c*origfaces + i);
                edgeSegs(c*origfaces + i, 1) = 2 * (c*origfaces + i) + 1;
                colors.row(c*origfaces + i) = Eigen::Vector3d(0, 0, 0).transpose();
            }
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

void CoverMesh::roundAntipodalCovers(int numISOLines)
{
    // create a mesh without deleted faces
    int undeletedFaces = fs->numUndeletedFaces();
    Eigen::MatrixXi undelF(undeletedFaces, 3);
    Eigen::VectorXi undelFaceMap(undeletedFaces);
    int fid=0;
    int globalfaces = fs->nFaces();
    for(int i=0; i<globalfaces; i++)
    {
        if(!fs->isFaceDeleted(i))
        {
            undelFaceMap[fid] = i;
            undelF.row(fid) = fs->data().F.row(i);
            fid++;
        }        
    }

    Eigen::MatrixXd prunedV;
    Eigen::MatrixXi prunedF;
    Eigen::VectorXi I;
    igl::remove_unreferenced(fs->data().V, undelF, prunedV, prunedF, I);

    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(prunedV, prunedF, L);
    
    int nfields = ncovers_ / 2;
    int origverts = originalSurf_->nVerts();
    int newverts = prunedV.rows();
    
    double phase = 2.0 * M_PI / double(numISOLines);
    double offset = M_PI / numISOLines;

    Eigen::VectorXd result(newverts + nfields * origverts);
    result.setZero();

    std::vector<Eigen::Triplet<double> > Ccoeffs;
    int row = 0;
    for (int i = 0; i < origverts; i++)
    {
        for (int j = 0; j < nfields; j++)
        {
            int v1 = j * origverts + i;
            int v2 = (j + nfields)*origverts + i;
            int covv1 = data_.splitToCoverVerts[v1];
            int covv2 = data_.splitToCoverVerts[v2];
            if (I[covv1] != -1 && I[covv2] != -1)
            {
                double thetadiff = theta[covv1] + theta[covv2];
                Ccoeffs.push_back(Eigen::Triplet<double>(row, I[covv1], 1.0));
                Ccoeffs.push_back(Eigen::Triplet<double>(row, I[covv2], 1.0));
                Ccoeffs.push_back(Eigen::Triplet<double>(row, newverts + j * origverts + i, phase));
                Ccoeffs.push_back(Eigen::Triplet<double>(row, newverts + nfields * origverts, offset + thetadiff));
                row++;
            }
        }
    }
    Eigen::SparseMatrix<double> C(row, newverts + nfields * origverts + 1);
    C.setFromTriplets(Ccoeffs.begin(), Ccoeffs.end());

    std::vector<Eigen::Triplet<double> > Acoeffs;
    for (int k = 0; k < L.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it)
        {
            Acoeffs.push_back(Eigen::Triplet<double>(it.row(), it.col(), -it.value()));
        }
    }
    for (int i = 0; i < newverts + nfields*origverts; i++)
        Acoeffs.push_back(Eigen::Triplet<double>(i, i, 1e-4));
    Eigen::SparseMatrix<double> A(newverts + nfields*origverts, newverts + nfields*origverts);
    A.setFromTriplets(Acoeffs.begin(), Acoeffs.end());

    Eigen::VectorXd rhs(newverts + nfields*origverts);
    rhs.setZero();

    Eigen::VectorXi toRound(nfields*origverts);
    for (int i = 0; i < nfields*origverts; i++)
        toRound[i] = newverts + i;

    ComisoWrapper(C, A, result, rhs, toRound, 1e-6);
    std::cout << "Residual: " << (A*result - rhs).norm() << std::endl;
    Eigen::VectorXd ctest(newverts + nfields * origverts + 1);
    ctest.segment(0, newverts + nfields * origverts) = result;
    ctest[newverts + nfields * origverts] = 1.0;
    std::cout << "Constraint residual: " << (C*ctest).norm() << std::endl;
    for (int i = 0; i < fs->nVerts(); i++)
    {
        if (I[i] != -1)
        {
            theta[i] += result[I[i]];
            theta[i] = std::remainder(theta[i], 2.0*M_PI);
        }
    }
}

void CoverMesh::integrateField(LocalFieldIntegration *lmethod, GlobalFieldIntegration *gmethod, double globalScale)
{
    int globalverts = fs->nVerts();
    theta.resize(globalverts);
    theta.setZero();
    scales.resize(fs->nFaces());
    scales.setZero();

    // create a mesh without deleted faces
    int undeletedFaces = fs->numUndeletedFaces();
    Eigen::MatrixXi undelF(undeletedFaces, 3);
    Eigen::VectorXi undelFaceMap(undeletedFaces);
    int fid=0;
    int globalfaces = fs->nFaces();
    for(int i=0; i<globalfaces; i++)
    {
        if(!fs->isFaceDeleted(i))
        {
            undelFaceMap[fid] = i;
            undelF.row(fid) = fs->data().F.row(i);
            fid++;
        }        
    }
    
    // separate cut mesh into connected components
    Eigen::VectorXi components;    
    igl::facet_components(undelF, components);
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
        Eigen::MatrixXd compField(componentsizes[component], 2);
        int idx=0;
        for(int i=0; i<components.size(); i++)
        {
            if(components[i] == component)
            {
                compFacesToGlobal[idx] = undelFaceMap[i];
                compF.row(idx) = undelF.row(i);
                Eigen::Vector2d vec = fs->v(undelFaceMap[i], 0); 
                double vecnorm = (fs->data().Bs[undelFaceMap[i]] * vec).norm();
                compField.row(idx) = vec.transpose()/vecnorm;
                idx++;
            }
        }
        
        Eigen::MatrixXd prunedV;
        Eigen::MatrixXi prunedF;
        Eigen::VectorXi I;
        igl::remove_unreferenced(fs->data().V, compF, prunedV, prunedF, I);
        // connected component surface
        Surface surf(prunedV, prunedF);
        
        std::cout << "Built connected component surface" << std::endl;

        // component theta and s
        Eigen::VectorXd compS;
        Eigen::VectorXd compTheta;
        lmethod->locallyIntegrateOneComponent(surf, compField, compS);
        
        double maxS = 0;
        for(int i=0; i<compS.size(); i++)
        {
            if ( fabs(compS[i]) > maxS ) 
            {
                maxS = fabs(compS[i]);
            }
        }

        double s_scale = 3.1415 / fs->data().averageEdgeLength / maxS;
        compS *= globalScale * s_scale;

        for (int i = 0; i < componentsizes[component]; i++)
        {
            scales[compFacesToGlobal[i]] = compS[i];
        }

        gmethod->globallyIntegrateOneComponent(surf, compField, compS, compTheta);
        
        
        // map component theta to the global theta vector
        for (int i = 0; i < globalverts; i++)
        {
            if (I[i] != -1)
                theta[i] = compTheta[I[i]];            
        }        
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

    int origverts = originalSurf_->nVerts();
    int origfaces = originalSurf_->nFaces();
    int newverts = ncovers_*origverts;
    int newfaces = ncovers_*origfaces;
    Eigen::MatrixXd V(newverts, 3);
    Eigen::MatrixXi F(newfaces, 3);
    renderScale_ = 1.0 / std::max(rows, meshesperrow);
    for (int i = 0; i < ncovers_; i++)
    {
        for (int j = 0; j < origverts; j++)
        {
            V.row(i*origverts + j) = data_.splitOffsets[i].transpose() + renderScale_ * originalSurf_->data().V.row(j);
        }
        for (int j = 0; j < origfaces; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                F(i*origfaces + j, k) = i*origverts + originalSurf_->data().F(j, k);
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
            if(edge == -1)
                continue;
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



