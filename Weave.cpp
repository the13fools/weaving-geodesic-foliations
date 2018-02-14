#include "Weave.h"
#include <math.h>
#include <igl/read_triangle_mesh.h>
#include <map>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/Eigenvalues> 
#include "Colors.h"
#include <deque>
#include <queue>
#include <algorithm>
#include <set>
#include <igl/remove_unreferenced.h>
#include <igl/writeOBJ.h>

typedef Eigen::Triplet<double> triplet;
# define M_PI           3.14159265358979323846

Weave::Weave(const std::string &objname, int m)
{
    Eigen::MatrixXd Vtmp;
    if (!igl::read_triangle_mesh(objname, Vtmp, F))
    {
        std::cerr << "Couldn't load mesh " << objname << std::endl;
        exit(-1);
    }
    if (Vtmp.cols() < 3)
    {
        std::cerr << "Mesh must 3D" << std::endl;
        exit(-1);
    }
    V.resize(Vtmp.rows(), 3);
    //wtf
    for (int i = 0; i < 3; i++)
    {
        V.col(i) = Vtmp.col(i);
    }

    centerAndScale();
    buildConnectivityStructures();
    buildGeometricStructures();

    // initialize vector fields
    nFields_ = m;
    vectorFields.resize(5*F.rows()*m);
    vectorFields.setZero();
    vectorFields.segment(0, 2 * F.rows()*m).setRandom();
    normalizeFields();

    // initialize permutation matrices
    int nedges = nEdges();
    Ps.resize(nedges);
    for (int i = 0; i < nedges; i++)
    {
        Ps[i].resize(m, m);
        Ps[i].setIdentity();
    }
    augmented = false;
}

Weave::~Weave()
{    
}

void Weave::centerAndScale()
{
    Eigen::Vector3d centroid(0, 0, 0);
    for (int i = 0; i < V.rows(); i++)
        centroid += V.row(i);
    centroid /= V.rows();

    double maxdist = 0;
    for (int i = 0; i < V.rows(); i++)
    {
        maxdist = std::max(maxdist, (V.row(i).transpose() - centroid).norm());
    }
    for (int i = 0; i < V.rows(); i++)
    {
        Eigen::Vector3d newpos = V.row(i).transpose() - centroid;
        V.row(i) = newpos / maxdist;
    }
}

void Weave::buildConnectivityStructures()
{
    std::map<std::pair<int, int>, Eigen::Vector2i > edgemap;
    int nfaces = nFaces();
    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int nextj = (j + 1) % 3;
            int v1 = F(i, j);
            int v2 = F(i, nextj);
            int idx = 0;
            if (v1 > v2)
            {
                idx = 1;
                std::swap(v1, v2);
            }
            std::pair<int, int> p(v1,v2);
            std::map<std::pair<int, int>, Eigen::Vector2i >::iterator it = edgemap.find(p);
            if(it == edgemap.end())
            {
                Eigen::Vector2i edge(-1,-1);
                edge[idx] = i;
                edgemap[p] = edge;
            }
            else
            {
                edgemap[p][idx] = i;
            }
        }
    }

    E.resize(edgemap.size(), 2);
    faceEdges.resize(nfaces, 3);
    edgeVerts.resize(edgemap.size(), 2);
    int idx = 0;
    for (std::map<std::pair<int, int>, Eigen::Vector2i >::iterator it = edgemap.begin(); it != edgemap.end(); ++it)
    {
        E(idx, 0) = it->second[0];
        E(idx, 1) = it->second[1];
        edgeVerts(idx, 0) = it->first.first;
        edgeVerts(idx, 1) = it->first.second;
        idx++;
    }
    faceEdges.resize(nfaces, 3);
    faceEdges.setConstant(-1);
    faceNeighbors.resize(nfaces, 3);
    faceNeighbors.setConstant(-1);
    faceWings.resize(nfaces, 3);
    faceWings.setConstant(-1);
    int nedges = nEdges();
    for (int edge = 0; edge < nedges; edge++)
    {
        for(int side = 0; side<2; side++)
        {
            if(E(edge,side) == -1)
                 continue;
            for(int j=0; j<3; j++)
                 if(F(E(edge,side), j) != edgeVerts(edge,0) && F(E(edge,side), j) != edgeVerts(edge,1))
                     faceEdges(E(edge, side), j) = edge;           
        }
        if(E(edge,0) == -1 || E(edge,1) == -1)
            continue;
        Eigen::Vector3i face1 = F.row(E(edge, 0));
        Eigen::Vector3i face2 = F.row(E(edge, 1));
        int idx1 = -1;
        for (int i = 0; i < 3; i++)
        {
            bool ok = true;
            for (int j = 0; j < 3; j++)
            {
                if (face2[j] == face1[i])
                    ok = false;
            }
            if (ok)
                idx1 = i;
        }
        int idx2 = -1;
        for (int i = 0; i < 3; i++)
        {
            bool ok = true;
            for (int j = 0; j < 3; j++)
            {
                if (face1[j] == face2[i])
                    ok = false;
            }
            if (ok)
                idx2 = i;
        }
        faceNeighbors(E(edge, 0), idx1) = E(edge, 1);
        faceNeighbors(E(edge, 1), idx2) = E(edge, 0);
        faceWings(E(edge, 0), idx1) = face2[idx2];
        faceWings(E(edge, 1), idx2) = face1[idx1];
    }

    vertEdges.resize(V.rows());
    for(int i=0; i<nedges; i++)
    {
        vertEdges[edgeVerts(i,0)].push_back(i);
        vertEdges[edgeVerts(i,1)].push_back(i);
    }
}

Eigen::Vector3d Weave::faceNormal(int face) const
{
    Eigen::Vector3d e1 = (V.row(F(face, 1)) - V.row(F(face, 0)));
    Eigen::Vector3d e2 = (V.row(F(face, 2)) - V.row(F(face, 0)));
    Eigen::Vector3d result = e1.cross(e2);
    result /= result.norm();
    return result;
}

void Weave::buildGeometricStructures()
{
    // compute barycentric matrices and Js
    int nfaces = nFaces();
    Bs.resize(nfaces);
    Js.resize(2 * nfaces, 2);

    averageEdgeLength = 0;

    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector3d v0 = V.row(F(i, 0));
        Eigen::Vector3d v1 = V.row(F(i, 1));
        Eigen::Vector3d v2 = V.row(F(i, 2));
        Bs[i].col(0) = v1 - v0;
        Bs[i].col(1) = v2 - v0;

        averageEdgeLength += (v2 - v1).norm();
        averageEdgeLength += (v2 - v0).norm();
        averageEdgeLength += (v1 - v0).norm();

        Eigen::Vector3d n = (v1 - v0).cross(v2 - v0);
        n /= n.norm();

        Eigen::Matrix2d BTB = Bs[i].transpose() * Bs[i];
        Eigen::Matrix<double, 3, 2> ncrossB;
        ncrossB.col(0) = n.cross(v1 - v0);
        ncrossB.col(1) = n.cross(v2 - v0);
        Js.block<2, 2>(2 * i, 0) = BTB.inverse() * Bs[i].transpose() * ncrossB;
    }

    averageEdgeLength /= (3.0 * nfaces);

    // compute cDiffs and transition matrices
    int nedges = nEdges();
    cDiffs.resize(2 * nedges, 2);
    Ts.resize(2 * nedges, 4);
    for (int edgeidx = 0; edgeidx < nedges; edgeidx++)
    {        
        //collect neighboring face indices
        int face1 = E(edgeidx, 0);
        int face2 = E(edgeidx, 1);

        if(face1 == -1 || face2 == -1)
            continue;

        int v1 = edgeVerts(edgeidx, 0);
        int v2 = edgeVerts(edgeidx, 1);

        Eigen::Vector3d n1 = faceNormal(face1);
        Eigen::Vector3d n2 = faceNormal(face2);
        // collect: (1) the midpoint of the common edge, (2) unit vector in direction of common edge,
        // (3) the face normals, (4) the centroid of neighboring faces

        Eigen::Vector3d midpt = 0.5 * (V.row(v1) + V.row(v2));
        Eigen::Vector3d commone = V.row(v2) - V.row(v1);
        commone /= commone.norm();
        Eigen::Vector3d centroids[2];
        centroids[0].setZero();
        centroids[1].setZero();

        for (int i = 0; i < 3; i++)
        {
            centroids[0] += V.row(F(face1, i));
            centroids[1] += V.row(F(face2, i));
        }

        centroids[0] /= 3.0;
        centroids[1] /= 3.0;

        //rotate each centroid into the plane of the opposite triangle and compute ci minus c

        Eigen::Vector3d t1 = n1.cross(commone);
        Eigen::Vector3d t2 = n2.cross(commone);
        Eigen::Vector3d diff2 = centroids[1] - midpt;
        double alpha = commone.dot(diff2);
        double beta = t2.dot(diff2);

        Eigen::Matrix2d BTB1 = Bs[face1].transpose() * Bs[face1];
        Eigen::Matrix2d BTB2 = Bs[face2].transpose() * Bs[face2];

        cDiffs.row(2 * edgeidx) = BTB1.inverse() * Bs[face1].transpose() * (midpt + alpha * commone + beta * t1 - centroids[0]);
        Eigen::Vector3d diff1 = centroids[0] - midpt;
        alpha = commone.dot(diff1);
        beta = t1.dot(diff1);
        cDiffs.row(2 * edgeidx + 1) = BTB2.inverse() * Bs[face2].transpose() * (midpt + alpha*commone + beta * t2 - centroids[1]);

        Eigen::Vector3d e1 = V.row(F(face1, 1)) - V.row(F(face1, 0));
        Eigen::Vector3d e2 = V.row(F(face1, 2)) - V.row(F(face1, 0));

        double alpha1 = commone.dot(e1);
        double beta1 = t1.dot(e1);
        Eigen::Vector3d newe1 = alpha1*commone + beta1 * t2;
        Ts.block<2, 1>(2 * edgeidx, 0) = BTB2.inverse() * Bs[face2].transpose() * newe1;

        double alpha2 = commone.dot(e2);
        double beta2 = t1.dot(e2);
        Eigen::Vector3d newe2 = alpha2*commone + beta2*t2;
        Ts.block<2, 1>(2 * edgeidx, 1) = BTB2.inverse() * Bs[face2].transpose() * newe2;

        e1 = V.row(F(face2, 1)) - V.row(F(face2, 0));
        e2 = V.row(F(face2, 2)) - V.row(F(face2, 0));

        alpha1 = commone.dot(e1);
        beta1 = t2.dot(e1);
        newe1 = alpha1 * commone + beta1 * t1;
        Ts.block<2, 1>(2 * edgeidx, 2) = BTB1.inverse() * Bs[face1].transpose() * newe1;

        alpha2 = commone.dot(e2);
        beta2 = t2.dot(e2);
        newe2 = alpha2*commone + beta2*t1;
        Ts.block<2, 1>(2 * edgeidx, 3) = BTB1.inverse() * Bs[face1].transpose() * newe2;
    }
}

int Weave::vidx(int face, int field) const
{
    return (2 * nFields()*face + 2 * field);
}

Eigen::Vector2d Weave::v(int face, int field) const
{
    return vectorFields.segment<2>(vidx(face,field));
}

int Weave::betaidx(int face, int field) const
{
    return 2 * nFields()*nFaces() + 2 * nFields()*face + 2 * field;
}
Eigen::Vector2d Weave::beta(int face, int field) const
{
    return vectorFields.segment<2>(betaidx(face,field));
}

int Weave::alphaidx(int face, int field) const
{
    return 4 * nFields()*nFaces() + nFields()*face + field;
}

double Weave::alpha(int face, int field) const
{
    return vectorFields[alphaidx(face,field)];
}

void Weave::createVisualizationEdges(Eigen::MatrixXd &edgePts, Eigen::MatrixXd &edgeVecs, Eigen::MatrixXi &edgeSegs, Eigen::MatrixXd &colors)
{
    int nfaces = nFaces();
    int m = nFields();
    int nhandles = nHandles();
    edgePts.resize(m*nfaces + nhandles, 3);
    edgeVecs.resize(m*nfaces + nhandles, 3);
    edgeVecs.setZero();
    edgeSegs.resize(m*nfaces + nhandles, 2);
    colors.resize(m*nfaces + nhandles, 3);
    
    Eigen::MatrixXd fcolors(m, 3);
    for (int i = 0; i < m; i++)
        fcolors.row(i).setZero();//heatmap(double(i), 0.0, double(m-1));

    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector3d centroid;
        centroid.setZero();
        for (int j = 0; j < 3; j++)
            centroid += V.row(F(i, j));
        centroid /= 3.0;

        for (int j = 0; j < m; j++)
        {
            edgePts.row(m*i + j) = centroid;
            edgeVecs.row(m*i + j) = Bs[i] * v(i, j);
            edgeSegs(m*i + j, 0) = 2 * (m*i + j);
            edgeSegs(m*i + j, 1) = 2 * (m*i + j) + 1;
            colors.row(m*i + j) = fcolors.row(j);
        }
    }

    for (int i = 0; i < nhandles; i++)
    {
        Eigen::Vector3d centroid;
        centroid.setZero();
        for (int j = 0; j < 3; j++)
            centroid += V.row(F(handles[i].face, j));
        centroid /= 3.0;

        Eigen::Vector3d white(1, 1, 1);
        edgePts.row(m*nfaces + i) = centroid;
        edgeVecs.row(m*nfaces + i) = Bs[handles[i].face] * handles[i].dir;
        edgeSegs(m*nfaces + i, 0) = 2 * m*nfaces + 2 * i;
        edgeSegs(m*nfaces + i, 1) = 2 * m*nfaces + 2 * i + 1;
        colors.row(m*nfaces + i).setConstant(1.0);
    }
}

bool Weave::addHandle(Handle h)
{
    if (h.face < 0 || h.face > nFaces())
        return false;
    if(h.field < 0 || h.field > nFields())
        return false;

    Eigen::Vector3d extrinsic = Bs[h.face] * h.dir;
    double mag = extrinsic.norm();
    h.dir /= mag;
    handles.push_back(h);
    return true;
}

void Weave::normalizeFields()
{
    int nfaces = nFaces();
    int m = nFields();
    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < m; j++)
        {
            Eigen::Vector2d vif = v(i, j);
            double norm = sqrt(vif.transpose() * Bs[i].transpose() * Bs[i] * vif);
            vectorFields.segment<2>(vidx(i, j)) /= norm;
        }
    }
}

using namespace std;

void Weave::removePointsFromMesh(std::vector<int> vIds)
{   
    std::set<int> facesToDelete;
    
    std::map<std::pair<int, int>, int> edgeMap;
    for (int e = 0; e < nEdges(); e++) 
    { 
        std::pair<int, int> p(edgeVerts(e,0), edgeVerts(e,1));
        edgeMap[p] = e; 
    }
   
    for (int v = 0; v < vIds.size(); v++)
    {
        for (int f = 0; f < nFaces(); f++)
        {
            for (int j = 0; j < 3; j++) 
            {
                if ( F(f, j) == vIds[v] ) 
                    facesToDelete.insert(f);       
            }
        }
    }

    std::vector<int> faceIds;
    for (std::set<int>::iterator it = facesToDelete.begin(); it != facesToDelete.end(); ++it)
        faceIds.push_back(*it);

    if (faceIds.empty())
        return;
    
    for (int i = 0; i < faceIds.size(); i++)
        cout << faceIds[i] << " " ;
    cout << "face ids \n";

    int fieldIdx = 0;
    int faceIdIdx = 0;
    int newNFaces = nFaces() - faceIds.size();
    Eigen::VectorXd vectorFields_clean = Eigen::VectorXd::Zero( 5*nFields()*newNFaces );
    Eigen::MatrixXi F_temp = Eigen::MatrixXi::Zero(newNFaces, 3); 
 
    cout << "fieldIdx \n";
    for (int i = 0; i < newNFaces; i++)   
    { 
        if ( fieldIdx == faceIds[faceIdIdx] )
        {
            fieldIdx++;
            faceIdIdx++;
            cout << fieldIdx << " ";
            i--;
            continue; // in case two ids to remove appear in sequence
        } 
        
        // vec field
        vectorFields_clean.segment(2*i*nFields(), 2*nFields()) = vectorFields.segment(2*fieldIdx*nFields(), 2*nFields());
        // beta
        vectorFields_clean.segment(2*i*nFields() + 2*newNFaces*nFields(), 2*nFields()) 
            = vectorFields.segment(2*fieldIdx*nFields() + 2*nFaces()*nFields(), 2*nFields() );
        // alpha
        vectorFields_clean.segment(i*nFields() + 4*newNFaces*nFields(), nFields()) 
            = vectorFields.segment(fieldIdx * nFields() + 4*nFaces()*nFields(), nFields() );
        // faces 
        F_temp.row(i) = F.row(fieldIdx);
        fieldIdx++;
    }

    vectorFields = vectorFields_clean;

    cout << faceIdIdx << " face id idx\n";
    

    for (int h = 0; h < handles.size(); h++)
    {
        int shift = 0;
        for (int f = 0; f < faceIds.size(); f++)
        {
            if ( handles[h].face > faceIds[f] )
                shift++;
        }
        handles[h].face = handles[h].face - shift;
    }

    Eigen::MatrixXd V_temp = V;
    Eigen::MatrixXd V_new;
    Eigen::MatrixXi F_new;
    Eigen::VectorXi marked; 
    Eigen::VectorXi vertMap; 
     
    igl::remove_unreferenced(V_temp, F_temp, V_new, F_new, marked, vertMap);
    V = V_new;
    F = F_new;


    buildConnectivityStructures();
    buildGeometricStructures();

    std::vector<Eigen::MatrixXi> Ps_new;
    for( int i = 0; i < nEdges(); i++) 
    {
        int v0 = vertMap( edgeVerts(i, 0) );
        int v1 = vertMap( edgeVerts(i, 1) );
        if ( v0 > v1 ) 
        { 
            std::swap(v0, v1);
        }
        std::pair<int, int> p(v0, v1);
        int oldEdge = edgeMap[p];
        Ps_new.push_back(Ps[oldEdge]);
    }
    Ps = Ps_new;
}

/*
 * Writes vector field to file. Format is:
 *
 * - the number of optimization variables, nvars (int)
 * - nvars doubles specifying the vector field variables, in the same format as Weave::vectorField
 * - nedges and nfields, two ints specifying the number of edges and vector fields per face
 * - nedges permutation matrices, each an nfields x nfields integer matrix, where the ith matrix corresponds to edge i
 * - the number of handles (int)
 * - for each handle, four numbers: the face of the handle (int), the field of the handle (int), and the direction of the handle,
 *   in the face's barycentric coordinates (two doubles)
 * 
 */
void Weave::serialize(const std::string &filename)
{
    std::string rawname = filename.substr(0, filename.find_last_of("."));
    std::ofstream ofs(rawname + ".relax");
    int nvars = vectorFields.size();
    ofs << nvars << std::endl;;
    for (int i = 0; i < nvars; i++)
    {
        ofs << vectorFields[i] << std::endl;;
    }

    int nedges = nEdges();
    int nfields = nFields();
    ofs << nedges << " " << nfields << std::endl;

    for (int i = 0; i < nedges; i++)
    {
        for (int j = 0; j < nfields; j++)
        {
            for (int k = 0; k < nfields; k++)
            {
                ofs << Ps[i](j, k) << " ";
            }
            ofs << std::endl;
        }
        ofs << std::endl;
    }

    int nhandles = nHandles();
    ofs << nhandles << std::endl;
    for (int i = 0; i < nhandles; i++)
    {
        ofs << handles[i].face << " " << handles[i].field << " " << handles[i].dir[0] << " " << handles[i].dir[1] << std::endl;
    }

    igl::writeOBJ(rawname + ".obj", V, F);

}


double Weave::barycentric(double val1, double val2, double target)
{
    return (target-val1) / (val2-val1);
}

bool Weave::crosses(double isoval, double val1, double val2, double minval, double maxval, double &bary)
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



int Weave::extractIsoline(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXi &faceNeighbors, const Eigen::VectorXd &func, double isoval, double minval, double maxval)
{
    int nfaces = F.rows();
    bool *visited = new bool[nfaces];
    for(int i=0; i<nfaces; i++)
        visited[i] = false;

    int ntraces = 0;

    for(int i=0; i<nfaces; i++)
    {
        if(visited[i])
            continue;
        visited[i] = true;
        std::vector<std::vector<Eigen::Vector3d> > traces;
        std::vector<std::vector<std::pair<int, int> > > traces_vids;
        for(int j=0; j<3; j++)
        {
            int vp1 = F(i, (j+1)%3);
            int vp2 = F(i, (j+2)%3);
            double bary;
            if(crosses(isoval, func[vp1], func[vp2], minval, maxval, bary))
            {
                std::vector<Eigen::Vector3d> trace;
                std::vector<std::pair<int, int> > trace_vid;
                trace.push_back( (1.0 - bary) * V.row(vp1) + bary * V.row(vp2) );
                trace_vid.push_back(std::make_pair(vp1, vp2));
                int prevface = i;
                int curface = faceNeighbors(i, j);
                while(curface != -1 && !visited[curface])
                {
                    visited[curface] = true;
                    bool found = false;
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
                            prevface = curface;
                            curface = faceNeighbors(curface, k);
                            found = true;
                            break;
                        }                       
                    }
                }
                traces.push_back(trace);
                traces_vids.push_back(trace_vid);
            }
        }
        assert(traces.size() == traces_vids.size());
        assert(traces.size() < 3);
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
        }
    }
    delete[] visited;
    return ntraces;
}

void Weave::drawISOLines(int numISOLines)
{
    double minval = -M_PI;
    double maxval = M_PI;
    double numlines = numISOLines;

    int nfaces = nFaces();
    int nverts = nVerts();
    
    std::map<std::pair<int, int>, Eigen::Vector2i > edgemap;
    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int nextj = (j + 1) % 3;
            int v1 = F(i, j);
            int v2 = F(i, nextj);
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
            int vp1 = F(i,(j+1)%3);
            int vp2 = F(i,(j+2)%3);
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
        ntraces += extractIsoline(V, F, faceNeighbors, 
            Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(theta.data(), theta.size()), 
            isoval, minval, maxval);
    }
    std::cout << ntraces << " 0 0 " << ntraces <<  " 0 0 " << std::endl;
}


vector<long> Weave::_BFS_adj_list(vector<vector<long> > adj_list, int startPoint)
{
    vector<long> traversed;
    queue<long> que;
    traversed.push_back(startPoint);
    que.push(startPoint);
    long cnt = 0;
    while (que.size() > 0)
    {
        long curPoint = que.front();
        que.pop();
        for (int j = 0; j < adj_list[curPoint].size(); j ++)
        {
            long to_add = adj_list[curPoint][j];
            bool visited = false;
            for (int i = 0; i < traversed.size(); i ++)
            {
                if (traversed[i] == to_add){
                    visited = true;
                    break;
                }
            }
            if (visited)
                continue;
            traversed.push_back(to_add);
            que.push(to_add);
        }
        cnt ++;
        if (cnt > 100000){
            cout << "++++++++++++++++++++++++++++";
            break;
        }
    }
    return traversed;
}

 std::vector<Eigen::MatrixXd> Weave::_augmentPs()
 {
    int nCover = nFields() * 2;
    int nfaces = nFaces();
    int nverts = nVerts();
    std::vector<Eigen::MatrixXd> perms;
    Eigen::MatrixXd perm;
    for (int e = 0; e < nEdges(); e++)
    {
        perm = Eigen::MatrixXd::Zero(nCover, nCover);
        for (int j = 0; j < nFields(); j++) 
        {
            for (int k = 0; k < nFields(); k++)
            {
                if( Ps[e](j, k) == 1 )
                {
                    perm(j,k) = 1;
                    perm(j+3, k+3) = 1;
                }
                if( Ps[e](j, k) == -1 )
                {
                    perm(j,k+3) = 1;
                    perm(j+3, k) = 1;
                }
            }
        }
        perms.push_back(perm);
    }
    return perms;
 }


void Weave::augmentField()
{
    std::cout << "Current fields number: " << nFields() << endl;
    int nCover = nFields() * 2;
    int nfaces = nFaces();
    int nverts = nVerts();
    std::vector<Eigen::MatrixXd> perms;
    Eigen::MatrixXd perm;
    perms = _augmentPs();
    Eigen::MatrixXd eye = Eigen::MatrixXd::Identity(nCover, nCover);

    // Compute points to glue
    vector<vector<long> > adj_list(nCover*nfaces*3);
    for (int e = 0; e < nEdges(); e++)
    {
        perm = perms[e];
        int f1Id = E(e, 0);
        int f2Id = E(e, 1);
        if(f1Id == -1 || f2Id == -1)
            continue;
        int v1ID = edgeVerts(e, 0);
        int v2ID = edgeVerts(e, 1);
        int v1f1 = -1, v2f1 = -1, v1f2 = -1, v2f2 = -1;
        for (int i = 0; i < 3; i ++)
        { // find the vid at face (0,1,or 2)
            if (F(f1Id,i) == v1ID) v1f1 = i;
            if (F(f1Id,i) == v2ID) v2f1 = i;
            if (F(f2Id,i) == v1ID) v1f2 = i;
            if (F(f2Id,i) == v2ID) v2f2 = i;
        }
        assert((v1f1 != -1) && (v2f1 != -1) && (v1f2 != -1) && (v2f2 != -1));
        if ((perm - eye).norm() == 0)
        { // perm == I case
            for (int l = 0; l < nCover; l ++)
            {
                long v1f1_idx = v1f1 + f1Id*3 + l*3*nfaces;
                long v2f1_idx = v2f1 + f1Id*3 + l*3*nfaces;
                long v1f2_idx = v1f2 + f2Id*3 + l*3*nfaces;
                long v2f2_idx = v2f2 + f2Id*3 + l*3*nfaces;
                adj_list[v1f1_idx].push_back(v1f2_idx);
                adj_list[v1f2_idx].push_back(v1f1_idx);
                adj_list[v2f1_idx].push_back(v2f2_idx);
                adj_list[v2f2_idx].push_back(v2f1_idx);
            }
        }
        else
        { // perm != I case
            cout << "Non Identity Perm!" << endl;
            for (int l1 = 0; l1 < nCover; l1 ++)
            {
                int l2 = -1;
                for (int j = 0; j < nCover; j ++)
                {
                    if (perm(l1, j) == 1)
                    {
                        l2 = j;
                        break;
                    }
                }
                long v1f1_idx = v1f1 + f1Id*3 + l1*3*nfaces;
                long v2f1_idx = v2f1 + f1Id*3 + l1*3*nfaces;
                long v1f2_idx = v1f2 + f2Id*3 + l2*3*nfaces;
                long v2f2_idx = v2f2 + f2Id*3 + l2*3*nfaces;
                adj_list[v1f1_idx].push_back(v1f2_idx);
                adj_list[v1f2_idx].push_back(v1f1_idx);
                adj_list[v2f1_idx].push_back(v2f2_idx);
                adj_list[v2f2_idx].push_back(v2f1_idx);
            }
        }
    }
    // Do some glueing
    vector<vector<long> > gluePointList;
    vector<bool> toSearchFlag(nCover*nfaces*3,1);
    for (int i = 0; i < nCover*nfaces*3; i ++)
    {
        if (i % 5000 == 0)
            cout << toSearchFlag[i] << " " << i << "/" << nCover*nfaces*3 << endl;
        if (toSearchFlag[i] == 0)
            continue;
        vector<long> gluePoint = _BFS_adj_list(adj_list, i);
        // cout << gluePoint.size() << endl;
        gluePointList.push_back(gluePoint);
        for (int j = 0; j < gluePoint.size(); j ++)
            toSearchFlag[gluePoint[j]] = 0;
    }
    int nNewPoints = gluePointList.size();
    Eigen::MatrixXd VAug = Eigen::MatrixXd::Zero(nNewPoints, 3); // |gluePointList| x 3
    vector<long> oldId2NewId(nCover*nverts);
    vector<long> encodeDOldId2NewId(nCover*3*nfaces);
    for (int i = 0; i < nNewPoints; i ++)
    { // Assign a new Vertex for each group of glue vetices
        long encodedVid = gluePointList[i][0];
        int layerId = floor(encodedVid / (nfaces*3));
        int atFace = floor((encodedVid - layerId*nfaces*3) / 3);
        int atVid = encodedVid - layerId*nfaces*3 - 3*atFace;
        int vid = F(atFace, atVid);
        for (int j = 0; j < 3; j ++)
            // VAug(i,j) = V(vid,j) + layerId;
            VAug(i,j) = V(vid,j);
        for (int j = 0; j < gluePointList[i].size(); j ++)
        { // Maintain a vid mapping
            encodedVid = gluePointList[i][j];
            layerId = floor(encodedVid / (nfaces*3));
            atFace = floor((encodedVid - layerId*nfaces*3) / 3);
            atVid = encodedVid - layerId*nfaces*3 - 3*atFace;
            assert(vid == F(atFace, atVid));
            oldId2NewId[vid + layerId*nverts] = i;
            encodeDOldId2NewId[gluePointList[i][j]] = i;
        }
    }
    Eigen::MatrixXi FAug = Eigen::MatrixXi::Zero(nCover*nfaces, 3);; // |gluePointList| x 3
    for (int cId = 0; cId < nCover; cId ++)
    {
        for (int fId = 0; fId < nfaces; fId ++)
        {
            int id0 = (fId + cId*nfaces) * 3;
            int id1 = (fId + cId*nfaces) * 3 + 1;
            int id2 = (fId + cId*nfaces) * 3 + 2;
            FAug(fId+cId*nfaces,0) = encodeDOldId2NewId[id0];
            FAug(fId+cId*nfaces,1) = encodeDOldId2NewId[id1];
            FAug(fId+cId*nfaces,2) = encodeDOldId2NewId[id2];
        }
    }
    igl::writeOBJ("debug.obj", VAug, FAug);
    cout << "finish augmenting the mesh" << endl;
    // VAugmented = VAug;
    // FAugmented = FAug;
    V = VAug;
    F = FAug;
    buildConnectivityStructures();
    buildGeometricStructures();
    augmented = true;
}

void Weave::computeFunc(double scalesInit)
{
    int nfaces = nFaces();
    int nverts = nVerts();
    cout << "nfaces: " << nfaces << endl;
    cout << "nverts: " << nverts << endl;
    vector<int> rowsL;
    vector<int> colsL;
    vector<double> difVecUnscaled;
    for (int fId = 0; fId < nfaces; fId ++)
    { // Compute rowsL, colsL, difVecUnscaled
        int vId0 = F(fId, 0);
        int vId1 = F(fId, 1);
        int vId2 = F(fId, 2);
        rowsL.push_back(vId0); rowsL.push_back(vId1); rowsL.push_back(vId2);
        colsL.push_back(vId1); colsL.push_back(vId2); colsL.push_back(vId0);
        Eigen::Vector3d p0 = V.row(vId0);
        Eigen::Vector3d p1 = V.row(vId1);
        Eigen::Vector3d p2 = V.row(vId2);
        Eigen::Vector3d e01 = p0 - p1;
        Eigen::Vector3d e12 = p1 - p2;
        Eigen::Vector3d e20 = p2 - p0;
        Eigen::Vector3d faceVec;
        if (augmented)
        {
            int oriFId = fId % (nfaces / (nFields() * 2));
            int layerId = fId / (nfaces / (nFields() * 2));
            assert((layerId >= 0) && (layerId < nFields() * 2));
            if (layerId >= 3)
                faceVec = - Bs[oriFId] * v(oriFId, layerId); // The original vec
            else
                faceVec = Bs[oriFId] * v(oriFId, layerId); // The original vec
        }
        else
            faceVec = Bs[fId] * v(fId, 0); // The original vec
        faceVec = faceVec.cross(faceNormal(fId));
        faceVec /= faceVec.norm();
        difVecUnscaled.push_back(e01.dot(faceVec));
        difVecUnscaled.push_back(e12.dot(faceVec));
        difVecUnscaled.push_back(e20.dot(faceVec));
    }
    assert((rowsL.size()==3*nfaces) && (colsL.size()==3*nfaces) && (difVecUnscaled.size()==3*nfaces));

    // Eigen::SparseMatrix<double> faceLapMat = faceLaplacian();
    Eigen::VectorXd scales(nfaces);
    scales.setConstant(scalesInit);
    int totalIter = 30;
    for (int iter = 0; iter < totalIter; iter ++)
    {
        vector<double> difVec;
        for (int i = 0; i < difVecUnscaled.size(); i ++)
            difVec.push_back(difVecUnscaled[i]*scales(i/3));
        std::vector<triplet> sparseContent;
        for (int i = 0; i < rowsL.size(); i ++)
            sparseContent.push_back(triplet(rowsL[i],colsL[i],1));
        Eigen::SparseMatrix<double> TP (nverts, nverts);
        Eigen::SparseMatrix<double> TPTran (nverts, nverts);
        TP.setFromTriplets(sparseContent.begin(),sparseContent.end());
        TPTran = TP.transpose();
        TP += TPTran;
        vector<int> degree;
        for (int i = 0; i < nverts; i ++)
            degree.push_back(TP.row(i).sum());
        std::vector<triplet> AContent;
        for (int i = 0; i < rowsL.size(); i ++)
        {
            double cVal = cos(difVec[i]);
            double sVal = sin(difVec[i]);
            AContent.push_back(triplet(2*rowsL[i], 2*colsL[i], cVal));
            AContent.push_back(triplet(2*rowsL[i], 2*colsL[i]+1, -sVal));
            AContent.push_back(triplet(2*rowsL[i]+1, 2*colsL[i], sVal));
            AContent.push_back(triplet(2*rowsL[i]+1, 2*colsL[i]+1, cVal));
        }
        Eigen::SparseMatrix<double> Amat (2*nverts, 2*nverts);
        Eigen::SparseMatrix<double> Amat_tran (2*nverts, 2*nverts);
        Amat.setFromTriplets(AContent.begin(),AContent.end());
        Amat_tran = Amat.transpose();
        Amat += Amat_tran;
        //
        std::vector<triplet> LContent;
        for (int i = 0; i < 2*nverts; i ++)
            LContent.push_back(triplet(i,i,degree[int(i/2)]));
        Eigen::SparseMatrix<double> Lmat (2*nverts, 2*nverts);
        Lmat.setFromTriplets(LContent.begin(), LContent.end());
        Lmat -= Amat;
        // Eigen Decompose
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solverL(Lmat);
        Eigen::VectorXd ones(Lmat.rows());
        ones.setConstant(1.0);
        ones /= ones.norm();
        Eigen::VectorXd eigenVec(Lmat.rows());
        eigenVec.setRandom();
        eigenVec += ones;
        eigenVec /= eigenVec.norm();
        // for (int k=0; k<Amat.outerSize(); ++k)
        //   for (Eigen::SparseMatrix<double>::InnerIterator it(Amat,k); it; ++it)
        //   {
        //     cout << it.row() << " ";
        //     cout << it.col() << " ";
        //     cout << it.value() << " ";
        //     cout << endl;
        //   }
        for(int i=0; i<30; i++)
        {
            eigenVec = solverL.solve(eigenVec);
            Eigen::VectorXd proj = eigenVec.dot(ones) * ones;
            eigenVec -= proj;
            eigenVec /= eigenVec.norm();
        }
        double eigenVal = eigenVec.transpose() * Lmat * eigenVec;
        cout << "Current iteration = " << iter  << " currents error is: " << eigenVal << endl;
        // Extract the function value
        vector<double> curTheta;
        for (int i = 0; i < nverts; i ++)
        {
            double curCos = eigenVec(2*i);
            double curSin = eigenVec(2*i+1);
            double normalizer = sqrt(curCos * curCos + curSin * curSin);
            double curFunc = acos(curCos / normalizer);
            if (curSin < 0)
                curFunc = -curFunc;
            curTheta.push_back(curFunc);
        }
        ////
        //// Re-compute face scales
        vector<double> difVecPred;
        for (int i = 0; i < rowsL.size(); i ++)
        {
            double curPred = curTheta[rowsL[i]] - curTheta[colsL[i]];
            if (curPred > M_PI) curPred -= 2*M_PI;
            if (curPred < -M_PI) curPred += 2*M_PI;
            difVecPred.push_back(curPred);
        }
        vector<double> bScales;
        vector<double> diagAScales;
        for (int i = 0; i < rowsL.size(); i=i+3)
        {
            double bVal = 0;
            double diagAVal = 0;
            for (int j = 0; j < 3; j ++)
            {
                bVal += difVecPred[i+j] * difVecUnscaled[i+j];
                diagAVal += difVecUnscaled[i+j] * difVecUnscaled[i+j];

            }
            bScales.push_back(bVal);
            diagAScales.push_back(diagAVal);
        }
        // Construct A
        // TODO mu and lambda
        std::vector<triplet> AScalesContent;
        for (int i = 0; i < nfaces; i ++)
            AScalesContent.push_back(triplet(i, i, diagAScales[i]));
        Eigen::SparseMatrix<double> AScalesMat (nfaces, nfaces);
        AScalesMat.setFromTriplets(AScalesContent.begin(),AScalesContent.end());
        // Solve for scale
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solverScales(AScalesMat);
        scales = solverScales.solve(
            Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(bScales.data(), bScales.size()));
        theta = curTheta;
    }
    for (int i = 0; i < nverts; i ++)
        cout << theta[i] << endl;;
}

Eigen::SparseMatrix<double> Weave::faceLaplacian()
{ // Only augment vector
    int nfaces = nFaces();
    // TODO: boundary
    // ids = find(min(adjFaces) > 0);
    // adjFaces = adjFaces(:, ids);
    std::vector<triplet> AContent;
    for (int i = 0; i < E.rows(); i ++)
        AContent.push_back(triplet(E(i,0), E(i,1), 1));
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


/*
 * Writes vector field to file. Format is:
 *
 * Number of vector fields (int) = |F|*m (m = fields per face). 
 * |F| x 3m matrix of vectors
 * Number of edges |E| (int) m
 * 1 row recording edge adjacency information (0,1 are adjacent faces, 2,3 are adjacent verts)
 * mxm permutation matrix
 * 
 */
void Weave::serialize_forexport(const std::string &filename)
{
    char buffer [100];
    sprintf(buffer, "%s.fields", filename.c_str());

    std::ofstream ofs(buffer);
    int nvars = vectorFields.size();
  //  ofs << Bs.size() << std::endl;
    for (int i = 0; i < nFaces(); i++)
    {
        for (int j = 0; j < nFields(); j++)
        {
            ofs << (Bs[i] * v(i, j)).transpose() << " "; 
        }
    
        for (int j = 0; j < nFields(); j++)
        {
            ofs << -(Bs[i] * v(i, j)).transpose() << " "; 
        }

        ofs << std::endl;
    }
    ofs.close();


    sprintf(buffer, "%s.edges", filename.c_str());
    std::ofstream ofs_edge(buffer);

    int nedges = nEdges();
    int nfields = nFields();
  //  ofs << nedges << " " << nfields << std::endl;

    for (int i = 0; i < nedges; i++)
    {
        ofs_edge << E(i, 0) + 1 << " " 
                 << E(i, 1) + 1 << " " 
                 << edgeVerts(i, 0) + 1  << " "
                 << edgeVerts(i, 1) + 1 <<  std::endl;
    }
    ofs_edge.close();
    sprintf(buffer, "%s.permmats", filename.c_str());
    std::ofstream ofs_mat(buffer);

    for (int i = 0; i < nedges; i++)
    {
        Eigen::MatrixXd perm = Eigen::MatrixXd::Zero(nfields * 2, nfields*2);
        for (int j = 0; j < nfields; j++) 
        {
            for (int k = 0; k < nfields; k++)
            {
                if( Ps[i](j, k) == 1 )
                {
                    perm(j,k) = 1;
                    perm(j+3, k+3) = 1;
                }

                if( Ps[i](j, k) == -1 )
                {
                    perm(j,k+3) = 1;
                    perm(j+3, k) = 1;
                }
            }
        }
 

        for (int j = 0; j < nfields * 2; j++)
        {
            for (int k = 0; k < nfields * 2; k++)
            {
                ofs_mat << perm(j, k) << " ";
            }
            ofs_mat << std::endl;
        }
    }
    ofs_mat.close();
}


void Weave::deserialize(const std::string &filename)
{
    std::ifstream ifs(filename);
    int nvars;
    ifs >> nvars;
    if (!ifs)
    {
        std::cerr << "Couldn't load vector field file " << filename << std::endl;
        return;
    }

    if (nvars != vectorFields.size())
    {
        std::cerr << "Vector field doesn't match mesh!" << std::endl;
        return;
    }

    for (int i = 0; i < nvars; i++)
        ifs >> vectorFields[i];

    int nedges, nfields;
    ifs >> nedges >> nfields;
    if (!ifs)
    {
        std::cerr << "Error reading vector field file " << filename << std::endl;
        return;
    }

    if (nedges != nEdges() && nfields != nFields())
    {
        std::cerr << "Vector field doesn't match mesh! edge/fields wrong." << std::endl;
        return;
    }

    for (int i = 0; i < nedges; i++)
    {
        for (int j = 0; j < nfields; j++)
        {
            for (int k = 0; k < nfields; k++)
            {
                ifs >> Ps[i](j, k);
            }
        }
    }

    int nhandles;
    ifs >> nhandles;
    if (!ifs)
    {
        std::cerr << "Error reading vector field file " << filename << std::endl;
        return;
    }
    handles.clear();

    for (int i = 0; i < nhandles; i++)
    {
        Handle h;
        ifs >> h.face >> h.field >> h.dir[0] >> h.dir[1];
        handles.push_back(h);
    }
    if (!ifs)
    {
        std::cerr << "Error reading the vector field file " << filename << std::endl;
    }
}

void Weave::shortestPath(int startVert, int endVert, std::vector<std::pair<int, int> > &path)
{
    int nverts = nVerts();
    path.clear();
    bool *visited = new bool[nverts];
    for(int i=0; i<nverts; i++)
        visited[i] = false;
    Eigen::VectorXi prev(nverts);
    Eigen::VectorXi prevOrient(nverts);
    Eigen::VectorXd prevEdge(nverts);
    prev.setConstant(-1);
    prevOrient.setConstant(-1);
    prevEdge.setConstant(-1);

    struct SearchNode
    {
        int next;
        int prev;
        int orient;
        int prevedge;
    };

    SearchNode start;
    start.next = startVert;
    start.prev = -1;
    start.orient = -1;
    start.prevedge = -1;
    
    std::deque<SearchNode> q;
    q.push_back(start);

    while(!q.empty())
    {
        SearchNode cur = q.front();
        q.pop_front();
        if(visited[cur.next])
            continue;
        visited[cur.next] = true;
        prev[cur.next] = cur.prev;
        prevOrient[cur.next] = cur.orient;
        prevEdge[cur.next] = cur.prevedge;
        
        if(cur.next == endVert)
        {
            int v = endVert;
            while(prev[v] != -1)
            {
                path.push_back(std::pair<int, int>(prevEdge[v], prevOrient[v]));
                v = prev[v];
            }
            std::reverse(path.begin(), path.end());
            delete[] visited;
            return;
        }

        int nbs = (int)vertEdges[cur.next].size();
        for(int i=0; i<nbs; i++)
        {
            int e = vertEdges[cur.next][i];
            int orient = (edgeVerts(e, 0) == cur.next) ? 0 : 1;
            int next = edgeVerts(e, 1-orient);
            if(visited[next])
                continue;
            SearchNode nextnode;
            nextnode.next = next;
            nextnode.prev = cur.next;
            nextnode.prevedge = e;
            nextnode.orient = orient;
            q.push_back(nextnode);
        }
    }

    delete[] visited;
}


void Weave::createVisualizationCuts(Eigen::MatrixXd &cutPts1, Eigen::MatrixXd &cutPts2)
{
    int totedges = 0;
    for (int i = 0; i < (int)cuts.size(); i++)
    {
        totedges += cuts[i].path.size();
    }
    cutPts1.resize(totedges, 3);
    cutPts2.resize(totedges, 3);
    int idx = 0;
    for (int i = 0; i < (int)cuts.size(); i++)
    {
        for (int j = 0; j < (int)cuts[i].path.size(); j++)
        {
            int f1 = E(cuts[i].path[j].first, 0);
            int f2 = E(cuts[i].path[j].first, 1);
            Eigen::Vector3d n1 = faceNormal(f1);
            Eigen::Vector3d n2 = faceNormal(f2);
            Eigen::Vector3d offset = 0.0001*(n1 + n2);
            cutPts1.row(idx) = V.row(edgeVerts(cuts[i].path[j].first, 0)) + offset.transpose();
            cutPts2.row(idx) = V.row(edgeVerts(cuts[i].path[j].first, 1)) + offset.transpose();
            idx++;
        }
    }
}

int Weave::numInteriorEdges() const
{
    int nedges = nEdges();
    int result = 0;
    for(int i=0; i<nedges; i++)
    {
        if(E(i,0) != -1 && E(i,1) != -1)
            result++;
    }
    return result;
}
