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
#include "Surface.h"
#include "CoverMesh.h"

Weave::Weave(const std::string &objname, int m)
{
    Eigen::MatrixXd Vtmp;
    Eigen::MatrixXi Ftmp;
    if (!igl::read_triangle_mesh(objname, Vtmp, Ftmp))
    {
        std::cerr << "Couldn't load mesh " << objname << std::endl;
        exit(-1);
    }
    if (Vtmp.cols() < 3)
    {
        std::cerr << "Mesh must 3D" << std::endl;
        exit(-1);
    }
    
    centerAndScale(Vtmp);
    fs = new FieldSurface(Vtmp, Ftmp, m);       
}

Weave::~Weave()
{    
    delete fs;
}

void Weave::centerAndScale(Eigen::MatrixXd &V)
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

void Weave::createVisualizationEdges(Eigen::MatrixXd &edgePts, Eigen::MatrixXd &edgeVecs, Eigen::MatrixXi &edgeSegs, Eigen::MatrixXd &colors)
{
    int nfaces = fs->nFaces();
    int m = fs->nFields();
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
            centroid += fs->data().V.row(fs->data().F(i, j));
        centroid /= 3.0;

        for (int j = 0; j < m; j++)
        {
            edgePts.row(m*i + j) = centroid;
            edgeVecs.row(m*i + j) = fs->data().Bs[i] * fs->v(i, j);
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
            centroid += fs->data().V.row(fs->data().F(handles[i].face, j));
        centroid /= 3.0;

        Eigen::Vector3d white(1, 1, 1);
        edgePts.row(m*nfaces + i) = centroid;
        edgeVecs.row(m*nfaces + i) = fs->data().Bs[handles[i].face] * handles[i].dir;
        edgeSegs(m*nfaces + i, 0) = 2 * m*nfaces + 2 * i;
        edgeSegs(m*nfaces + i, 1) = 2 * m*nfaces + 2 * i + 1;
        colors.row(m*nfaces + i).setConstant(1.0);
    }
}

bool Weave::addHandle(Handle h)
{
    if (h.face < 0 || h.face > fs->nFaces())
        return false;
    if(h.field < 0 || h.field > fs->nFields())
        return false;

    Eigen::Vector3d extrinsic = fs->data().Bs[h.face] * h.dir;
    double mag = extrinsic.norm();
    h.dir /= mag;
    handles.push_back(h);
    return true;
}

using namespace std;

void Weave::removePointsFromMesh(std::vector<int> vIds)
{   
    std::map<int, int> facemap;
    FieldSurface *newfs = fs->removePointsFromMesh(vIds, facemap);

    std::vector<Handle> newhandles;
    for (int h = 0; h < handles.size(); h++)
    {
        Handle hand = handles[h];

        auto it = facemap.find(hand.face);
        if (it != facemap.end())
        {
            hand.face = it->second;
            newhandles.push_back(hand);
        }
    }
    delete fs;
    fs = newfs;
    handles = newhandles;

    // for now, clear (now stale) cuts
    cuts.clear();
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
void Weave::serialize(std::ostream &ofs)
{
    fs->serialize(ofs);

    int nhandles = nHandles();
    ofs.write((char *)&nhandles, sizeof(int));
    for (int i = 0; i < nhandles; i++)
    {
        ofs.write((char *)&handles[i].face, sizeof(int));
        ofs.write((char *)&handles[i].field, sizeof(int));
        double d0 = handles[i].dir[0];
        ofs.write((char *)&d0, sizeof(double));
        double d1 = handles[i].dir[1];
        ofs.write((char *)&d1, sizeof(double));        
    }    

    int ncuts = cuts.size();
    ofs.write((char *)&ncuts, sizeof(int));
    for (int i = 0; i < ncuts; i++)
    {
        int nlinks = cuts[i].path.size();
        ofs.write((char *)&nlinks, sizeof(int));
        for (int j = 0; j < nlinks; j++)
        {
            ofs.write((char *)&cuts[i].path[j].first, sizeof(int));
            ofs.write((char *)&cuts[i].path[j].second, sizeof(int));
        }
    }
}

void Weave::deserialize(std::istream &ifs)
{
    FieldSurface *newfs = FieldSurface::deserialize(ifs);
    if (!newfs)
        return;

    delete fs;
    fs = newfs;


    handles.clear();
    int nhandles;
    ifs.read((char *)&nhandles, sizeof(int));
    for (int i = 0; i < nhandles; i++)
    {
        Handle h;
        ifs.read((char *)&h.face, sizeof(int));
        ifs.read((char *)&h.field, sizeof(int));
        double d0;
        double d1;
        ifs.read((char *)&d0, sizeof(double));
        ifs.read((char *)&d1, sizeof(double));
        h.dir[0] = d0;
        h.dir[1] = d1;
        handles.push_back(h);
    }

    cuts.clear();
    int ncuts;
    ifs.read((char *)&ncuts, sizeof(int));
    for (int i = 0; i < ncuts; i++)
    {
        int nlinks;
        ifs.read((char *)&nlinks, sizeof(int));
        Cut c;
        for (int j = 0; j < nlinks; j++)
        {
            int f;
            ifs.read((char *)&f, sizeof(int));
            int s;
            ifs.read((char *)&s, sizeof(int));
            c.path.push_back(std::pair<int, int>(f, s));
        }
        cuts.push_back(c);
    }
}

void Weave::deserializeOldRelaxFile(std::istream &ifs)
{    
    int nvars;
    ifs >> nvars;
    if (!ifs)
    {
        std::cerr << "Couldn't load vector field file" << std::endl;
        return;
    }

    if (nvars != fs->vectorFields.size())
    {
        std::cerr << "Vector field doesn't match mesh!" << std::endl;
        return;
    }

    for (int i = 0; i < nvars; i++)
        ifs >> fs->vectorFields[i];

    int nedges, nfields;
    ifs >> nedges >> nfields;
    if (!ifs)
    {
        std::cerr << "Error reading vector field file" << std::endl;
        return;
    }

    if (nedges != fs->nEdges() && nfields != fs->nFields())
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
                ifs >> fs->Ps_[i](j, k);
            }
        }
    }

    int nhandles;
    ifs >> nhandles;
    if (!ifs)
    {
        std::cerr << "Error reading vector field file" << std::endl;
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
        std::cerr << "Error reading the vector field file" << std::endl;
    }
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
            int f1 = fs->data().E(cuts[i].path[j].first, 0);
            int f2 = fs->data().E(cuts[i].path[j].first, 1);
            Eigen::Vector3d n1 = fs->faceNormal(f1);
            Eigen::Vector3d n2 = fs->faceNormal(f2);
            Eigen::Vector3d offset = 0.0001*(n1 + n2);
            cutPts1.row(idx) = fs->data().V.row(fs->data().edgeVerts(cuts[i].path[j].first, 0)) + offset.transpose();
            cutPts2.row(idx) = fs->data().V.row(fs->data().edgeVerts(cuts[i].path[j].first, 1)) + offset.transpose();
            idx++;
        }
    }
}


CoverMesh *Weave::createCover() const
{
    int nCover = fs->nFields() * 2;
    int nfaces = fs->nFaces();
    int nverts = fs->nVerts();
    std::vector<Eigen::MatrixXd> perms;
    Eigen::MatrixXd perm;
    perms = _augmentPs();
    Eigen::MatrixXd eye = Eigen::MatrixXd::Identity(nCover, nCover);
    // Compute points to glue
    vector<vector<long> > adj_list(nCover*nfaces*3);
    for (int e = 0; e < fs->nEdges(); e++)
    {
        perm = perms[e];
        int f1Id = fs->data().E(e, 0);
        int f2Id = fs->data().E(e, 1);
        if(f1Id == -1 || f2Id == -1)
            continue;
        int v1ID = fs->data().edgeVerts(e, 0);
        int v2ID = fs->data().edgeVerts(e, 1);
        int v1f1 = -1, v2f1 = -1, v1f2 = -1, v2f2 = -1;
        for (int i = 0; i < 3; i ++)
        { // find the vid at face (0,1,or 2)
            if (fs->data().F(f1Id,i) == v1ID) v1f1 = i;
            if (fs->data().F(f1Id,i) == v2ID) v2f1 = i;
            if (fs->data().F(f2Id,i) == v1ID) v1f2 = i;
            if (fs->data().F(f2Id,i) == v2ID) v2f2 = i;
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
            for (int l1 = 0; l1 < nCover; l1 ++)
            {
                int l2 = -1;
                for (int j = 0; j < nCover; j ++)
                    if (perm(l1, j) == 1){ l2 = j; break; }
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
        gluePointList.push_back(gluePoint);
        for (int j = 0; j < gluePoint.size(); j ++)
            toSearchFlag[gluePoint[j]] = 0;
    }
    int nNewPoints = gluePointList.size();
    Eigen::MatrixXd VAug = Eigen::MatrixXd::Zero(nNewPoints, 3); // |gluePointList| x 3
    Eigen::VectorXi oldId2NewId(nCover*nverts);
    oldId2NewId.setConstant(-1);
    vector<long> encodeDOldId2NewId(nCover*3*nfaces);
    for (int i = 0; i < nNewPoints; i ++)
    { // Assign a new Vertex for each group of glue vetices
        long encodedVid = gluePointList[i][0];
        int layerId = floor(encodedVid / (nfaces*3));
        int atFace = floor((encodedVid - layerId*nfaces*3) / 3);
        int atVid = encodedVid - layerId*nfaces*3 - 3*atFace;
        int vid = fs->data().F(atFace, atVid);
        for (int j = 0; j < 3; j ++)
            VAug(i,j) = fs->data().V(vid,j);
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
    Eigen::MatrixXd flattenedField(nCover*nfaces, 2);
    for (int cId = 0; cId < nCover; cId++)
    {
        for (int fId = 0; fId < nfaces; fId++)
        {
            int field = cId % fs->nFields();
            double sign = (cId < fs->nFields() ? 1.0 : -1.0);
            Eigen::Vector2d vec = sign*fs->v(fId, field).transpose();
            Eigen::Vector3d embvec = fs->data().Bs[fId] * vec;
            double norm = embvec.norm();
            flattenedField.row(fId + cId * nfaces) = vec/norm;
        }
    }
    //igl::writeOBJ("debug.obj", VAug, FAug);
    cout << "finish augmenting the mesh" << endl;
    CoverMesh *ret = new CoverMesh(*this, VAug, FAug, oldId2NewId, flattenedField, nCover); 
    return ret;
}


std::vector<long> Weave::_BFS_adj_list(std::vector<std::vector<long> > & adj_list, int startPoint) const
{
    vector<long> traversed;
    queue<long> que;
    traversed.push_back(startPoint);
    que.push(startPoint);
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
    }
    return traversed;
}

std::vector<Eigen::MatrixXd> Weave::_augmentPs() const
{
    int nCover = fs->nFields() * 2;
    int nfaces = fs->nFaces();
    int nverts = fs->nVerts();
    std::vector<Eigen::MatrixXd> perms;
    Eigen::MatrixXd perm;
    for (int e = 0; e < fs->nEdges(); e++)
    {
        perm = Eigen::MatrixXd::Zero(nCover, nCover);
        for (int j = 0; j < fs->nFields(); j++) 
        {
            for (int k = 0; k < fs->nFields(); k++)
            {
                if( fs->Ps(e)(j, k) == 1 )
                {
                    perm(j,k) = 1;
                    perm(j+3, k+3) = 1;
                }
                if( fs->Ps(e)(j, k) == -1 )
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
