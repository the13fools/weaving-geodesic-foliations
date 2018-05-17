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
    delete fs;
    fs = FieldSurface::deserialize(ifs);

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
