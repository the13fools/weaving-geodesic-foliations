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
#include "RoSyUtils.h"
#include <Eigen/Geometry>
#include "Permutations.h"

Weave::Weave(const std::string &objname, int m)
{
    Eigen::MatrixXd Vtmp;
    Eigen::MatrixXi Ftmp;
    if (!igl::read_triangle_mesh(objname, Vtmp, Ftmp))
    {    
        std::cerr << "Couldn't load mesh " << objname << std::endl;
        std::string modname = "../" + objname;
        std::cerr << "Trying " << modname << " instead" << std::endl;
        if (!igl::read_triangle_mesh(modname, Vtmp, Ftmp))
        {
            std::cerr << "Couldn't load mesh " << modname << " either" << std::endl;
            exit(-1);
        }
    }
    if (Vtmp.cols() < 3)
    {
        std::cerr << "Mesh must 3D" << std::endl;
        exit(-1);
    }
    
    centerAndScale(Vtmp);
    fs = new FieldSurface(Vtmp, Ftmp, m);       
}

Weave::Weave(Eigen::MatrixXd Vtmp, Eigen::MatrixXi Ftmp, int m)
{
    centerAndScale(Vtmp);
    fs = new FieldSurface(Vtmp, Ftmp, m);       
}

Weave::Weave(const Weave &w)
{
    fs = nullptr;
    operator=(w);
}

Weave& Weave::operator=(const Weave &w)
{
    delete fs;
    handles = w.handles;
    cuts = w.cuts;
    fixFields = w.fixFields;
    fs = new FieldSurface(w.fs->data().V, w.fs->data().F, w.fs->nFields());
    fs->vectorFields = w.fs->vectorFields;
    fs->Ps_ = w.fs->Ps_;

    return *this;
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

void Weave::createVisualizationEdges(
    Eigen::MatrixXd &edgePts, 
    Eigen::MatrixXi &edgeSegs, 
    Eigen::MatrixXd &colors,
    VectorVisualizationMode mode,
    bool normalizeVectors,
    double baseVectorLength) // ignored if normalizeVectors=true
{
    int nfaces = fs->nFaces();
    int m = fs->nFields();
    int nhandles = nHandles();
    int vecsperface = 0;
    if (mode == VMM_VF || mode == VMM_VFPLUSDELTA)
        vecsperface = 1;
    else if (mode == VMM_VFANDDELTA)
        vecsperface = 2;
    edgePts.resize(2 * m * nfaces * vecsperface + 2 * nhandles, 3);
    edgePts.setZero();
    edgeSegs.resize(m*nfaces * vecsperface + nhandles, 2);
    edgeSegs.setZero();
    colors.resize(m*nfaces * vecsperface + nhandles, 3);

    Eigen::Vector3d deltacolor = Eigen::Vector3d(1,0,0);
    //Eigen::VectorXd dual = fs->vectorFields.segment(2*nfaces*m, 2*nfaces*m);
    
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
            // regular vector field
            if (mode == VMM_VF || mode == VMM_VFANDDELTA || mode == VMM_VFPLUSDELTA)
            {
                Eigen::Vector3d vec = fs->data().Bs[i] * fs->v(i, j);
                if(mode == VMM_VFPLUSDELTA)
                    vec += fs->data().Bs[i] * fs->beta(i, j); // this is actually delta now...
                if (normalizeVectors)
                    vec.normalize();
                else
                {
                    vec *= baseVectorLength;
                }
                vec *= fs->data().averageEdgeLength * sqrt(3.0) / 6.0 * 0.75;

                edgePts.row(2 * m * i + 2 * j) = centroid.transpose();
                edgePts.row(2 * m * i + 2 * j + 1) = (centroid + vec).transpose();
                edgeSegs(m*i + j, 0) = 2 * (m*i + j);
                edgeSegs(m*i + j, 1) = 2 * (m*i + j) + 1;
                colors.row(m*i + j) = fcolors.row(j);
            }
            if (mode == VMM_VFANDDELTA)
            {
                // include delta vectors
                Eigen::Vector3d delta = fs->data().Bs[i] * fs->beta(i, j); // this is actually delta now...
                if (normalizeVectors)
                    delta.normalize();
                else
                {
                    delta *= baseVectorLength;
                }
                delta *= fs->data().averageEdgeLength * sqrt(3.0) / 6.0 * 0.75;

                edgePts.row(2 * m * nfaces + 2 * m*i + 2 * j) = centroid.transpose();
                edgePts.row(2 * m * nfaces + 2 * m*i + 2 * j + 1) = (centroid + delta).transpose();
                edgeSegs(m*nfaces + m * i + j, 0) = 2 * (m*i + j) + 2 * m*nfaces;
                edgeSegs(m*nfaces + m * i + j, 1) = 2 * (m*i + j) + 1 + 2 * m*nfaces;
                colors.row(m*i + j + m * nfaces) = deltacolor;
            }
        }
    }

    for (int i = 0; i < nhandles; i++)
    {
        Eigen::Vector3d centroid;
        centroid.setZero();
        for (int j = 0; j < 3; j++)
            centroid += fs->data().V.row(fs->data().F(handles[i].face, j));
        centroid /= 3.0;

        Eigen::Vector3d dir = fs->data().Bs[handles[i].face] * handles[i].dir;
        dir *= fs->data().averageEdgeLength * sqrt(3.0) / 6.0 * 0.75  / dir.norm();
        Eigen::Vector3d white(1, 1, 1);
        edgePts.row(2 * m*nfaces*vecsperface + 2 * i) = centroid.transpose();
        edgePts.row(2 * m*nfaces*vecsperface + 2 * i + 1) = (centroid+dir).transpose();
        edgeSegs(m*nfaces*vecsperface + i, 0) = 2 * m*nfaces *vecsperface + 2 * i;
        edgeSegs(m*nfaces*vecsperface + i, 1) = 2 * m*nfaces *vecsperface + 2 * i + 1;
        colors.row(m*nfaces*vecsperface + i).setConstant(1.0);        
    }
}

void Weave::createVisualizationEdges(
    Eigen::MatrixXd &edgePts, 
    Eigen::MatrixXi &edgeSegs, 
    Eigen::MatrixXd &colors,
    RoSyVisualizationMode mode,
    bool normalizeVectors,
    double baseVectorLength,
    int rosyN)
{
    int nfaces = fs->nFaces();
    int nhandles = nHandles();
    int vecsperface = 0;
    if (mode == RVM_REPVEC)
        vecsperface = 1;
    else if (mode == RVM_ROSY)
        vecsperface = rosyN;
    edgePts.resize(2 * nfaces * vecsperface + 2 * nhandles, 3);
    edgePts.setZero();
    edgeSegs.resize(nfaces * vecsperface + nhandles, 2);
    edgeSegs.setZero();
    colors.resize(nfaces * vecsperface + nhandles, 3);

    Eigen::Vector3d fcolor(0,0,0);
    
    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector3d centroid;
        centroid.setZero();
        for (int j = 0; j < 3; j++)
            centroid += fs->data().V.row(fs->data().F(i, j));
        centroid /= 3.0;

        // regular vector field
        if (mode == RVM_REPVEC)
        {
            Eigen::Vector3d vec = fs->data().Bs[i] * fs->v(i, 0);
            if (normalizeVectors)
                vec.normalize();
            else
            {
                vec *= baseVectorLength;
            }
            vec *= fs->data().averageEdgeLength * sqrt(3.0) / 6.0 * 0.75;

            edgePts.row(2 * i) = centroid.transpose();
            edgePts.row(2 * i + 1) = (centroid + vec).transpose();
            edgeSegs(i, 0) = 2 * i;
            edgeSegs(i, 1) = 2 * i + 1;
            colors.row(i) = fcolor;
        }
        if (mode == RVM_ROSY)
        {
            std::vector<Eigen::Vector2d> vecs;
            repVecToRoSy(*fs, i, fs->v(i, 0), vecs, rosyN);
            for (int j = 0; j < rosyN; j++)
            {
                Eigen::Vector3d extv = fs->data().Bs[i] * vecs[j];

                if (normalizeVectors)
                    extv.normalize();
                else
                {
                    extv *= baseVectorLength;
                }
                extv *= fs->data().averageEdgeLength * sqrt(3.0) / 6.0 * 0.75;

                edgePts.row(2 * rosyN * i + 2 * j) = centroid.transpose();
                edgePts.row(2 * rosyN * i + 2 * j + 1) = (centroid + extv).transpose();
                edgeSegs(rosyN * i + j, 0) = 2 * (rosyN * i + j);
                edgeSegs(rosyN * i + j, 1) = 2 * (rosyN * i + j) + 1;
                colors.row(rosyN * i + j) = fcolor;
            }
        }
    }

    for (int i = 0; i < nhandles; i++)
    {
        Eigen::Vector3d centroid;
        centroid.setZero();
        for (int j = 0; j < 3; j++)
            centroid += fs->data().V.row(fs->data().F(handles[i].face, j));
        centroid /= 3.0;

        Eigen::Vector3d dir = fs->data().Bs[handles[i].face] * handles[i].dir;
        dir *= fs->data().averageEdgeLength * sqrt(3.0) / 6.0 * 0.75  / dir.norm();
        Eigen::Vector3d white(1, 1, 1);
        edgePts.row(2 * nfaces*vecsperface + 2 * i) = centroid.transpose();
        edgePts.row(2 * nfaces*vecsperface + 2 * i + 1) = (centroid+dir).transpose();
        edgeSegs(nfaces*vecsperface + i, 0) = 2 * nfaces *vecsperface + 2 * i;
        edgeSegs(nfaces*vecsperface + i, 1) = 2 * nfaces *vecsperface + 2 * i + 1;
        colors.row(nfaces*vecsperface + i).setConstant(1.0);        
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

void Weave::deserializePaulFile(std::ifstream &ifs)
{
    FieldSurface *newfs = new FieldSurface(fs->data().V, fs->data().F, 2);    
    int nfaces = fs->nFaces();
    for(int i=0; i<nfaces; i++)
    {
        char comma;
        Eigen::Vector3d vecs[3];
        for(int j=0; j<3; j++)
        {
            ifs >> vecs[0][j];
            ifs >> comma;
        }
        for(int j=0; j<3; j++)
        {
            ifs >> vecs[1][j];
            ifs >> comma;
        }
        for(int j=0; j<3; j++)
        {
            ifs >> vecs[2][j];
            if(j != 2)
                ifs >> comma;
        }
        Eigen::Matrix<double, 3, 2> B = fs->data().Bs[i];
        Eigen::Vector3d n = fs->faceNormal(i);
        
        Eigen::Matrix2d BTB = B.transpose() * B;
        
        int idx=0;
        for(int j=0; j<3 && idx<2; j++)
        {
            if(fabs(n.dot(vecs[j])) > 0.5)
                continue;
                
            Eigen::Vector2d vint = BTB.inverse() * B.transpose() * vecs[j];
            newfs->vectorFields.segment<2>(newfs->vidx(i,idx)) = vint;
            idx++;
        }
    }
    if(ifs)
    {
        delete fs;
        fs = newfs;
    }
    else
    {
        delete newfs;
    }
}

static Eigen::Vector3d parallelTransport(const Eigen::Vector3d &v, const Eigen::Vector3d &e1, const Eigen::Vector3d &e2)
{
    Eigen::Vector3d t1 = e1 / e1.norm();
    Eigen::Vector3d t2 = e2 / e2.norm();
    Eigen::Vector3d n = t1.cross(t2);
    if (n.norm() < 1e-8)
        return v;
    n /= n.norm();
    Eigen::Vector3d p1 = n.cross(t1);
    Eigen::Vector3d p2 = n.cross(t2);
    return v.dot(n)*n + v.dot(t1)*t2 + v.dot(p1)*p2;
}

static double angle(Eigen::Vector3d v, Eigen::Vector3d w, Eigen::Vector3d n)
{
    return 2.0 * atan2( v.cross(w).dot(n), v.norm() * w.norm() + v.dot(w));
}

static void alignFrame(std::vector<Eigen::Vector3d> &from, const std::vector<Eigen::Vector3d> &to, const Eigen::Vector3d &n)
{
    int m = to.size();
    std::vector<int> cand;
    std::vector<int> bestperm;
    double bestangle = std::numeric_limits<double>::infinity();
    for(int i=0; i<m; i++)
        cand.push_back(i);
    do 
    {
        double totangle = 0;
        for(int j=0; j<m; j++)
        {
            Eigen::Vector3d v1 = to[j];
            Eigen::Vector3d v2 = from[cand[j]];
            double theta = angle(v1, v2, n);
            totangle += theta*theta;            
        }
        if(totangle < bestangle)
        {
            bestangle = totangle;
            bestperm = cand;
        }
    } 
    while ( std::next_permutation(cand.begin(), cand.end()) );
    
    std::vector<Eigen::Vector3d> newvecs;
    for(int i=0; i<m; i++)
        newvecs.push_back(from[bestperm[i]]);
    from = newvecs;
}

void Weave::deserializeVertexFile(std::ifstream &ifs)
{
    int numfields;
    ifs >> numfields;
    FieldSurface *newfs = new FieldSurface(fs->data().V, fs->data().F, numfields);
    
    int nfaces = fs->nFaces();
    int nverts = fs->nVerts();
    Eigen::MatrixXd normals(nverts, 3);
    Eigen::MatrixXd vfields(nverts, 3*numfields);
    
    for(int i=0; i<nverts; i++)
    {
        ifs >> normals(i,0) >> normals(i,1) >> normals(i,2);
        for(int j=0; j<3*numfields; j++)
        {
            ifs >> vfields(i, j);
        }
    }
    
    for(int i=0; i<nfaces; i++)
    {
        Eigen::Matrix<double, 3, 2> B = fs->data().Bs[i];
        Eigen::Vector3d n = fs->faceNormal(i);
        
        std::vector<Eigen::Vector3d> alignedfields[3];
        for(int j=0; j<3; j++)
        {
            int vid = fs->data().F(i,j);
            Eigen::Vector3d vertnormal = normals.row(vid).transpose();
            for(int k=0; k<numfields; k++)
            {
                Eigen::Vector3d v = vfields.block(vid, 3*k, 1, 3).transpose();
                Eigen::Vector3d xported = parallelTransport(v, vertnormal, n);
                alignedfields[j].push_back(xported);
            }
        }
        
        alignFrame(alignedfields[1], alignedfields[0], n);
        alignFrame(alignedfields[2], alignedfields[0], n);
        
        Eigen::Matrix2d BTB = B.transpose() * B;
        
        for(int j=0; j<numfields; j++)
        {
            Eigen::Vector3d avv = alignedfields[0][j] + alignedfields[1][j] + alignedfields[2][j];
            avv /= 3.0;
        
            Eigen::Vector2d vint = BTB.inverse() * B.transpose() * avv;
            newfs->vectorFields.segment<2>(newfs->vidx(i,j)) = vint;
        }
    }
    if(ifs)
    {
        delete fs;
        fs = newfs;
    }
    else
    {
        delete newfs;
    }
}

void Weave::deserializeQixingFile(std::ifstream &ifs)
{
    FieldSurface *newfs = new FieldSurface(fs->data().V, fs->data().F, 1);    
    int nfaces = fs->nFaces();
    for(int i=0; i<nfaces; i++)
    {
        char comma;
        Eigen::Vector3d vecs;
        for(int j=0; j<3; j++)
        {
            ifs >> vecs[j];
        }

        Eigen::Matrix<double, 3, 2> B = fs->data().Bs[i];
        Eigen::Vector3d n = fs->faceNormal(i);
        
        Eigen::Matrix2d BTB = B.transpose() * B;
        Eigen::Vector2d vint = BTB.inverse() * B.transpose() * vecs;
        newfs->vectorFields.segment<2>(newfs->vidx(i,0)) = vint;                
    }
    if(ifs)
    {
        delete fs;
        fs = newfs;
    }
    else
    {
        delete newfs;
    }
}

void Weave::deserializeTransportFile(std::ifstream &ifs)
{
    FieldSurface *newfs = new FieldSurface(fs->data().V, fs->data().F, 3);    
    int nfaces = fs->nFaces();
    for(int j = 0; j < 2; j++)
    {
        for(int i=0; i<nfaces; i++)
        {
            double dummy;
            char comma;
            Eigen::Vector3d vecs;
            for(int j=0; j<3; j++)
            {
                ifs >> vecs[j];
            }

            Eigen::Matrix<double, 3, 2> B = fs->data().Bs[i];
            Eigen::Vector3d n = fs->faceNormal(i);
            
            Eigen::Matrix2d BTB = B.transpose() * B;
            Eigen::Vector2d vint = BTB.inverse() * B.transpose() * vecs;
            newfs->vectorFields.segment<2>(newfs->vidx(i,j*2)) = vint;                
        }
    }
    if(ifs)
    {
        delete fs;
        fs = newfs;
    }
    else
    {
        delete newfs;
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
            Eigen::Vector3d n(0,0,0);
            if(f1 != -1)
                n += fs->faceNormal(f1);
            if(f2 != -1)
                n += fs->faceNormal(f2);
            Eigen::Vector3d offset = 0.0001*n;
            cutPts1.row(idx) = fs->data().V.row(fs->data().edgeVerts(cuts[i].path[j].first, 0)) + offset.transpose();
            cutPts2.row(idx) = fs->data().V.row(fs->data().edgeVerts(cuts[i].path[j].first, 1)) + offset.transpose();
            idx++;
        }
    }
}


CoverMesh *Weave::createCover(const std::vector<std::pair<int, int> > &singularities) const
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
    
    /*igl::writeOBJ("debug_single.obj",fs->data().V,fs->data().F);
    std::ofstream debugField("debug.field");
    for (int cId = 0; cId < nCover; cId++)
    {
        for (int fId = 0; fId < nfaces; fId++)
        {
            int field = cId % fs->nFields();
            double sign = (cId < fs->nFields() ? 1.0 : -1.0);
            Eigen::Vector2d vec = sign*fs->v(fId, field).transpose();
            Eigen::Vector3d embvec = fs->data().Bs[fId] * vec;
            debugField << embvec.normalized().transpose() << "\n";
        }
    }*/

    CoverMesh *ret = new CoverMesh(*this->fs, VAug, FAug, oldId2NewId, flattenedField, nCover); 
    // transfer deleted faces
    for(int i=0; i<nfaces; i++)
    {
        for(int j=0; j<nCover; j++)
        {
            ret->fs->setFaceDeleted(j*nfaces + i, fs->isFaceDeleted(i));
        }
    }
    // also remove faces around singularities
    std::map<int, std::set<int> > delvertcovers;
    for(auto it : singularities)
    {
        delvertcovers[it.first].insert(it.second);
    }
    for(int i=0; i<nfaces; i++)
    {
        for(int j=0; j<3; j++)
        {
            int vert = fs->data().F(i,j);
            for(int k=0; k<fs->nFields(); k++)
            {
                if(delvertcovers[vert].count(k))
                {
                    ret->fs->setFaceDeleted(k*nfaces + i, true);
                    ret->fs->setFaceDeleted((k+fs->nFields())*nfaces + i, true);
                }               
            }
        }
    }
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
    int nfields = fs->nFields();
    std::vector<Eigen::MatrixXd> perms;
    Eigen::MatrixXd perm;
    for (int e = 0; e < fs->nEdges(); e++)
    {
        perm = Eigen::MatrixXd::Zero(nCover, nCover);
        for (int j = 0; j < nfields; j++) 
        {
            for (int k = 0; k < nfields; k++)
            {
                if( fs->Ps(e)(j, k) == 1 )
                {
                    perm(j,k) = 1;
                    perm(j+nfields, k+nfields) = 1;
                }
                if( fs->Ps(e)(j, k) == -1 )
                {
                    perm(j,k+nfields) = 1;
                    perm(j+nfields, k) = 1;
                }
            }
        }
        perms.push_back(perm);
    }
    return perms;
}

void Weave::convertToRoSy(int rosyN)
{
    assert(fs->nFields() == 1);        

    int nfaces = fs->nFaces();
    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector2d srcvec = fs->v(i, 0);
        double theta = vectorAngle(*fs, i, srcvec);
        theta *= double(rosyN);

        Eigen::Matrix3d rot = Eigen::AngleAxisd(theta, fs->faceNormal(i)).toRotationMatrix();
        Eigen::Matrix<double, 3, 2> B = fs->data().Bs[i];
        Eigen::Vector3d newvecext = rot * B.col(0);
        newvecext.normalize();
        Eigen::Vector2d newvec = (B.transpose()*B).inverse() * B.transpose() * newvecext;
        int vidx = fs->vidx(i, 0);
        fs->vectorFields.segment<2>(vidx) = newvec;
    }
}


void Weave::transportToRoSy(int rosyN, double factor)
{     
    FieldSurface *newfs = new FieldSurface(fs->data().V, fs->data().F, 1);    
    int nfaces = fs->nFaces();
    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector2d srcvec = fs->v(i, 0);
        Eigen::Vector2d theta;
        theta(0) = vectorAngle(*fs, i, srcvec);
        srcvec = fs->v(i, 2);
        theta(1) = vectorAngle(*fs, i, srcvec);
        theta *= double(rosyN);

        Eigen::Matrix<double, 2, 3> inps;
        Eigen::Matrix<double, 3, 2> B = fs->data().Bs[i];
        
        for (int j = 0; j < 2; j++)
        {
            Eigen::Matrix3d rot = Eigen::AngleAxisd(theta(j), fs->faceNormal(i)).toRotationMatrix();
            Eigen::Vector3d newvecext = rot * B.col(0);
            newvecext.normalize();
            inps.row(j) = newvecext;
        }

        Eigen::Vector3d interp = (1 - factor) * inps.row(0) + factor * inps.row(1);
        Eigen::Vector2d newvec = (B.transpose()*B).inverse() * B.transpose() * interp;
        int vidx = newfs->vidx(i, 0);
        newfs->vectorFields.segment<2>(vidx) = newvec;
    }
    delete fs;
    fs = newfs;
}

Weave *Weave::splitFromRosy(int rosyN)
{
    if (fs->nFields() != 1)
        return NULL;

    bool isodd = (rosyN % 2);
    int m = (isodd ? rosyN : rosyN / 2);

    Weave *result = new Weave(fs->data().V, fs->data().F, m);

    // set vector fields    
    int nfaces = fs->nFaces();
    std::vector<bool> visited(nfaces, false);
    // bread-first search of faces
    for (int i = 0; i < nfaces; i++)
    {
        if (visited[i])
            continue;
        
        struct Visit
        {
            int from;
            int to;
        };

        std::deque<Visit> q;
        q.push_back(Visit{ -1,i });
        while (!q.empty())
        {
            Visit vis = q.front();
            q.pop_front();
            if (visited[vis.to])
                continue;
            visited[vis.to] = true;

            std::vector<Eigen::Vector2d> rv;
            repVecToRoSy(*fs, vis.to, fs->v(vis.to, 0), rv, rosyN);
            
            if (vis.from == -1)
            {
                for (int j = 0; j < m; j++)
                {
                    int vidx = result->fs->vidx(vis.to, j);
                    result->fs->vectorFields.segment<2>(vidx) = rv[j];
                }
            }
            else
            {
                // find optimal cyclic permutation

                // vector 0 on neighbor face
                Eigen::Vector2d v0 = result->fs->v(vis.from, 0);
                // tranport to current face
                int edge = -1;
                int fromside = -1;
                for (int j = 0; j < 3; j++)
                {
                    int eid = result->fs->data().faceEdges(vis.to, j);
                    int opp = 0;
                    if(result->fs->data().E(eid, opp) == vis.to)
                        opp = 1;
                    if (result->fs->data().E(eid, opp) == vis.from)
                    {
                        edge = eid;
                        fromside = opp;
                    }
                }
                assert(edge != -1);
                assert(result->fs->data().E(edge, fromside) == vis.from);
                assert(result->fs->data().E(edge, 1-fromside) == vis.to);
                Eigen::Vector2d xportv0 = result->fs->data().Ts.block<2, 2>(2 * edge, 2 * fromside) * v0;
                int bestidx = 0;
                double bestdot = -1.0;
                for (int j = 0; j < rosyN; j++)
                {
                    Eigen::Vector3d vec1 = result->fs->data().Bs[vis.to] * xportv0;
                    vec1.normalize();
                    Eigen::Vector3d vec2 = result->fs->data().Bs[vis.to] * rv[j];
                    vec2.normalize();
                    double curdot = (vec1).dot(vec2);
                    if (curdot > bestdot)
                    {
                        bestidx = j;
                        bestdot = curdot;
                    }
                }

                // set vectors
                for (int j = 0; j < m; j++)
                {
                    int vidx = result->fs->vidx(vis.to, j);
                    result->fs->vectorFields.segment<2>(vidx) = rv[(j + bestidx) % rosyN];

                }
            }
            
            // queue neighbors
            for (int j = 0; j < 3; j++)
            {
                int nb = result->fs->data().faceNeighbors(vis.to, j);
                if (nb != -1 && !visited[nb])
                    q.push_back(Visit{ vis.to, nb });
            }
        }        
    }

    //reassignAllPermutations(*result);

    return result;
}
