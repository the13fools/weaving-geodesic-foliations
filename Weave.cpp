#include "Weave.h"
#include <igl/read_triangle_mesh.h>
#include <map>
#include <Eigen/Dense>
#include "Colors.h"
#include <deque>
#include <algorithm>

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
    std::ofstream ofs(filename);
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
		 << edgeVerts(i, 0) << " "
		 << edgeVerts(i, 1) <<  std::endl;
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
 

        for (int j = 0; j < nfields; j++)
        {
            for (int k = 0; k < nfields; k++)
            {
                ofs_mat << perm(j, k) << " ";
            }
            ofs_mat << std::endl;
        }
//        ofs << std::endl;
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
        std::cerr << "Vector field doesn't match mesh!" << std::endl;
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
