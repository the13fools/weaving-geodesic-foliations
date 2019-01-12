#include "Surface.h"
#include <map>
#include <queue>
#include <Eigen/Dense>

#include <iostream>

Surface::Surface(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
    data_.V = V;
    data_.F = F;
    buildConnectivityStructures();
    buildGeometricStructures();
}


void Surface::buildConnectivityStructures()
{
    std::map<std::pair<int, int>, Eigen::Vector2i > edgemap;
    int nfaces = nFaces();
    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int nextj = (j + 1) % 3;
            int v1 = data_.F(i, j);
            int v2 = data_.F(i, nextj);
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

    data_.E.resize(edgemap.size(), 2);
    data_.faceEdges.resize(nfaces, 3);
    data_.edgeVerts.resize(edgemap.size(), 2);
    int idx = 0;
    for (std::map<std::pair<int, int>, Eigen::Vector2i >::iterator it = edgemap.begin(); it != edgemap.end(); ++it)
    {
        data_.E(idx, 0) = it->second[0];
        data_.E(idx, 1) = it->second[1];
        data_.edgeVerts(idx, 0) = it->first.first;
        data_.edgeVerts(idx, 1) = it->first.second;
        idx++;
    }
    data_.faceEdges.resize(nfaces, 3);
    data_.faceEdges.setConstant(-1);
    data_.faceNeighbors.resize(nfaces, 3);
    data_.faceNeighbors.setConstant(-1);
    data_.faceWings.resize(nfaces, 3);
    data_.faceWings.setConstant(-1);
    int nedges = nEdges();
    for (int edge = 0; edge < nedges; edge++)
    {
        for(int side = 0; side<2; side++)
        {
            if(data_.E(edge,side) == -1)
                continue;
            for(int j=0; j<3; j++)
                if(data_.F(data_.E(edge,side), j) != data_.edgeVerts(edge,0) && data_.F(data_.E(edge,side), j) != data_.edgeVerts(edge,1))
                    data_.faceEdges(data_.E(edge, side), j) = edge;           
        }
        if(data_.E(edge,0) == -1 || data_.E(edge,1) == -1)
            continue;
        Eigen::Vector3i face1 = data_.F.row(data_.E(edge, 0));
        Eigen::Vector3i face2 = data_.F.row(data_.E(edge, 1));
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
        data_.faceNeighbors(data_.E(edge, 0), idx1) = data_.E(edge, 1);
        data_.faceNeighbors(data_.E(edge, 1), idx2) = data_.E(edge, 0);
        data_.faceWings(data_.E(edge, 0), idx1) = face2[idx2];
        data_.faceWings(data_.E(edge, 1), idx2) = face1[idx1];
    }

    data_.vertEdges.resize(data_.V.rows());
    for(int i=0; i<nedges; i++)
    {
        data_.vertEdges[data_.edgeVerts(i,0)].push_back(i);
        data_.vertEdges[data_.edgeVerts(i,1)].push_back(i);
    }
}

Eigen::Vector3d Surface::faceNormal(int face) const
{
    Eigen::Vector3d e1 = (data_.V.row(data_.F(face, 1)) - data_.V.row(data_.F(face, 0)));
    Eigen::Vector3d e2 = (data_.V.row(data_.F(face, 2)) - data_.V.row(data_.F(face, 0)));
    Eigen::Vector3d result = e1.cross(e2);
    result /= result.norm();
    return result;
}


double Surface::faceArea(int face) const
{
    Eigen::Vector3d e1 = (data_.V.row(data_.F(face, 1)) - data_.V.row(data_.F(face, 0)));
    Eigen::Vector3d e2 = (data_.V.row(data_.F(face, 2)) - data_.V.row(data_.F(face, 0)));
    Eigen::Vector3d result = e1.cross(e2);
    return 0.5 * result.norm();
}

void Surface::buildGeometricStructures()
{
    // compute barycentric matrices and Js
    int nfaces = nFaces();
    data_.Bs.resize(nfaces);
    data_.Js.resize(2 * nfaces, 2);

    data_.averageEdgeLength = 0;
    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector3d v0 = data_.V.row(data_.F(i, 0));
        Eigen::Vector3d v1 = data_.V.row(data_.F(i, 1));
        Eigen::Vector3d v2 = data_.V.row(data_.F(i, 2));
        data_.Bs[i].col(0) = v1 - v0;
        data_.Bs[i].col(1) = v2 - v0;

        data_.averageEdgeLength += (v1 - v0).norm();
        data_.averageEdgeLength += (v2 - v1).norm();
        data_.averageEdgeLength += (v0 - v2).norm();

        Eigen::Vector3d n = (v1 - v0).cross(v2 - v0);
        n /= n.norm();

    //    data_.Bs[i].col(1) = (v1 - v0).cross(n); // switch to orthogonal local coordinate system.

        Eigen::Matrix2d BTB = data_.Bs[i].transpose() * data_.Bs[i];
        Eigen::Matrix<double, 3, 2> ncrossB;
        ncrossB.col(0) = n.cross(v1 - v0);
        ncrossB.col(1) = n.cross(v2 - v0);
        data_.Js.block<2, 2>(2 * i, 0) = BTB.inverse() * data_.Bs[i].transpose() * ncrossB;
    }

    data_.averageEdgeLength /= 3.0 * nfaces;

    // compute cDiffs and transition matrices
    int nedges = nEdges();
    data_.cDiffs.resize(2 * nedges, 2);
    data_.Ts.resize(2 * nedges, 4);
    data_.Ts_rosy.resize(2 * nedges, 4);
    for (int edgeidx = 0; edgeidx < nedges; edgeidx++)
    {        
        //collect neighboring face indices
        int face1 = data_.E(edgeidx, 0);
        int face2 = data_.E(edgeidx, 1);

        if(face1 == -1 || face2 == -1)
            continue;

        int v1 = data_.edgeVerts(edgeidx, 0);
        int v2 = data_.edgeVerts(edgeidx, 1);

        Eigen::Vector3d n1 = faceNormal(face1);
        Eigen::Vector3d n2 = faceNormal(face2);
        // collect: (1) the midpoint of the common edge, (2) unit vector in direction of common edge,
        // (3) the face normals, (4) the centroid of neighboring faces

        Eigen::Vector3d midpt = 0.5 * (data_.V.row(v1) + data_.V.row(v2));
        Eigen::Vector3d commone = data_.V.row(v2) - data_.V.row(v1);
        commone /= commone.norm();
        Eigen::Vector3d centroids[2];
        centroids[0].setZero();
        centroids[1].setZero();

        for (int i = 0; i < 3; i++)
        {
            centroids[0] += data_.V.row(data_.F(face1, i));
            centroids[1] += data_.V.row(data_.F(face2, i));
        }

        centroids[0] /= 3.0;
        centroids[1] /= 3.0;

        //rotate each centroid into the plane of the opposite triangle and compute ci minus c

        Eigen::Vector3d t1 = n1.cross(commone);
        Eigen::Vector3d t2 = n2.cross(commone);
        Eigen::Vector3d diff2 = centroids[1] - midpt;
        double alpha = commone.dot(diff2);
        double beta = t2.dot(diff2);

        Eigen::Matrix2d BTB1 = data_.Bs[face1].transpose() * data_.Bs[face1];
        Eigen::Matrix2d BTB2 = data_.Bs[face2].transpose() * data_.Bs[face2];

        data_.cDiffs.row(2 * edgeidx) = BTB1.inverse() * data_.Bs[face1].transpose() * (midpt + alpha * commone + beta * t1 - centroids[0]);
        Eigen::Vector3d diff1 = centroids[0] - midpt;
        alpha = commone.dot(diff1);
        beta = t1.dot(diff1);
        data_.cDiffs.row(2 * edgeidx + 1) = BTB2.inverse() * data_.Bs[face2].transpose() * (midpt + alpha*commone + beta * t2 - centroids[1]);

        Eigen::Vector3d e1 = data_.V.row(data_.F(face1, 1)) - data_.V.row(data_.F(face1, 0));
        Eigen::Vector3d e2 = data_.V.row(data_.F(face1, 2)) - data_.V.row(data_.F(face1, 0));

        double alpha1 = commone.dot(e1);
        double beta1 = t1.dot(e1);
        Eigen::Vector3d newe1 = alpha1*commone + beta1 * t2;
        data_.Ts.block<2, 1>(2 * edgeidx, 0) = BTB2.inverse() * data_.Bs[face2].transpose() * newe1;

        double alpha2 = commone.dot(e2);
        double beta2 = t1.dot(e2);
        Eigen::Vector3d newe2 = alpha2*commone + beta2*t2;
        data_.Ts.block<2, 1>(2 * edgeidx, 1) = BTB2.inverse() * data_.Bs[face2].transpose() * newe2;

        e1 = data_.V.row(data_.F(face2, 1)) - data_.V.row(data_.F(face2, 0));
        e2 = data_.V.row(data_.F(face2, 2)) - data_.V.row(data_.F(face2, 0));

        alpha1 = commone.dot(e1);
        beta1 = t2.dot(e1);
        newe1 = alpha1 * commone + beta1 * t1;
        data_.Ts.block<2, 1>(2 * edgeidx, 2) = BTB1.inverse() * data_.Bs[face1].transpose() * newe1;

        alpha2 = commone.dot(e2);
        beta2 = t2.dot(e2);
        newe2 = alpha2*commone + beta2*t1;
        data_.Ts.block<2, 1>(2 * edgeidx, 3) = BTB1.inverse() * data_.Bs[face1].transpose() * newe2;


        Eigen::Vector2d vec(1, 0); 
        // use the computed transition matricies to compute the tranition angle in local coordinates
        Eigen::Vector3d basis1 = data_.V.row(data_.F(face1, 1)) - data_.V.row(data_.F(face1, 0));
        Eigen::Vector3d basis2 = data_.V.row(data_.F(face2, 1)) - data_.V.row(data_.F(face2, 0));
        Eigen::Vector3d b1_f2 =  data_.Bs[face2] * data_.Ts.block<2, 2>(2 * edgeidx, 0) * vec;
        Eigen::Vector3d b2_f1 =  data_.Bs[face1] * data_.Ts.block<2, 2>(2 * edgeidx, 2) * vec;
        basis1.normalize();
        basis2.normalize();
        b1_f2.normalize();
        b2_f1.normalize();

        Eigen::Matrix3d source1;
        source1.col(0) = b1_f2;
        source1.col(1) = n2.cross(b1_f2);
        source1.col(2) = n2;

        Eigen::Matrix3d target1;
        target1.col(0) = basis2;
        target1.col(1) = n2.cross(basis2);
        target1.col(2) = n2;

        Eigen::Matrix3d source2;
        source2.col(0) = b2_f1;
        source2.col(1) = n1.cross(b2_f1);
        source2.col(2) = n1;

        Eigen::Matrix3d target2;
        target2.col(0) = basis1;
        target2.col(1) = n1.cross(basis1);
        target2.col(2) = n1;

        // note the index change because we apply the rotation on the target face
        Eigen::Matrix3d R1 = target1 * source1.inverse();
        Eigen::Matrix3d R2 = target2 * source2.inverse();

     //   std::cout << R1.determinant() << " " << R2.determinant() << std::endl;

        data_.Ts_rosy.block<2, 2>(2 * edgeidx, 0) = BTB2.inverse() * data_.Bs[face2].transpose() * R1 * data_.Bs[face2];
        data_.Ts_rosy.block<2, 2>(2 * edgeidx, 2) = BTB1.inverse() * data_.Bs[face1].transpose() * R2 * data_.Bs[face1];      
    }
}

int Surface::numInteriorEdges() const
{
    int nedges = nEdges();
    int result = 0;
    for(int i=0; i<nedges; i++)
    {
        if(data_.E(i,0) != -1 && data_.E(i,1) != -1)
            result++;
    }
    return result;
}


void Surface::shortestPath(int startVert, int endVert, std::vector<std::pair<int, int> > &path) const
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

        int nbs = (int)data_.vertEdges[cur.next].size();
        for(int i=0; i<nbs; i++)
        {
            int e = data_.vertEdges[cur.next][i];
            int orient = (data_.edgeVerts(e, 0) == cur.next) ? 0 : 1;
            int next = data_.edgeVerts(e, 1-orient);
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
