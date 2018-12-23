#include "BommesFieldIntegration.h"
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <deque>
#include <igl/remove_unreferenced.h>
#include <igl/cotmatrix.h>
#include <igl/writeOBJ.h>
#include "CoMISoWrapper.h"

static void findCuts(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
    std::vector<std::vector<int> > &cuts)
{
    cuts.clear();

    int nfaces = F.rows();
    int nverts = V.rows();

    if (nfaces == 0)
        return;

    std::map<std::pair<int, int>, std::vector<int> > edges;
    // build edges

    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int v0 = F(i, j);
            int v1 = F(i, (j + 1) % 3);
            std::pair<int, int> e;
            e.first = std::min(v0, v1);
            e.second = std::max(v0, v1);
            edges[e].push_back(i);
        }
    }

    int nedges = edges.size();
    Eigen::MatrixXi edgeVerts(nedges,2);
    Eigen::MatrixXi edgeFaces(nedges,2);
    Eigen::MatrixXi faceEdges(nfaces, 3);
    std::set<int> boundaryEdges;
    std::map<std::pair<int, int>, int> edgeidx;
    int idx = 0;
    for (auto it : edges)
    {
        edgeidx[it.first] = idx;
        edgeVerts(idx, 0) = it.first.first;
        edgeVerts(idx, 1) = it.first.second;
        edgeFaces(idx, 0) = it.second[0];
        if (it.second.size() > 1)
        {
            edgeFaces(idx, 1) = it.second[1];
        }
        else
        {
            edgeFaces(idx, 1) = -1;
            boundaryEdges.insert(idx);
        }
        idx++;
    }
    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int v0 = F(i, j);
            int v1 = F(i, (j + 1) % 3);
            std::pair<int, int> e;
            e.first = std::min(v0, v1);
            e.second = std::max(v0, v1);
            faceEdges(i, j) = edgeidx[e];
        }
    }

    bool *deleted = new bool[nfaces];
    for (int i = 0; i < nfaces; i++)
        deleted[i] = false;

    std::set<int> deletededges;

    // loop over faces
    for (int face = 0; face < nfaces; face++)
    {
        // stop at first undeleted face
        if (deleted[face])
            continue;
        deleted[face] = true;
        std::deque<int> processEdges;
        for (int i = 0; i < 3; i++)
        {
            int e = faceEdges(face, i);
            if (boundaryEdges.count(e))
                continue;
            int ndeleted = 0;
            if (deleted[edgeFaces(e, 0)])
                ndeleted++;
            if (deleted[edgeFaces(e, 1)])
                ndeleted++;
            if (ndeleted == 1)
                processEdges.push_back(e);
        }
        // delete all faces adjacent to edges with exactly one adjacent face
        while (!processEdges.empty())
        {
            int nexte = processEdges.front();
            processEdges.pop_front();
            int todelete = -1;
            if (!deleted[edgeFaces(nexte, 0)])
                todelete = edgeFaces(nexte, 0);
            if (!deleted[edgeFaces(nexte, 1)])
                todelete = edgeFaces(nexte, 1);
            if (todelete != -1)
            {
                deletededges.insert(nexte);
                deleted[todelete] = true;
                for (int i = 0; i < 3; i++)
                {
                    int e = faceEdges(todelete, i);
                    if (boundaryEdges.count(e))
                        continue;
                    int ndeleted = 0;
                    if (deleted[edgeFaces(e, 0)])
                        ndeleted++;
                    if (deleted[edgeFaces(e, 1)])
                        ndeleted++;
                    if (ndeleted == 1)
                        processEdges.push_back(e);
                }
            }
        }
    }
    delete[] deleted;

    // accumulated non-deleted edges
    std::vector<int> leftedges;
    for (int i = 0; i < nedges; i++)
    {
        if (!deletededges.count(i))
            leftedges.push_back(i);
    }

    deletededges.clear();
    // prune spines
    std::map<int, std::vector<int> > spinevertedges;
    for (int i : leftedges)
    {
        spinevertedges[edgeVerts(i, 0)].push_back(i);
        spinevertedges[edgeVerts(i, 1)].push_back(i);
    }

    std::deque<int> vertsProcess;
    std::map<int, int> spinevertnbs;
    for (auto it : spinevertedges)
    {
        spinevertnbs[it.first] = it.second.size();
        if (it.second.size() == 1)
            vertsProcess.push_back(it.first);
    }
    while (!vertsProcess.empty())
    {
        int vert = vertsProcess.front();
        vertsProcess.pop_front();
        for (int e : spinevertedges[vert])
        {
            if (!deletededges.count(e))
            {
                deletededges.insert(e);
                for (int j = 0; j < 2; j++)
                {
                    spinevertnbs[edgeVerts(e, j)]--;
                    if (spinevertnbs[edgeVerts(e, j)] == 1)
                    {
                        vertsProcess.push_back(edgeVerts(e, j));
                    }
                }
            }
        }
    }
    std::vector<int> loopedges;
    for (int i : leftedges)
        if (!deletededges.count(i))
            loopedges.push_back(i);

    int nloopedges = loopedges.size();
    if (nloopedges == 0)
        return;

    std::map<int, std::vector<int> > loopvertedges;
    for (int e : loopedges)
    {
        loopvertedges[edgeVerts(e, 0)].push_back(e);
        loopvertedges[edgeVerts(e, 1)].push_back(e);
    }

    std::set<int> usededges;
    for (int e : loopedges)
    {
        // make a cycle or chain starting from this edge
        while (!usededges.count(e))
        {
            std::vector<int> cycleverts;
            std::vector<int> cycleedges;
            cycleverts.push_back(edgeVerts(e, 0));
            cycleverts.push_back(edgeVerts(e, 1));
            cycleedges.push_back(e);

            std::map<int, int> cycleidx;
            cycleidx[cycleverts[0]] = 0;
            cycleidx[cycleverts[1]] = 1;

            int curvert = edgeVerts(e, 1);
            int cure = e;
            bool foundcycle = false;
            while (curvert != -1 && !foundcycle)
            {
                int nextvert = -1;
                int nexte = -1;
                for (int cande : loopvertedges[curvert])
                {
                    if (!usededges.count(cande) && cande != cure)
                    {
                        int vidx = 0;
                        if (curvert == edgeVerts(cande, vidx))
                            vidx = 1;
                        nextvert = edgeVerts(cande, vidx);
                        nexte = cande;
                        break;
                    }
                }
                if (nextvert != -1)
                {
                    auto it = cycleidx.find(nextvert);
                    if (it != cycleidx.end())
                    {
                        // we've hit outselves
                        std::vector<int> cut;
                        for (int i = it->second; i < cycleverts.size(); i++)
                        {
                            cut.push_back(cycleverts[i]);
                        }
                        cut.push_back(nextvert);
                        cuts.push_back(cut);
                        for (int i = it->second; i < cycleedges.size(); i++)
                        {
                            usededges.insert(cycleedges[i]);
                        }
                        usededges.insert(nexte);
                        foundcycle = true;
                    }
                    else
                    {
                        cycleidx[nextvert] = cycleverts.size();
                        cycleverts.push_back(nextvert);
                        cycleedges.push_back(nexte);                        
                    }
                }                
                curvert = nextvert;
                cure = nexte;
            }
            if (!foundcycle)
            {
                // we've hit a dead end. reverse and try the other direction
                std::reverse(cycleverts.begin(), cycleverts.end());
                std::reverse(cycleedges.begin(), cycleedges.end());
                cycleidx.clear();
                for (int i = 0; i < cycleverts.size(); i++)
                {
                    cycleidx[cycleverts[i]] = i;
                }
                curvert = cycleverts.back();
                cure = cycleedges.back();
                while (curvert != -1 && !foundcycle)
                {
                    int nextvert = -1;
                    int nexte = -1;
                    for (int cande : loopvertedges[curvert])
                    {
                        if (!usededges.count(cande) && cande != cure)
                        {
                            int vidx = 0;
                            if (curvert == edgeVerts(cande, vidx))
                                vidx = 1;
                            nextvert = edgeVerts(cande, vidx);
                            nexte = cande;
                            break;
                        }
                    }
                    if (nextvert != -1)
                    {
                        auto it = cycleidx.find(nextvert);
                        if (it != cycleidx.end())
                        {
                            // we've hit outselves
                            std::vector<int> cut;
                            for (int i = it->second; i < cycleverts.size(); i++)
                            {
                                cut.push_back(cycleverts[i]);
                            }
                            cut.push_back(nextvert);
                            cuts.push_back(cut);
                            for (int i = it->second; i < cycleedges.size(); i++)
                            {
                                usededges.insert(cycleedges[i]);
                            }
                            usededges.insert(nexte);
                            foundcycle = true;
                        }
                        else
                        {
                            cycleidx[nextvert] = cycleverts.size();
                            cycleverts.push_back(nextvert);
                            cycleedges.push_back(nexte);                        
                        }
                    }                
                    curvert = nextvert;
                    cure = nexte;
                }
                if (!foundcycle)
                {
                    // we've found a chain
                    std::vector<int> cut;
                    for (int i = 0; i < cycleverts.size(); i++)
                    {
                        cut.push_back(cycleverts[i]);
                    }
                    cuts.push_back(cut);
                    for (int i = 0; i < cycleedges.size(); i++)
                    {
                        usededges.insert(cycleedges[i]);
                    }                    
                }
            }
        }
    }
}


static void cutMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
    // list of cuts, each of which is a list (in order) of vertex indices of one cut.
    // Cuts can be closed loops (in which case the last vertex index should equal the
    // first) or open (in which case the two endpoint vertices should be distinct).
    // Multiple cuts can cross but there may be strange behavior if cuts share endpoint
    // vertices, or are non-edge-disjoint.
    const std::vector<std::vector<int> > &cuts,
    // new vertices and faces
    // **DO NOT ALIAS V OR F!**
    Eigen::MatrixXd &newV,
    Eigen::MatrixXi &newF
)
{
    int ncuts = (int)cuts.size();

    // junction vertices that lie on multiple cuts
    std::set<int> junctions;
    std::set<int> seenverts;
    for (int i = 0; i < ncuts; i++)
    {
        std::set<int> seenincut;
        for (int j = 0; j < cuts[i].size(); j++)
        {
            if (seenverts.count(cuts[i][j]))
                junctions.insert(cuts[i][j]);
            seenincut.insert(cuts[i][j]);
        }
        for (int v : seenincut)
            seenverts.insert(v);
    }

    // "interior" cut vertices: vertices that are part of a cut but not a cut endpoint
    // or junction vertex
    std::vector<std::set<int> > cutints;
    cutints.resize(ncuts);
    for (int i = 0; i < ncuts; i++)
    {
        if (cuts[i].empty())
            continue;
        if (cuts[i].front() == cuts[i].back())
        {
            // closed loop
            for (int v : cuts[i])
            {
                if(!junctions.count(v))
                    cutints[i].insert(v);                
            }
        }
        else
        {
            // open cut
            for (int j = 1; j < cuts[i].size() - 1; j++)
            {
                if(!junctions.count(cuts[i][j]))
                    cutints[i].insert(cuts[i][j]);                
            }
        }
    }

    struct edge
    {
        std::pair<int, int> verts;
        edge(int v1, int v2)
        {
            verts.first = std::min(v1, v2);
            verts.second = std::max(v1, v2);
        }

        bool operator<(const edge &other) const
        {
            return verts < other.verts;
        }
    };

    // maps each edge to incident triangles
    std::map<edge, std::vector<int> > edgeTriangles;
    for (int i = 0; i < F.rows(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            edge e(F(i, j), F(i, (j + 1) % 3));
            edgeTriangles[e].push_back(i);
        }
    }

    // have we visited this face yet?
    bool *visited = new bool[F.rows()];

    // edges that form part of a cut
    std::set<edge> forbidden;
    for (int i = 0; i < ncuts; i++)
    {        
        for (int j = 0; j < (int)cuts[i].size() - 1; j++)
        {
            // works for both open and closed curves
            edge e(cuts[i][j], cuts[i][j + 1]);
            forbidden.insert(e);
        }
    }

    // connected components of faces adjacent to the cuts
    std::vector<std::vector<std::vector<int> > > components;
    components.resize(ncuts);

    // for each cut
    for (int cut = 0; cut < ncuts; cut++)
    {
        for (int i = 0; i < (int)F.rows(); i++)
            visited[i] = false;

        // find a face we haven't visited yet
        for (int i = 0; i < F.rows(); i++)
        {
            if (visited[i]) continue;
            bool found = false;
            for (int j = 0; j < 3; j++)
            {
                if (cutints[cut].count(F(i, j)))
                {
                    found = true;
                }
            }

            if (found)
            {
                // run a BFS along the cut edges, accumulating one connected component
                // cross only edges that contain a vertex in cutints[cut], but are not forbidden
                std::deque<int> q;
                std::vector<int> component;
                q.push_back(i);
                while (!q.empty())
                {
                    int next = q.front();
                    q.pop_front();
                    if (visited[next])
                        continue;
                    visited[next] = true;
                    component.push_back(next);
                    for (int j = 0; j < 3; j++)
                    {
                        int v1 = F(next, j);
                        int v2 = F(next, (j + 1) % 3);
                        edge e(v1, v2);
                        if (cutints[cut].count(v1) == 0 && cutints[cut].count(v2) == 0)
                        {
                            continue;
                        }
                        if (forbidden.count(e))
                            continue;
                        for (int nb : edgeTriangles[e])
                        {
                            if (!visited[nb])
                            {
                                q.push_back(nb);
                            }
                        }
                    } // end BFS
                }
                components[cut].push_back(component);
            } // end if found
        } // end loop over all faces
    } // end loop over cuts

    std::map<int, std::vector<std::vector<int> > > junctcomponents;

    // for each junction
    for (int junc : junctions)
    {
        for (int i = 0; i < (int)F.rows(); i++)
            visited[i] = false;

        // find a face we haven't visited yet
        for (int i = 0; i < F.rows(); i++)
        {
            if (visited[i]) continue;
            bool found = false;
            for (int j = 0; j < 3; j++)
            {
                if (junc == F(i, j))
                {
                    found = true;
                }
            }

            if (found)
            {
                // run a BFS along the cut edges, accumulating one connected component
                // cross only edges that contain the junction, but are not forbidden
                std::deque<int> q;
                std::vector<int> component;
                q.push_back(i);
                while (!q.empty())
                {
                    int next = q.front();
                    q.pop_front();
                    if (visited[next])
                        continue;
                    visited[next] = true;
                    component.push_back(next);
                    for (int j = 0; j < 3; j++)
                    {
                        int v1 = F(next, j);
                        int v2 = F(next, (j + 1) % 3);
                        edge e(v1, v2);
                        if (v1 != junc && v2 != junc)
                        {
                            continue;
                        }
                        if (forbidden.count(e))
                            continue;
                        for (int nb : edgeTriangles[e])
                        {
                            if (!visited[nb])
                            {
                                q.push_back(nb);
                            }
                        }
                    } // end BFS
                }
                junctcomponents[junc].push_back(component);
            } // end if found
        } // end loop over all faces
    } // end loop over cuts

    int vertstoadd = 0;
    // create a copy of each vertex for each component of each cut
    for (int i = 0; i < ncuts; i++)
    {
        vertstoadd += components[i].size() * cutints[i].size();
    }
    // create a copy of each junction point for each component of each junction
    for (int v : junctions)
    {
        vertstoadd += junctcomponents[v].size();
    }

    Eigen::MatrixXd augV(V.rows() + vertstoadd, 3);
    for (int i = 0; i < V.rows(); i++)
    {
        augV.row(i) = V.row(i);
    }

    // create new faces
    Eigen::MatrixXi augF = F;

    // duplicate vertices and reindex faces
    int idx = V.rows();
    for (int cut = 0; cut < ncuts; cut++)
    {
        for (int i = 0; i < components[cut].size(); i++)
        {
            // duplicate vertices
            std::map<int, int> idxmap;
            for (int v : cutints[cut])
            {
                idxmap[v] = idx;
                augV.row(idx) = V.row(v);
                idx++;
            }
            for (int f : components[cut][i])
            {
                for (int j = 0; j < 3; j++)
                {
                    int v = augF(f, j);
                    if (cutints[cut].count(v))
                        augF(f, j) = idxmap[v];
                }
            }
        }        
    }

    for (int junc : junctions)
    {
        for (int i = 0; i < junctcomponents[junc].size(); i++)
        {

            augV.row(idx) = V.row(junc);
            int newidx = idx;
            idx++;

            for (int f : junctcomponents[junc][i])
            {
                for (int j = 0; j < 3; j++)
                {
                    int v = augF(f, j);
                    if (v == junc)
                        augF(f, j) = newidx;
                }
            }
        }
    }

    assert(idx == augV.rows());

    // some duplicated vertices will not have been used
    Eigen::MatrixXi I;
    igl::remove_unreferenced(augV, augF, newV, newF, I);
    delete[] visited;
}

void BommesFieldIntegration::integrateOneComponent(const Surface &surf, const Eigen::MatrixXd &v, Eigen::VectorXd &theta)
{
    int nverts = surf.nVerts();
    theta.resize(nverts);
    theta.setZero();
    std::vector<std::vector<int> > cuts;
    findCuts(surf.data().V, surf.data().F, cuts);
    std::cout << "Used " << cuts.size() << " cuts" << std::endl;

    std::map<std::pair<int, int>, int> cutpairs;

    // convert cuts to edge indices
    for (int i = 0; i < cuts.size(); i++)
    {        
        int len = cuts[i].size();
        for (int j = 0; j < len - 1; j++)
        {
            cutpairs[std::pair<int, int>(cuts[i][j], cuts[i][j + 1])] = i;
        }
    }

    Eigen::MatrixXd newV;
    Eigen::MatrixXi newF;

    cutMesh(surf.data().V, surf.data().F, cuts, newV, newF);
    
    Surface s(newV, newF);
    int newverts = s.nVerts();
    int intdofs = cuts.size();

    std::vector<Eigen::Triplet<double> > Ccoeffs;
    int row = 0;
    for (int i = 0; i < surf.nEdges(); i++)
    {
        int f1 = surf.data().E(i, 0);
        int f2 = surf.data().E(i, 1);
        if (f1 == -1 || f2 == -1)
            continue;
        int v1 = surf.data().edgeVerts(i, 0);
        int v2 = surf.data().edgeVerts(i, 1);
        double sign = 1.0;
        auto it = cutpairs.find(std::pair<int, int>(v1, v2));
        if (it == cutpairs.end())
        {
            it = cutpairs.find(std::pair<int, int>(v2, v1));
            sign = -1.0;
        }
        if (it == cutpairs.end())
            continue;

        int newv1 = -1;
        int newv2 = -1;
        int neww1 = -1;
        int neww2 = -1;
        for (int j = 0; j < 3; j++)
        {
            if (surf.data().F(f1, j) == v1)
                newv1 = newF(f1, j);
            if (surf.data().F(f2, j) == v1)
                newv2 = newF(f2, j);
            if (surf.data().F(f1, j) == v2)
                neww1 = newF(f1, j);
            if (surf.data().F(f2, j) == v2)
                neww2 = newF(f2, j);
        }

        Ccoeffs.push_back(Eigen::Triplet<double>(row, newv1, sign));
        Ccoeffs.push_back(Eigen::Triplet<double>(row, newv2, -sign));
        Ccoeffs.push_back(Eigen::Triplet<double>(row, newverts + it->second, 1.0));
        Ccoeffs.push_back(Eigen::Triplet<double>(row + 1, neww1, sign));
        Ccoeffs.push_back(Eigen::Triplet<double>(row + 1, neww2, -sign));
        Ccoeffs.push_back(Eigen::Triplet<double>(row + 1, newverts + it->second, 1.0));
        row += 2;
    }
    Eigen::SparseMatrix<double> C(row, newverts + intdofs + 1);
    C.setFromTriplets(Ccoeffs.begin(), Ccoeffs.end());
    //igl::writeOBJ("test.obj", newV, newF);

    std::vector<Eigen::Triplet<double> > Acoeffs;
    int nfaces = s.nFaces();
    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Acoeffs.push_back(Eigen::Triplet<double>(3 * i + j, s.data().F(i, j), 1.0));
        }
    }
    Eigen::SparseMatrix<double> A(3 * nfaces, newverts + intdofs);
    A.setFromTriplets(Acoeffs.begin(), Acoeffs.end());

    std::vector<Eigen::Triplet<double> > Dcoeffs;
    for (int i = 0; i < nfaces; i++)
    {
        Dcoeffs.push_back(Eigen::Triplet<double>(2 * i, 3 * i + 1, 1.0));
        Dcoeffs.push_back(Eigen::Triplet<double>(2 * i, 3 * i, -1.0));
        Dcoeffs.push_back(Eigen::Triplet<double>(2 * i + 1, 3 * i + 2, 1.0));
        Dcoeffs.push_back(Eigen::Triplet<double>(2 * i + 1, 3 * i, -1.0));
    }
    Eigen::SparseMatrix<double> D(2 * nfaces, 3 * nfaces);
    D.setFromTriplets(Dcoeffs.begin(), Dcoeffs.end());

    std::vector<Eigen::Triplet<double> > Mcoeffs;
    double totarea = 0;
    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Matrix<double, 3, 2> B = s.data().Bs[i];
        Eigen::Vector3d vec = B*v.row(i).transpose();
        vec.normalize();
        Eigen::Vector3d n = s.faceNormal(i);
        double area = s.faceArea(i);
        totarea += area;
        Eigen::Vector3d vecperp = n.cross(vec);
        Eigen::Matrix3d metric = vec * vec.transpose() + aniso_*vecperp * vecperp.transpose();
        Eigen::Matrix2d BTB = area * B.transpose()*metric *B;
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
            {
                Mcoeffs.push_back(Eigen::Triplet<double>(2 * i + j, 2 * i + k, BTB(j, k)));
            }
    }
    Eigen::SparseMatrix<double> M(2 * nfaces, 2 * nfaces);
    M.setFromTriplets(Mcoeffs.begin(), Mcoeffs.end());
    Eigen::VectorXd projvf(2 * nfaces);
    double scale = surf.data().averageEdgeLength * globalScale_;
    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector3d vec = surf.data().Bs[i] * v.row(i).transpose();
        double fac = scale / vec.norm();
        projvf.segment<2>(2 * i) = fac * surf.data().Js.block<2,2>(2*i,0) * v.row(i).transpose();
    }

    Eigen::VectorXd rhs = A.transpose() * D.transpose() * M * projvf;
    Eigen::SparseMatrix<double> Mat = A.transpose()  * D.transpose() * M * D * A;

    Eigen::SparseMatrix<double> L;
    igl::cotmatrix(s.data().V, s.data().F, L);
    std::vector<Eigen::Triplet<double> > Laugcoeffs;
    for (int k = 0; k < L.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it)
        {
            Laugcoeffs.push_back(Eigen::Triplet<double>(it.row(), it.col(), -smoothreg_ * totarea * it.value()));
        }
    }
    
    for (int i = 0; i < intdofs; i++)
    {
        Laugcoeffs.push_back(Eigen::Triplet<double>(newverts + i, newverts + i, 1e-6));
    }

    Eigen::SparseMatrix<double> Laug(newverts + intdofs, newverts + intdofs);
    Laug.setFromTriplets(Laugcoeffs.begin(), Laugcoeffs.end());
    Mat += Laug;

    Eigen::VectorXd result;
    Eigen::VectorXi toRound(intdofs);
    for (int i = 0; i < intdofs; i++)
    {
        toRound[i] = newverts + i;
    }
    ComisoWrapper(C, Mat, result, rhs, toRound, 1e-6);
    std::cout << "residual: " << (Mat*result - rhs).norm() << std::endl;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver(Mat);
    Eigen::VectorXd unconstr = solver.solve(rhs);
    //result = unconstr;
    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            double intpart;
            double fracpart = std::modf(result[newF(i, j)], &intpart);
            if (fracpart < 0)
                fracpart = 1 + fracpart;
            fracpart -= 0.5;
            theta[surf.data().F(i, j)] = 2.0 * M_PI * fracpart;
        }
    }
    std::cout << "Integer vars: ";
    std::cout << result.segment(newverts, intdofs).transpose() << std::endl;

    //std::cout << result << std::endl;
}