#include "FaceBased.h"
#include <Eigen/Geometry>

Eigen::Vector3d faceNormal(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, int faceidx)
{
    Eigen::Vector3d p0 = V.row(F(faceidx, 0));
    Eigen::Vector3d p1 = V.row(F(faceidx, 1));
    Eigen::Vector3d p2 = V.row(F(faceidx, 2));
    Eigen::Vector3d n = (p1 - p0).cross(p2 - p0);
    n /= n.norm();
    return n;
}

static void oneFaceGradientMatrix(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const Eigen::MatrixXi &F_edges, int faceidx, std::vector<Eigen::Triplet<double> > &Mcoeffs)
{
    Mcoeffs.clear();
    //collect neighboring face indices
    int nbs[3];
    for (int i = 0; i < 3; i++)
    {
        const Eigen::Vector4i &einfo = E.row(F_edges(faceidx,i));
        nbs[i] = (einfo[2] == faceidx ? einfo[3] : einfo[2]);
        // if face faceidx is on the boundary, return empty matrix
        if (nbs[i] == -1)
            return;
    }
    Eigen::Vector3d n = faceNormal(F, V, faceidx);
    // collect for each neighboring face: (1) the midpoint of the common edge, (2) unit vector in direction of common edge,
    // (3) the face normal, (4) the centroid of neighboring face
    Eigen::Vector3d midpts[3];
    Eigen::Vector3d commones[3];
    Eigen::Vector3d nbns[3];
    Eigen::Vector3d centroids[3];
    for (int i = 0; i < 3; i++)
    {
        int v1 = (i + 1) % 3;
        int v2 = (i + 2) % 3;
        midpts[i] = 0.5 * (V.row(F(faceidx, v1)) + V.row(F(faceidx, v2)));
        commones[i] = V.row(F(faceidx, v2)) - V.row(F(faceidx, v1));
        commones[i] /= commones[i].norm();
        nbns[i] = faceNormal(F, V, nbs[i]);
        centroids[i].setZero();
        for (int j = 0; j < 3; j++)
        {
            centroids[i] += V.row(F(nbs[i], j));            
        }
        centroids[i] /= 3.0;
    }

    // main face centroid
    Eigen::Vector3d c;
    c.setZero();
    for (int i = 0; i < 3; i++)
    {
        c += V.row(F(faceidx, i));
    }
    c /= 3.0;

    double centroidDist = 0.;
    double edgeLen = 0.;
    for (int i = 0; i < 3; i++)
    { 
        int v1 = (i + 1) % 3;
        int v2 = (i + 2) % 3;
    	centroidDist += (centroids[i] - c).norm();
	edgeLen += (V.row(F(faceidx, v2)) - V.row(F(faceidx, v1))).norm();
    }

    //rotate each centroid into the plane of the triangle faceidx and compute ci minus c
    Eigen::Vector3d ciminusc[3];
    for (int i = 0; i < 3; i++)
    {
        Eigen::Vector3d t1 = n.cross(commones[i]);
        Eigen::Vector3d t2 = nbns[i].cross(commones[i]);
        Eigen::Vector3d diff = centroids[i] - midpts[i];
        double alpha = commones[i].dot(diff);
        double beta = t2.dot(diff);
        ciminusc[i] = midpts[i] + alpha * commones[i] + beta * t1 - c;
     
        int v1 = (i + 1) % 3;
        int v2 = (i + 2) % 3;	
	double edge = (V.row(F(faceidx, v2)) - V.row(F(faceidx, v1))).norm();
	ciminusc[i] *= (edge / edgeLen) * (centroids[i] - c).norm() / centroidDist;
    }

    // compute N
    Eigen::Vector3d N[3];
    for (int i = 0; i < 3; i++)
    {
        int v1 = (i + 2) % 3;
        int v2 = (i + 1) % 3;
        N[i] = 1.0 / 3.0 * (n.cross(ciminusc[v1]) / (ciminusc[v1].cross(ciminusc[i]).dot(n)) - n.cross(ciminusc[v2]) / (ciminusc[i].cross(ciminusc[v2]).dot(n)));
    }

    // finally, set up M
    for (int i = 0; i < 3; i++)
    {
        Mcoeffs.push_back(Eigen::Triplet<double>(i, faceidx, -N[0][i] - N[1][i] - N[2][i]));
    }
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Mcoeffs.push_back(Eigen::Triplet<double>(j, nbs[i], N[i][j]));
        }
    }
}

void computeGradientMatrices(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V, const Eigen::MatrixXi &E, const Eigen::MatrixXi &F_edges, std::vector<Eigen::SparseMatrix<double> > &Ms)
{
    int nfaces = (int)F.rows();
    Ms.resize(nfaces);
    for (int i = 0; i < nfaces; i++)
    {
        std::vector<Eigen::Triplet<double> > Mcoeffs;
        oneFaceGradientMatrix(F, V, E, F_edges, i, Mcoeffs);
        Ms[i].resize(3, nfaces);
        Ms[i].setFromTriplets(Mcoeffs.begin(), Mcoeffs.end());
    }
}
