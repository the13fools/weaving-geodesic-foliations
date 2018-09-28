#include "FieldSurface.h"
#include "Surface.h"
#include <set>
#include <map>
#include <igl/remove_unreferenced.h>
#include <Eigen/Dense>

FieldSurface::FieldSurface(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, int numFields) : Surface(V,F), nFields_(numFields)
{
    int nfaces = data().F.rows();
    // initialize vector fields
    vectorFields.resize(6*nfaces*numFields);
    vectorFields.setZero();
    vectorFields.segment(0, 2 * nfaces*numFields).setRandom();

    Eigen::Vector3d target;
    target << 1, 1, 1;
    for(int i = 0; i < nfaces; i++)
    {
        Eigen::Matrix<double, 3, 2> B = this->data().Bs[i];
        // To see this, start with Bv = target vector (tv).
        //    B^TBv = B^T tv >==> v = (B^TB)^-1B tv
        Eigen::MatrixXd inverse_trans = (B.transpose() * B).inverse() * B.transpose();
        std::cout << inverse_trans << std::endl;
        for(int m = 0; m < numFields; m++)
        {
            vectorFields.segment<2>(2*numFields*i + m) = inverse_trans * target;
        } 


    }

    vectorFields.segment(5*nfaces*numFields, nfaces*numFields) = Eigen::VectorXd::Constant(nfaces * numFields, 1.);
    normalizeFields();

    // initialize permutation matrices
    int nedges = nEdges();
    Ps_.resize(nedges);
    for (int i = 0; i < nedges; i++)
    {
        Ps_[i].resize(numFields, numFields);
        Ps_[i].setIdentity();
    }
    
    faceDeleted_.resize(nfaces);
    for(int i=0; i<nfaces; i++)
        faceDeleted_[i] = false;
}


int FieldSurface::vidx(int face, int field) const
{
    return (2 * nFields()*face + 2 * field);
}

Eigen::Vector2d FieldSurface::v(int face, int field) const
{
    return vectorFields.segment<2>(vidx(face,field));
}

int FieldSurface::betaidx(int face, int field) const
{
    return 2 * nFields()*nFaces() + 2 * nFields()*face + 2 * field;
}
Eigen::Vector2d FieldSurface::beta(int face, int field) const
{
    return vectorFields.segment<2>(betaidx(face,field));
}

int FieldSurface::alphaidx(int face, int field) const
{
    return 4 * nFields()*nFaces() + nFields()*face + field;
}

int FieldSurface::sidx(int face, int field) const
{
    return 5 * nFields()*nFaces() + nFields()*face + field;
}

double FieldSurface::alpha(int face, int field) const
{
    return vectorFields[alphaidx(face,field)];
}

double FieldSurface::sval(int face, int field) const
{
    return vectorFields[sidx(face,field)];
}


void FieldSurface::normalizeFields()
{
    int nfaces = nFaces();
    int m = nFields();
    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < m; j++)
        {
            Eigen::Vector2d vif = v(i, j);
            double norm = sqrt(vif.transpose() * data().Bs[i].transpose() * data().Bs[i] * vif);
            vectorFields.segment<2>(vidx(i, j)) /= norm;
        }
    }
}

FieldSurface *FieldSurface::removeDeletedFacesFromMesh(std::map<int, int> &faceMap, std::map<int, int> &vertMap) const
{
    faceMap.clear();
    vertMap.clear();
    
    std::set<int> facesToDelete;
    
    std::map<std::pair<int, int>, int> edgeMap;
    for (int e = 0; e < nEdges(); e++) 
    { 
        std::pair<int, int> p(data().edgeVerts(e,0), data().edgeVerts(e,1));
        edgeMap[p] = e; 
    }

    for(int i=0; i<nFaces(); i++)
    {
        if(isFaceDeleted(i))
            facesToDelete.insert(i);
    }

    std::vector<int> faceIds;
    for (std::set<int>::iterator it = facesToDelete.begin(); it != facesToDelete.end(); ++it)
        faceIds.push_back(*it);

    if (faceIds.empty())
    {
        // nothing to do
        FieldSurface *ret = new FieldSurface(data().V, data().F, nFields());
        ret->vectorFields = vectorFields;
        ret->Ps_ = Ps_;
        for(int i=0; i<nFaces(); i++)
            faceMap[i] = i;
        for(int i=0; i<nVerts(); i++)
            vertMap[i] = i;
        return ret;
    }
    
    int fieldIdx = 0;
    int newNFaces = nFaces() - faceIds.size();
    Eigen::VectorXd vectorFields_clean = Eigen::VectorXd::Zero( 6*nFields()*newNFaces );
    Eigen::MatrixXi F_temp = Eigen::MatrixXi::Zero(newNFaces, 3); 

    std::vector<bool> newdeleted(newNFaces);
    
    for (int i = 0; i < nFaces(); i++)   
    { 
        if (facesToDelete.count(i))
            continue;
         // TODO ADD THE S VARIABLE HERE
        // vec field
        vectorFields_clean.segment(2*fieldIdx*nFields(), 2*nFields()) = vectorFields.segment(2*i*nFields(), 2*nFields());
        // beta
        vectorFields_clean.segment(2*fieldIdx*nFields() + 2*newNFaces*nFields(), 2*nFields()) 
            = vectorFields.segment(2*i*nFields() + 2*nFaces()*nFields(), 2*nFields() );
        // alpha
        vectorFields_clean.segment(fieldIdx*nFields() + 4*newNFaces*nFields(), nFields()) 
            = vectorFields.segment(i * nFields() + 4*nFaces()*nFields(), nFields() );
        // faces 
        F_temp.row(fieldIdx) = data().F.row(i);
        newdeleted[fieldIdx] = faceDeleted_[i];
        faceMap[i] = fieldIdx;
        fieldIdx++;
    }

    Eigen::MatrixXd V_new;
    Eigen::MatrixXi F_new;
    Eigen::VectorXi marked; 
    Eigen::VectorXi vertMapVec; 

    igl::remove_unreferenced(data().V, F_temp, V_new, F_new, marked, vertMapVec);

    FieldSurface *result = new FieldSurface(V_new, F_new, nFields());
    result->vectorFields = vectorFields_clean;
    result->faceDeleted_ = newdeleted;

    for(int i=0; i<vertMapVec.rows(); i++)
    {
        vertMap[vertMapVec(i)] = i;
    }
           
    std::vector<Eigen::MatrixXi> Ps_new;
    for( int i = 0; i < result->nEdges(); i++) 
    {
        int v0 = vertMapVec( result->data().edgeVerts(i, 0) );
        int v1 = vertMapVec( result->data().edgeVerts(i, 1) );
        if ( v0 > v1 ) 
        { 
            std::swap(v0, v1);
        }
        std::pair<int, int> p(v0, v1);
        int oldEdge = edgeMap[p];
        Ps_new.push_back(Ps_[oldEdge]);
    }
    result->Ps_ = Ps_new;

    return result;
}

const Eigen::MatrixXi FieldSurface::Ps(int edge) const
{
    return Ps_[edge];
}

void FieldSurface::serialize(std::ostream &os) const
{
    int magic = -777;
    os.write((char *)&magic, sizeof(int));
    int version = 1;
    os.write((char *)&version, sizeof(int));
    int nverts = nVerts();
    os.write((char *)&nverts, sizeof(int));
    int nfaces = nFaces();
    os.write((char *)&nfaces, sizeof(int));
    int nfields = nFields();
    os.write((char *)&nfields, sizeof(int));
    int nperms = Ps_.size();
    os.write((char *)&nperms, sizeof(int));
    for (int i = 0; i < nverts; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            double tmp = data().V(i, j);
            os.write((char *)&tmp, sizeof(double));
        }
    }
    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int tmp = data().F(i, j);
            os.write((char *)&tmp, sizeof(int));
        }
    }
    for (int i = 0; i < 6 * nfaces*nfields; i++)
    {
        double tmp = vectorFields[i];
        os.write((char *)&tmp, sizeof(double));
    }
    for (int i = 0; i < nperms; i++)
    {
        for (int j = 0; j < nfields; j++)
        {
            for (int k = 0; k < nfields; k++)
            {
                int tmp = Ps_[i](j, k);
                os.write((char *)&tmp, sizeof(int));
            }
        }
    }
    for(int i=0; i<nfaces; i++)
    {
        int deleted = faceDeleted_[i] ? 1 : 0;
        os.write((char *)&deleted, sizeof(int));
    }
}

FieldSurface *FieldSurface::deserialize(std::istream &is)
{
    int version = 0;
    int magic;
    int nverts;
    is.read((char *)&magic, sizeof(int));
    if(magic == -777)
    {
        is.read((char *)&version, sizeof(int));
        is.read((char *)&nverts, sizeof(int));
    }
    else
    {
        // old format
        nverts = magic;
    }
    Eigen::MatrixXd V(nverts, 3);
    int nfaces;
    is.read((char *)&nfaces, sizeof(int));
    Eigen::MatrixXi F(nfaces, 3);
    int nfields;
    is.read((char *)&nfields, sizeof(int));
    int nperms;
    is.read((char *)&nperms, sizeof(int));
    if (!is)
        return NULL;
    for (int i = 0; i < nverts; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            double tmp;
            is.read((char *)&tmp, sizeof(double));
            V(i, j) = tmp;
        }
    }
    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int tmp;
            is.read((char *)&tmp, sizeof(int));
            F(i, j) = tmp;
        }
    }
    if (!is)
        return NULL;
    FieldSurface *ret = new FieldSurface(V, F, nfields);
    for (int i = 0; i < 5 * nfaces*nfields; i++)
    {
        double tmp;
        is.read((char *)&tmp, sizeof(double));
        ret->vectorFields[i] = tmp;
    }
    for (int i = 0; i < nperms; i++)
    {
        for (int j = 0; j < nfields; j++)
        {
            for (int k = 0; k < nfields; k++)
            {
                int tmp;
                is.read((char *)&tmp, sizeof(int));
                ret->Ps_[i](j, k) = tmp;
            }
        }
    }
    if(version > 0)
    {
        for(int i=0; i<nfaces; i++)
        {
            int deleted;
            is.read((char *)&deleted, sizeof(int));
            ret->faceDeleted_[i] = deleted;
        }
    }
    if (!is)
    {
        delete ret;
        return NULL;
    }
    return ret;
}


void FieldSurface::connectionEnergy(Eigen::VectorXd &energies)
{
    energies.resize(nFaces());
    energies.setZero();

    int nedges = nEdges();
    int nfields = nFields();
    for(int i=0; i<nedges; i++)
    {
        if(data().E(i,0) == -1 || data().E(i,1) == -1)
            continue;

        int face = data().E(i,0);
        double facearea = faceArea(face);
        int opp = data().E(i,1);
        double opparea = faceArea(opp);

        for(int j=0; j<nfields; j++)
        {
            Eigen::Vector2d vec = v(face, j);
            Eigen::Vector2d oppvec(0,0);
            for(int k=0; k<nfields; k++)
                oppvec += Ps_[i](j,k)*v(opp,k);
            Eigen::Vector2d mappedvec = data().Ts.block<2,2>(2*i,0) * vec;
            // mappedvec and oppvec now both live on face opp.
            // compute the angle between them

            Eigen::Vector3d v1 = data().Bs[opp]*mappedvec;
            Eigen::Vector3d v2 = data().Bs[opp]*oppvec;
            Eigen::Vector3d n = faceNormal(opp);
            double angle = 2.0 * atan2(v1.cross(v2).dot(n), v1.norm() * v2.norm() + v1.dot(v2));
            energies[face] += facearea*fabs(angle);
            energies[opp] += opparea*fabs(angle);
        }
    }
}

void FieldSurface::deleteVertex(int vid)
{
    int nfaces = data().F.rows();
    for(int i=0; i<nfaces; i++)
    {
        for(int j=0; j<3; j++)
        {
            if(data().F(i,j) == vid)
                faceDeleted_[i] = true;
        }
    }
}


void FieldSurface::undeleteAllFaces()
{
    int nfaces = data().F.rows();
    for(int i=0; i<nfaces; i++)
    {
        faceDeleted_[i] = false;
    }
}

void FieldSurface::setFaceDeleted(int fid, bool newstatus)
{
    faceDeleted_[fid] = newstatus;
}

int FieldSurface::numUndeletedFaces() const
{
    int ret = 0;
    int nfaces = data().F.rows();
    for(int i=0; i<nfaces; i++)
        if(!faceDeleted_[i])
            ret++;
    return ret;
}
