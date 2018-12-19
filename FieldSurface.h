#ifndef FIELDSURFACE_H
#define FIELDSURFACE_H

#include <Eigen/Core>
#include <vector>
#include <iostream>
#include <map>
#include "Surface.h"

#include "GaussNewton.h"

/*
 * Surface on one of more families of vector fields live (with permutations mapping between the along edges).
 */

class FieldSurface : public Surface
{
public:
    FieldSurface(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, int numFields);

    int nFields() const { return nFields_; }

    // index of entry in vectorFields of various quantities on a face
    int vidx(int face, int field) const;                
    int betaidx(int face, int field) const;
    int alphaidx(int face, int field) const;
    
    // accessors into vectorFields
    Eigen::Vector2d v(int face, int field) const;
    Eigen::Vector2d beta(int face, int field) const;
    double alpha(int face, int field) const;
    double getGeodesicEnergy(SolverParams params);

    const Eigen::MatrixXi Ps(int edge) const;

    void normalizeFields(); // make all vectors unit-length

    // returns copy of the surface with all deleted faces completely removed, and all appropriate data structures reindexed
    // faceMap maps old face ids to new face ids (deleted faces are not in the map).
    // likewise for vertMap
    FieldSurface *removeDeletedFacesFromMesh(std::map<int, int> &faceMap, std::map<int, int> &vertMap) const; 
    
    void serialize(std::ostream &ofs) const;
    static FieldSurface *deserialize(std::istream &is);    
    
    // mark a vertex and all neighboring faces as deleted
    void deleteVertex(int vid);
    
    void setFaceDeleted(int fid, bool newstatus);
    void undeleteAllFaces();    
    bool isFaceDeleted(int fid) const {return faceDeleted_[fid];}
    int numUndeletedFaces() const;

    // compute an energy on each face that measures the failure of the vector fields on that face to parallel transport to the
    // equivalent vector on the neighboring faces, using the trivial connection between faces
    // each value will be in the range [0, 3*m*PI].
    void connectionEnergy(Eigen::VectorXd &energies, double threshold, SolverParams params);


    Eigen::VectorXd vectorFields;    // the vector fields, unrolled into a vector:
                                     // first 2m|F| entries: first vector on face 1, second vector on face 1, third vector on face 1, ..., last vector on face m
                                     // next 2m|F| entries: first beta vector on face 1, ...
                                     // next m|F| entries: first alpha on face 1, ...
    std::vector<Eigen::MatrixXi> Ps_; // for each edge i, maps indices from triangle E(i,1) to indices in triangle E(i,0), with sign. I.e. the vector on E(i,1) corresponding to vector j on E(i,0) is \sum_k Ps[i](j,k) v(E(i,1),k)
    int nFields_;

private:
    double geodesicEnergy_;
    std::vector<bool> faceDeleted_;
};

#endif
