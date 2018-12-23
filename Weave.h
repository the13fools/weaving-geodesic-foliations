#ifndef WEAVE_H
#define WEAVE_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include "FieldSurface.h"

enum VectorVisualizationMode
{
    VMM_NOTHING,
    VMM_VF,
    VMM_VFANDDELTA,
    VMM_VFPLUSDELTA
};

enum RoSyVisualizationMode
{
    RVM_NOTHING,
    RVM_ROSY,
    RVM_REPVEC
};


// handles on the mesh
struct Handle
{
    int face; // which face the handle is on
    int field; // which of the m fields we're prescribing
    Eigen::Vector2d dir; // the vector itself in triangle barycentric coordinates
};

struct Cut
{ 
    std::vector<std::pair<int, int> > path; // (edge, orientation) list        
};

class CoverMesh;

class Weave
{
public:
    Weave(const std::string &objname, int m);
    Weave(const Eigen::MatrixXd Vtmp, const Eigen::MatrixXi Ftmp, int m);
    Weave(const Weave &w);
    ~Weave();

    Weave &operator=(const Weave &w);
   
    FieldSurface *fs;
    
    std::vector<Handle> handles; // handles on the vector fields
 
    std::vector<Cut> cuts; // list of cuts 

    bool addHandle(Handle h);  // this method will add a handle, also taking care to normalize the handle vector length

    int nHandles() const { return handles.size(); }    

    
    bool fixFields;  // Do not allow vectors to change in optimization.  

    void createVisualizationEdges(
        Eigen::MatrixXd &edgePts, 
        Eigen::MatrixXi &edgeSegs, 
        Eigen::MatrixXd &colors,
        VectorVisualizationMode mode,
        bool normalizeVectors,
        double baseVectorLength // ignored if normalizeVectors=true
    );

    void createVisualizationEdges(
        Eigen::MatrixXd &edgePts, 
        Eigen::MatrixXi &edgeSegs, 
        Eigen::MatrixXd &colors,
        RoSyVisualizationMode mode,
        bool normalizeVectors,
        double baseVectorLength // ignored if normalizeVectors=true
    );

    void createVisualizationCuts(Eigen::MatrixXd &cutPts1, Eigen::MatrixXd &cutPts2);
    
    CoverMesh *createCover() const;

    void serialize(std::ostream &os);    
    void deserialize(std::istream &is);
    void deserializeOldRelaxFile(std::istream &is);

    void convertToRoSy();

private:
    // scale mesh to unit size
    void centerAndScale(Eigen::MatrixXd &V);  

    std::vector<long> _BFS_adj_list(std::vector<std::vector<long> > & relaxadj_list, int i) const;
    std::vector<Eigen::MatrixXd> _augmentPs() const;
};

#endif
