#ifndef TRACE_H
#define TRACE_H

#include <Eigen/Core>
#include <Eigen/Sparse>

class Weave;





class Trace
{
public:
    Trace();
    ~Trace();

    std::vector< Eigen::MatrixXd > curves;
    std::vector< Eigen::MatrixXd > normals;    
    std::vector< Trace_Mode > modes;
    std::vector< Eigen::VectorXd > bending;

    void loadSampledCurves(const std::string &filename);
    void loadGeneratedCurves(std::vector<std::vector<Eigen::Vector3d> > isoLines, 
                                std::vector<std::vector<Eigen::Vector3d> > isoNormal);

    
    
    void logRibbonsToFile(std::string foldername, std::string filename, const Weave &wv);
    void computeIntersections(int curveIdx1, int curveIdx2, 
                              const std::vector<Eigen::MatrixXd> &splitcurves, 
                                    std::vector<Collision> &collisions);
    void save(const std::string &filename);
    void load(const std::string &filename);



};


#endif
