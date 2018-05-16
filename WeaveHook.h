#ifndef WEAVEHOOK_H
#define WEAVEHOOK_H

#include "PhysicsHook.h"
#include "Weave.h"
#include "GaussNewton.h"
#include "Trace.h"
#include <string>

#include <igl/unproject_onto_mesh.h>

enum Shading_Enum {
    NONE = 0,
    F1_ENERGY,
    F2_ENERGY,
    F3_ENERGY,
    TOT_ENERGY,
    FUN_VAL,
    CONNECTION_ENERGY
};

class WeaveHook : public PhysicsHook
{
public:
    WeaveHook() : PhysicsHook(), weave(NULL), vectorScale(1.0), normalizeVectors(true)
    {
        // meshName = "meshes/bunny_coarser.obj";
        meshName = "meshes/tet.obj";
        // vectorFieldName = "bunny_coarser_nosing";
        vectorFieldName = "tet.relax";
        traceFile = "example.tr";
        params.lambdacompat = 100;
        params.lambdareg = 1e-3;

        traceIdx = 0;
        traceSign = 1;
        traceSteps = 100;
        traceFaceId = 0;
        isDrawTrace = false;

        hideVectors = false;
        showBending = false;
        showSingularities = false;
    
        handleLocation = 0;
        handleParams = Eigen::VectorXd::Zero(6);
        handleParams(0) = 0;
        handleParams(1) = 1;
        handleParams(2) = 1;
        handleParams(3) = 0;
        handleParams(4) = 1;
        handleParams(5) = -1;

        drawLine = false;
        trace = new Trace();
    }

    virtual void drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu);
    virtual bool mouseClicked(igl::opengl::glfw::Viewer &viewer, int button);

    void reassignPermutations();
    void normalizeFields();
    void serializeVectorField();
    void deserializeVectorField();
    void exportVectorField();
    void augmentField();
    void computeFunc();
    void drawISOLines();
    void resetCutSelection();
    void addCut();
    void removeSingularities();
    void removePrevCut(); 
    void saveTraces();
    void loadTraces();
    void loadSampledTraces();
    
    
    virtual void initSimulation();

    virtual void updateRenderGeometry()
    {
        renderQ = weave->V;
        renderF = weave->F;        
        weave->createVisualizationEdges(edgePts, edgeVecs, edgeSegs, edgeColors);
        weave->createVisualizationCuts(cutPos1, cutPos2);
        faceColors.resize(weave->nFaces(), 3);
        faceColors.setConstant(0.3);
        baseLength = weave->averageEdgeLength;
        curFaceEnergies = tempFaceEnergies;


    weave->handles[0].face = handleLocation; 
    weave->handles[0].dir(0) = handleParams(0); 
    weave->handles[0].dir(1) = handleParams(1); 
    weave->handles[1].face = handleLocation; 
    weave->handles[1].dir(0) = handleParams(2); 
    weave->handles[1].dir(1) = handleParams(3); 
    weave->handles[2].face = handleLocation; 
    weave->handles[2].dir(0) = handleParams(4); 
    weave->handles[2].dir(1) = handleParams(5); 

    }

    virtual bool simulateOneStep();    

    virtual void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer);    

    void setFaceColors(igl::opengl::glfw::Viewer &viewer);
 
    void drawTraceCenterlines(igl::opengl::glfw::Viewer &viewer);
    void drawCuts(igl::opengl::glfw::Viewer &viewer);

    void showCutVertexSelection(igl::opengl::glfw::Viewer &viewer);
    void updateSingularVerts(igl::opengl::glfw::Viewer &viewer);
private:
    std::string meshName;
    Weave *weave;
    SolverParams params;

    std::string traceFile;
    Trace *trace;

    std::vector<std::pair<int, int > > selectedVertices; // (face, vert) pairs
    
    double vectorScale;
    double baseLength;

    Eigen::VectorXd handleParams;
    int handleLocation;

    Eigen::MatrixXd faceColors;
    Eigen::MatrixXd curFaceEnergies;
    Eigen::MatrixXd tempFaceEnergies;
    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
    Eigen::MatrixXd edgePts;
    Eigen::MatrixXd edgeVecs;
    Eigen::MatrixXi edgeSegs;
    Eigen::MatrixXd edgeColors;    
    std::vector<Eigen::Vector3d> renderSelectedVertices; // teal selected vertex spheres
    bool normalizeVectors;
    bool hideVectors;

    Shading_Enum shading_state = Shading_Enum::NONE;
    Trace_Mode trace_state = Trace_Mode::GEODESIC;
    
    int traceIdx;
    int traceSign;
    int traceFaceId;
    int traceSteps;
    // nutton Variables for weave hook
    bool isDrawTrace;
    bool isDeleteLastTrace;
    bool isSaveTrace;
    
    bool showBending;
    bool showSingularities;
    Eigen::MatrixXd singularVerts_topo;
    Eigen::MatrixXd singularVerts_geo;
    Eigen::MatrixXd nonIdentityEdges;
    Eigen::MatrixXd cutPos1; // endpoints of cut edges
    Eigen::MatrixXd cutPos2;

    std::string vectorFieldName;

    bool drawLine;
    std::vector<Eigen::MatrixXd> paths;
    int numISOLines;
    double scalesInit;
};

#endif
