#ifndef WEAVEHOOK_H
#define WEAVEHOOK_H

#include "PhysicsHook.h"
#include "Weave.h"
#include "GaussNewton.h"
#include "LinearSolver.h"
#include "Traces.h"
#include <string>
#include "Surface.h"
#include <igl/unproject_onto_mesh.h>

class CoverMesh;

enum WeaveShading_Enum {
    WS_NONE = 0,
    F1_ENERGY,
    F2_ENERGY,
    F3_ENERGY,
    TOT_ENERGY,
    WS_CONNECTION_ENERGY
};

enum LocalFieldIntegration_Enum {
    LFI_TRIVIAL,     // simple normalization
    LFI_CURLCORRECT, // local curl correction
    LFI_SPECTRAL     // our spectral approach
};

enum GlobalFieldIntegration_Enum {
    GFI_GN,     // our Gauss-Newton code
    GFI_MI      // Bommes et al mixed-integer
};

enum CoverShading_Enum {
    CS_NONE = 0,
    CS_S_VAL,
    FUN_VAL,
    CS_CONNECTION_ENERGY,
    CS_GRAD_DEVIATION
};

enum GUIMode_Enum {
    WEAVE = 0,
    COVER
};

enum Solver_Enum {
    CURLFREE = 0,
    SMOOTH
};

class WeaveHook : public PhysicsHook
{
public:
    WeaveHook() : PhysicsHook(), weave(NULL), cover(NULL), vectorScale(1.0), normalizeVectors(true)
    {
        gui_mode = GUIMode_Enum::WEAVE;
        solver_mode = Solver_Enum::CURLFREE;
        weave_shading_state = WeaveShading_Enum::WS_NONE;
        cover_shading_state = CoverShading_Enum::CS_NONE;
        // meshName = "meshes/bunny_coarser.obj";
//        meshName = "../meshes/sphere_small_regular.obj";
        meshName = "meshes/tet.obj";
        // vectorFieldName = "bunny_coarser_nosing";
        vectorFieldName = "tet.rlx";
        rodFilename = "example.rod";
        exportPrefix = "../final/export/example";
        params.lambdacompat = 100;
        params.lambdareg = 1e-3;
        params.softHandleConstraint = true;
        params.disableCurlConstraint = false;

        params.vizVectorCurl = 1.; // in field surface, vizualization variable
        params.vizCorrectionCurl = 0. ; // in field surface, vizualization variable
        params.vizNormalizeVecs = false;
        params.vizShowCurlSign = false;

     //   ls = new LinearSolver();

        traceIdx = 0;
        traceSign = 1;
        traceSteps = 100;
        traceFaceId = 0;
        
        vectorVisMode = VMM_VFANDDELTA;
        rosyVisMode = RVM_ROSY;
        showSingularities = false;
        wireframe = false;

        targetResolution = 1000;
    
        fieldCount = 1;
        handleLocation = Eigen::VectorXi::Zero(2);
        handleParams = Eigen::VectorXd::Zero(3);
        handleParams(0) = 1;
        handleParams(1) = 1;
        handleParams(2) = 1;

        showCoverCuts = true;
        numISOLines = 1;
        
        local_field_integration_method = LFI_SPECTRAL;
        global_field_integration_method = GFI_GN;
        bommesAniso = 1.0;
        initSReg = 1e-4;
        globalSScale = 0.5;
        globalThetaReg = 1e-4;
        globalAlternations = 10;
        globalPowerIters = 10;

        showTraces = true;
        showRatTraces = true;
        extendTrace = 0.;
        segLen = 0.005;
        maxCurvature = 0.5;
        minRodLen = .1;

        hideCoverVectors = false;
        rosyN = 0;
        desiredRoSyN = 6;
        
        numRandomTraces = 100;
    }

    virtual void drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu);
    virtual bool mouseClicked(igl::opengl::glfw::Viewer &viewer, int button);

    void reassignPermutations();
    void normalizeFields();
    void serializeVectorField();
    void deserializeVectorField();    
    void deserializeVectorFieldOld();
    void deserializePaulField();
    void augmentField();
    void computeFunc();
    void roundCovers();
    void drawISOLines();
    void resetCutSelection();
    void addCut();
    void clearCuts();
    void resample();
    void addHandle();
    void removeHandle();
    void removePrevCut(); 
    void clearTraces();
    void deleteLastTrace();
    void computeTrace();   
    void computeRandomTraces(int numtraces);   
    void rationalizeTraces();
    void saveRods();
    void exportForRendering();
    void convertToRoSy();
    void splitFromRoSy();
    
    virtual void initSimulation();

    virtual void updateRenderGeometry();

    virtual bool simulateOneStep();    

    virtual void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer);    

    void setFaceColorsWeave(igl::opengl::glfw::Viewer &viewer);
    void setFaceColorsCover(igl::opengl::glfw::Viewer &viewer);
 
    void drawCuts(igl::opengl::glfw::Viewer &viewer);

    void showCutVertexSelection(igl::opengl::glfw::Viewer &viewer);
    void updateSingularVerts(igl::opengl::glfw::Viewer &viewer);

private:
    void clear();
    std::string meshName;
    Weave *weave;
    CoverMesh *cover;
    SolverParams params;

    TraceSet traces;

    std::vector<std::pair<int, int > > selectedVertices; // (face, vert) pairs
    
    double vectorScale;
    
    Eigen::VectorXd handleParams;
    Eigen::VectorXi handleLocation;

    LinearSolver ls;

    Eigen::MatrixXd curFaceEnergies;
    Eigen::MatrixXd tempFaceEnergies;
    Eigen::MatrixXd renderQWeave;
    Eigen::MatrixXi renderFWeave;
    Eigen::MatrixXd edgePtsWeave;
    Eigen::MatrixXi edgeSegsWeave;
    Eigen::MatrixXd edgeColorsWeave;    
    Eigen::MatrixXd edgePtsCover;
    Eigen::MatrixXi edgeSegsCover;
    Eigen::MatrixXd edgeColorsCover;    
    std::vector<Eigen::Vector3d> renderSelectedVertices; // teal selected vertex spheres
    VectorVisualizationMode vectorVisMode;
    RoSyVisualizationMode rosyVisMode;
    bool normalizeVectors;
    bool showCoverCuts;
    bool wireframe;

    Eigen::MatrixXd renderQCover;
    Eigen::MatrixXi renderFCover;
    
    Solver_Enum solver_mode;
    GUIMode_Enum gui_mode;
    WeaveShading_Enum weave_shading_state;
    CoverShading_Enum cover_shading_state;
    Trace_Mode trace_state = Trace_Mode::GEODESIC;
    
    int traceIdx;
    int traceSign;
    int traceFaceId;
    int traceSteps;
    int targetResolution;
    
    bool showSingularities;
    Eigen::MatrixXd singularVerts_topo;
    Eigen::MatrixXd singularVerts_geo;
    Eigen::MatrixXd nonIdentity1Weave;
    Eigen::MatrixXd nonIdentity2Weave;
    Eigen::MatrixXd cutPos1Weave; // endpoints of cut edges
    Eigen::MatrixXd cutPos2Weave;
    Eigen::MatrixXd cutPos1Cover;
    Eigen::MatrixXd cutPos2Cover;
    Eigen::MatrixXd cutColorsCover;

    std::string vectorFieldName;
    std::string exportPrefix;

    bool showTraces;
    bool showRatTraces;
    double extendTrace;
    double segLen;
    double maxCurvature;
    double minRodLen;
    // isolines on the split mesh
    Eigen::MatrixXd pathstarts;
    Eigen::MatrixXd pathends;
    // traces on the single mesh
    Eigen::MatrixXd tracestarts;
    Eigen::MatrixXd traceends;
    Eigen::MatrixXd tracecolors;
    
    std::string rodFilename;

    Eigen::MatrixXd rattracestarts;
    Eigen::MatrixXd rattraceends;
    Eigen::MatrixXd ratcollisions;

    LocalFieldIntegration_Enum local_field_integration_method;
    GlobalFieldIntegration_Enum global_field_integration_method;

    int numISOLines;
    double bommesAniso;
    double initSReg;
    double globalSScale;
    double globalThetaReg;
    int globalAlternations;
    int globalPowerIters;

    int fieldCount;
    bool hideCoverVectors;
    
    int numRandomTraces;

    int rosyN;
    int desiredRoSyN;
};

#endif
