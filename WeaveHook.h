#ifndef WEAVEHOOK_H
#define WEAVEHOOK_H

#include "PhysicsHook.h"
#include "Weave.h"
#include "GaussNewton.h"
#include "Trace.h"
#include <string>

enum Shading_Enum {
    NONE = 0,
    F1_ENERGY,
    F2_ENERGY,
    F3_ENERGY,
    TOT_ENERGY
};

class WeaveHook : public PhysicsHook
{
public:
    WeaveHook() : PhysicsHook(), weave(NULL), vectorScale(1.0), normalizeVectors(true)
    {
        meshName = "meshes/torus.obj";
        params.lambdacompat = 100;
        params.lambdareg = 1e-3;

        traceIdx = 0;
        traceSign = 1;
        traceSteps = 100;
        traceFaceId = 0;
        isDrawTrace = false;

        hideVectors = false;
        showBending = false;

        trace = new Trace();
    }

    virtual void initGUI(igl::viewer::Viewer &viewer)
    {
        viewer.ngui->addVariable("Mesh", meshName);
        viewer.ngui->addGroup("Visualization");
        viewer.ngui->addVariable("Vector Scale", vectorScale);
        viewer.ngui->addVariable("Normalize Vectors", normalizeVectors);
        viewer.ngui->addVariable("Hide Vectors", hideVectors);
        viewer.ngui->addVariable("Shading", shading_state, true)
                   ->setItems({"None", "F1 Energy", "F2 Energy", "F3 Energy", "Total Energy"});
        
        viewer.ngui->addGroup("Solver Parameters");
        viewer.ngui->addVariable("Compatilibity Lambda", params.lambdacompat);
        viewer.ngui->addVariable("Tikhonov Reg", params.lambdareg);
        viewer.ngui->addButton("Reassign Permutations", std::bind(&WeaveHook::reassignPermutations, this));
        
//	viewer.ngui->addVariable("Trace Field/Geo", isTraceField);
	
        viewer.ngui->addGroup("Tracing Controls");
        viewer.ngui->addVariable("Trace Face", traceFaceId);	
        viewer.ngui->addVariable("Trace Steps", traceSteps);	
        viewer.ngui->addVariable("Trace Field", traceIdx);
        viewer.ngui->addVariable("Trace Sign", traceSign);
        viewer.ngui->addVariable("Trace Mode", trace_state, true)
                   ->setItems({"Geodesic", "Field"}); 
        viewer.ngui->addVariable("Show Bending", showBending);
    }

    void reassignPermutations();
    void normalizeFields();

    virtual void initSimulation()
    {
        if (weave)
            delete weave;
        weave = new Weave(meshName, 3);     
        Handle h;
        h.face = 0;
        h.dir << 1, 0;
        h.field = 2;
        weave->addHandle(h);
        h.face = 0;
        h.dir << 0, 1;
        h.field = 1;
        weave->addHandle(h);
        h.face = 0;
        h.dir << 1, -1;
        h.field = 0;
        weave->addHandle(h);
	curFaceEnergies = Eigen::MatrixXd::Zero(3,3);
    }

    virtual void updateRenderGeometry()
    {
        renderQ = weave->V;
        renderF = weave->F;        
        weave->createVisualizationEdges(edgePts, edgeVecs, edgeSegs, edgeColors);
        faceColors.resize(weave->nFaces(), 3);
        faceColors.setConstant(0.3);
        baseLength = weave->averageEdgeLength;
        curFaceEnergies = tempFaceEnergies;
    }

    virtual bool simulateOneStep();    

    virtual void renderRenderGeometry(igl::viewer::Viewer &viewer);    

    void setFaceColors(igl::viewer::Viewer &viewer);
 
    void drawTraceCenterlines(igl::viewer::Viewer &viewer);
private:
    std::string meshName;
    Weave *weave;
    SolverParams params;
    Trace *trace;

    double vectorScale;
    double baseLength;

    Eigen::MatrixXd faceColors;
    Eigen::MatrixXd curFaceEnergies;
    Eigen::MatrixXd tempFaceEnergies;
    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
    Eigen::MatrixXd edgePts;
    Eigen::MatrixXd edgeVecs;
    Eigen::MatrixXi edgeSegs;
    Eigen::MatrixXd edgeColors;    
    bool normalizeVectors;
    bool hideVectors;

    Shading_Enum shading_state = Shading_Enum::NONE;
    Trace_Mode trace_state = Trace_Mode::GEODESIC;

    bool isTraceField; // this controls if we trace a geodesic or the field

    int traceIdx;
    int traceSign;
    int traceFaceId;
    int traceSteps;
    
    bool showBending;
};

#endif
