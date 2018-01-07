#ifndef WEAVEHOOK_H
#define WEAVEHOOK_H

#include "PhysicsHook.h"
#include "Weave.h"
#include "GaussNewton.h"
#include <string>

class WeaveHook : public PhysicsHook
{
public:
    WeaveHook() : PhysicsHook(), weave(NULL), vectorScale(1.0), normalizeVectors(true)
    {
        meshName = "meshes/torus.obj";
        params.lambdacompat = 100;
        params.lambdareg = 1e-3;
    }

    virtual void initGUI(igl::viewer::Viewer &viewer)
    {
        viewer.ngui->addVariable("Mesh", meshName);
        viewer.ngui->addGroup("Visualization");
        viewer.ngui->addVariable("Vector Scale", vectorScale);
        viewer.ngui->addVariable("Normalize Vectors", normalizeVectors);
        viewer.ngui->addGroup("Solver Parameters");
        viewer.ngui->addVariable("Compatilibity Lambda", params.lambdacompat);
        viewer.ngui->addVariable("Tikhonov Reg", params.lambdareg);
    }

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
    }

    virtual void updateRenderGeometry()
    {
        renderQ = weave->V;
        renderF = weave->F;        
        weave->createVisualizationEdges(edgePts, edgeVecs, edgeSegs, edgeColors);
        faceColors.resize(weave->nFaces(), 3);
        faceColors.setConstant(0.5);
        baseLength = weave->averageEdgeLength;
    }

    virtual bool simulateOneStep();    

    virtual void renderRenderGeometry(igl::viewer::Viewer &viewer);    

private:
    std::string meshName;
    Weave *weave;
    SolverParams params;

    double vectorScale;
    double baseLength;

    Eigen::MatrixXd faceColors;
    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
    Eigen::MatrixXd edgePts;
    Eigen::MatrixXd edgeVecs;
    Eigen::MatrixXi edgeSegs;
    Eigen::MatrixXd edgeColors;    
    bool normalizeVectors;
};

#endif
