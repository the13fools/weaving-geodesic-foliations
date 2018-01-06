#ifndef WEAVEHOOK_H
#define WEAVEHOOK_H

#include "PhysicsHook.h"
#include "Weave.h"

class WeaveHook : public PhysicsHook
{
public:
    WeaveHook() : PhysicsHook(), weave(NULL), vectorScale(1.0) {}

    virtual void initGUI(igl::viewer::Viewer &viewer)
    {
        viewer.ngui->addGroup("Visualization");
        viewer.ngui->addVariable("Vector Scale", vectorScale);
    }

    virtual void initSimulation()
    {
        if (weave)
            delete weave;
        weave = new Weave("meshes/torus.obj", 3);     
        Handle h;
        h.face = 0;
        h.dir << 1, 0;
        h.field = 0;
        weave->handles.push_back(h);
    }

    virtual void updateRenderGeometry()
    {
        renderQ = weave->V;
        renderF = weave->F;        
        weave->createVisualizationEdges(edgePts, edgeVecs, edgeSegs, edgeColors);
    }

    virtual bool simulateOneStep();    

    virtual void renderRenderGeometry(igl::viewer::Viewer &viewer)
    {
        viewer.data.clear();
        viewer.data.set_mesh(renderQ, renderF);
        Eigen::MatrixXd renderPts = edgePts + vectorScale * edgeVecs;
        viewer.data.set_edges(renderPts, edgeSegs, edgeColors);        
    }

private:
    Weave *weave;

    double vectorScale;

    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;
    Eigen::MatrixXd edgePts;
    Eigen::MatrixXd edgeVecs;
    Eigen::MatrixXi edgeSegs;
    Eigen::MatrixXd edgeColors;    
};

#endif
