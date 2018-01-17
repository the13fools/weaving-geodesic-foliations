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
    TOT_ENERGY
};

class WeaveHook : public PhysicsHook
{
public:
    WeaveHook() : PhysicsHook(), weave(NULL), vectorScale(1.0), normalizeVectors(true)
    {
        meshName = "meshes/bunny.obj";
        vectorFieldName = "artery.txt";
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

        viewer.ngui->addGroup("Save/Load Field");
        viewer.ngui->addVariable("Filename", vectorFieldName);
        viewer.ngui->addButton("Save Field", std::bind(&WeaveHook::serializeVectorField, this));
        viewer.ngui->addButton("Load Field", std::bind(&WeaveHook::deserializeVectorField, this)); 
        viewer.ngui->addButton("Export Field", std::bind(&WeaveHook::exportVectorField, this));

	viewer.ngui->addGroup("Add Cut");

        viewer.ngui->addWindow(Eigen::Vector2i(300, 10), "Manipulate");	
        viewer.ngui->addGroup("Tracing Controls");
        viewer.ngui->addVariable("Trace Face", traceFaceId);	
        viewer.ngui->addVariable("Trace Steps", traceSteps);	
        viewer.ngui->addVariable("Trace Field", traceIdx);
        viewer.ngui->addVariable("Trace Sign", traceSign);
        viewer.ngui->addVariable("Trace Mode", trace_state, true)
                   ->setItems({"Geodesic", "Field"}); 
        viewer.ngui->addVariable("Show Bending", showBending);
   //     viewer.ngui->addVariable("Show Singularities", showSingularities);

        viewer.callback_mouse_down =
	    [this](igl::viewer::Viewer& viewer, int, int)->bool
	  {
	    int fid;
	    Eigen::Vector3f bc;
	    // Cast a ray in the view direction starting from the mouse position
	    double x = viewer.current_mouse_x;
	    double y = viewer.core.viewport(3) - viewer.current_mouse_y;
	    if(igl::unproject_onto_mesh(Eigen::Vector2f(x,y), viewer.core.view * viewer.core.model,
	      viewer.core.proj, viewer.core.viewport, this->weave->V, this->weave->F, fid, bc))
	    {

	      clicked.row(fid)<<1,0,0;
	      vertexSelect(fid) = (vertexSelect(fid) + 1) % 3;
              std::cout << fid << " " << vertexSelect(fid) << "\n"; 
		// paint hit red
//	      clicked.row(fid)<<1,0,0;
//	      edgeSelect(fid) = (edgeSelect(fid) + 1) % 3;
//     	      params.edgeWeights(this->weave->faceEdges(fid, edgeSelect(fid))) = 0;    
//     	      params.edgeWeights(this->weave->faceEdges(fid, ( edgeSelect(fid) + 2 ) % 3 )) = 1;    
//	      std::cout << this->weave->faceEdges(fid, edgeSelect(fid)) << " " << this->weave->Ps.size() << "\n";
//              this->weave->Ps[this->weave->faceEdges(fid, edgeSelect(fid) )] = this->m;
//              this->weave->Ps[this->weave->faceEdges(fid, ( edgeSelect(fid) + 2 ) % 3 )] = this->idmat;
	      return true;
	    }
	    return false;
	  };
	  std::cout<<R"(Usage:
	  [click]  Pick face on shape

	)";

    }

    void reassignPermutations();
    void normalizeFields();
    void serializeVectorField();
    void deserializeVectorField();
    void exportVectorField();
    
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
        clicked = Eigen::MatrixXd::Constant(weave->nFaces(), 3, .7);
	params.edgeWeights = Eigen::VectorXd::Constant(weave->nEdges(), 1);
	edgeSelect = Eigen::VectorXi::Constant(weave->nFaces(), -1);
	vertexSelect = Eigen::VectorXi::Constant(weave->nFaces(), -1);
    }

    virtual bool simulateOneStep();    

    virtual void renderRenderGeometry(igl::viewer::Viewer &viewer);    

    void setFaceColors(igl::viewer::Viewer &viewer);
 
    void drawTraceCenterlines(igl::viewer::Viewer &viewer);
    void showCutVertexSelection(igl::viewer::Viewer &viewer);
    void updateSingularVerts(igl::viewer::Viewer &viewer);
private:
    std::string meshName;
    Weave *weave;
    SolverParams params;
    Trace *trace;

    Eigen::VectorXi edgeSelect;
    Eigen::VectorXi vertexSelect;

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
    Eigen::MatrixXd clicked;    
    bool normalizeVectors;
    bool hideVectors;

    Shading_Enum shading_state = Shading_Enum::NONE;
    Trace_Mode trace_state = Trace_Mode::GEODESIC;

    int traceIdx;
    int traceSign;
    int traceFaceId;
    int traceSteps;
    
    bool showBending;
    bool showSingularities;
    Eigen::MatrixXd singularVerts_topo;
    Eigen::MatrixXd singularVerts_geo;
    Eigen::MatrixXd nonIdentityEdges;

    std::string vectorFieldName;
};

#endif
