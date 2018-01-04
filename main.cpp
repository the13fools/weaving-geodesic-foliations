#include <igl/viewer/Viewer.h>
#include <igl/avg_edge_length.h>

#include <math.h>
#include <iostream>
#include <fstream>

// For making directories
#include <sys/stat.h>
// #include <direct.h>

#include "InitField.h"
#include "Covariant.h"
#include "FaceBased.h"
#include "FieldOptimization.h"
#include "RelaxViewer.h"
#include "Physics.h"
#include "VectorUtils.h"

MeshData *curMesh;
PhysicsData phydata;

igl::viewer::Viewer *viewer;
Weights w;

double px = 2;
double py = 0;
int desc_steps = 1; // Inner Loop
int desc_loops = 1; // Outer Loop
char fileName[50] = "0-0-10ksteps";
std::string folderName = "logging_location";
std::string fieldName = "field_dt.txt";

int idx0 = 0;
int idx1 = 1;
int idx2 = 2;
int idx3 = 3;

double pu0 = 1.;
double pv0 = 1.;
double pu1 = 1.;
double pv1 = 1.;
double pu2 = 1.;
double pv2 = 1.;
double pu3 = 1.;
double pv3 = 1.;

double operator_scale = 1.;
int opt_step = 0;


Eigen::MatrixXd forceField; // temp global for now

void setHandleWeights(Weights &weight)
{
    int nfaces = curMesh->F.rows();
    w.handleWeights.resize(nfaces);
    w.handleWeights.setConstant(0.0);
    
    int idx[] = {idx0, idx1, idx2, idx3};
    double pu[] = {pu0, pu1, pu2, pu3};
    double pv[] = {pv0, pv1, pv2, pv3};
       
    for (int i = 0; i < 4; i++)
    {
        w.handleWeights(idx[i]) = 1.;
        Eigen::Vector3d u = curMesh->V.row(curMesh->F(idx[i], 0))
            - curMesh->V.row(curMesh->F(idx[i], 1));
        Eigen::Vector3d n = faceNormal(curMesh->F, curMesh->V, idx[i]);
        u.normalize();
        Eigen::Vector3d v = u.cross(n);
        curMesh->v0.row(idx[i]) = u * pu[i] + v * pv[i];
    }
}


int descentStep = 0;
void takeGradientDescentStep()
{
    int nfaces = curMesh->F.rows();
    setHandleWeights(w);

    for (int loops = 0; loops < desc_loops; loops++) {

        Eigen::MatrixXd Op_Grad(curMesh->F.rows(), 3);
        double t = 1.0;
        Eigen::MatrixXd del_W_V;
        // Not effecient, but will make it feel more correct to update, then show
        for (int i = 0; i < desc_steps; i++)
        {
            alternatingMinimization(*curMesh, curMesh->vs, w, curMesh->optVars);
    	}
        descentStep++;
    }
    Eigen::VectorXd wf(curMesh->F.rows());
    wf.setConstant(1.0);
    physicalForces(*curMesh, phydata, wf, curMesh->optVars.vbar, 1.0, 1.0 / 3.0, forceField);    
    updateView(curMesh, viewer);
}

void showVectorField()
{
    curMesh->v0 = Eigen::MatrixXd::Zero(curMesh->F.rows(), 3);
    setHandleWeights(w);
    propogateField(curMesh->F, curMesh->V, curMesh->E, curMesh->F_edges, curMesh->v0);
    initOptVars(*curMesh, curMesh->v0, curMesh->optVars);
    
    Eigen::VectorXd wf(curMesh->F.rows());
    wf.setConstant(1.0);
    physicalForces(*curMesh, phydata, wf, curMesh->optVars.vbar, 1.0, 1.0 / 3.0, forceField);    

    // test the force code using finite difference
    for (int i = 0; i < 3; i++)
    {
        // test bend force
        testForces(*curMesh, phydata, wf, curMesh->optVars.vbar, 100.0, 0.0, 100, i);
        // test twist force
        testForces(*curMesh, phydata, wf, curMesh->optVars.vbar, 0.0, 100.0, 100, i);
    }
    
    energy(curMesh->optVars, *curMesh, w, curMesh->vs);

    descentStep = 1;
    updateView(curMesh, viewer);

}


void addNoiseToField() 
{
    int nvectors = curMesh->v0.rows();
    if (nvectors == 0)
        return;
    double eps = .1;
    Eigen::MatrixXd noise = Eigen::MatrixXd::Random( nvectors, 3 ) * eps;
    curMesh->v0 += noise;
    for (int i = 0; i < nvectors; i++)
    {
        curMesh->v0.row(i) = projectOntoFace(curMesh->v0.row(i), curMesh->F, curMesh->V, i);
    }
    initOptVars(*curMesh, curMesh->v0, curMesh->optVars);
  
    descentStep = 1;
    Eigen::VectorXd wf(curMesh->F.rows());
    wf.setConstant(1.0);
    physicalForces(*curMesh, phydata, wf, curMesh->optVars.vbar, 1.0, 1.0 / 3.0, forceField);    

    updateView(curMesh, viewer);
}


int main(int argc, char *argv[])
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    if (!igl::readOBJ("../models/torus.obj", V, F))
        return -1;
    curMesh = new MeshData(V, F);
    physicsDataFromMesh(*curMesh, phydata);

    w.handleWeights.resize(F.rows());
    w.handleWeights.setConstant(1.0);
    w.lambdaDreg = 1e-4;
    w.lambdaGeodesic = 1000;
    w.lambdaVD = 1000;
    w.lambdaVW = 1000;
    w.lambdaunit = 10;

    curMesh->vs.physicsEnergyArrowScale = 1;

    // Plot the mesh  
    viewer = new igl::viewer::Viewer();
    viewer->data.set_mesh(curMesh->V, curMesh->F);
    viewer->data.set_face_based(true);
    viewer->callback_init = [&](igl::viewer::Viewer& viewer)
    {

        viewer.ngui->window()->setVisible(false);
        viewer.ngui->addWindow(Eigen::Vector2i(10, 10), "Weaving");
        // Add new group
        viewer.ngui->addGroup("Vector Field Options");

        // Expose a variable
        viewer.ngui->addVariable("Operator Scale", operator_scale);
        viewer.ngui->addVariable("Center X", px);
        viewer.ngui->addVariable("Center Y", py);
        viewer.ngui->addVariable("Descent Steps", desc_steps);
        viewer.ngui->addVariable("Descent Outer Loops", desc_loops);

        // Add a button
        viewer.ngui->addButton("Add Noise to Field", addNoiseToField);
        viewer.ngui->addButton("Compute Initial Field", showVectorField);
        viewer.ngui->addButton("Grad Descent Step", takeGradientDescentStep);

        viewer.ngui->addVariable("Log Folder", folderName);


        viewer.ngui->addGroup("Weights");
        viewer.ngui->addVariable("Geodesic", w.lambdaGeodesic);
        viewer.ngui->addVariable("VW", w.lambdaVW);
        viewer.ngui->addVariable("D Compatibility", w.lambdaVD);
        viewer.ngui->addVariable("D Regularization", w.lambdaDreg);
        viewer.ngui->addVariable("Unit Length", w.lambdaunit);

      viewer.ngui->addVariable("Norm Vectors", curMesh->vs.physicsEnergyArrowScale);


      viewer.ngui->addWindow(Eigen::Vector2i(1000, 10), "Handles");
      viewer.ngui->addVariable("Norm Vectors", curMesh->vs.normFaceVectors);
      viewer.ngui->addVariable("idx0", idx0);
      viewer.ngui->addVariable("pu",pu0);
      viewer.ngui->addVariable("pv",pv0);
      viewer.ngui->addVariable("idx1", idx1);
      viewer.ngui->addVariable("pu",pu1);
      viewer.ngui->addVariable("pv",pv1);
      viewer.ngui->addVariable("idx2", idx2);
      viewer.ngui->addVariable("pu",pu2);
      viewer.ngui->addVariable("pv",pv2);
      viewer.ngui->addVariable("idx3", idx3);
      viewer.ngui->addVariable("pu",pu3);
      viewer.ngui->addVariable("pv",pv3);

              // call to generate menu
        viewer.screen->performLayout();
        return false;
    };

    viewer->launch();
}
