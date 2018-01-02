#include <igl/viewer/Viewer.h>
#include <igl/avg_edge_length.h>

#include <math.h>
#include <iostream>
#include <fstream>

// For making directories
#include <sys/stat.h>
// #include <direct.h>

#include "DataLoad.h"
#include "Covariant.h"
#include "FaceBased.h"
#include "FieldOptimization.h"
#include "RelaxViewer.h"
MeshData *curMesh;

igl::viewer::Viewer *viewer;
Weights w;

double px = 2;
double py = 0;
int desc_steps = 1; // Inner Loop
int desc_loops = 1; // Outer Loop
char fileName[50] = "0-0-10ksteps";
std::string folderName = "logging_location";
std::string fieldName = "field_dt.txt";

int idx0 = 890;
int idx1 = 4064;
int idx2 = 4537;
int idx3 = 9243;

double lambda = 10.;

double px0 = 10.;
double py0 = 0.;
double px1 = 10.;
double py1 = 0.;
double px2 = 0.;
double py2 = 10.;
double px3 = 0.;
double py3 = 10.;

double operator_scale = 1.;
int opt_step = 0;




int descentStep = 0;
void takeGradientDescentStep()
{
    int nfaces = curMesh->F.rows();
  //  Weights w;
    w.handleWeights.resize(nfaces);
    w.handleWeights.setConstant(1.0);
    w.lambdaDreg = 1;
    w.lambdaGeodesic = lambda;
    w.lambdaVD = lambda;
    w.lambdaVW = lambda;

    for (int loops = 0; loops < desc_loops; loops++) {

        Eigen::MatrixXd Op_Grad(curMesh->F.rows(), 3);
        double t = 1.0;
        Eigen::MatrixXd del_W_V;
        // Not effecient, but will make it feel more correct to update, then show
        for (int i = 0; i < desc_steps; i++)
        {
            alternatingMinimization(*curMesh, w, curMesh->optVars);
    	}
        descentStep++;
        updateView(curMesh, viewer);
    }
}

void showVectorField()
{
    Eigen::Vector3d p(px, py,0);
    computeDistanceField(p, curMesh->centroids_F, curMesh->v0);
    initOptVars(*curMesh, curMesh->v0, curMesh->optVars);

//    computeDistanceField(p, centroids_F, W_init);
//    computeWhirlpool(p, centroids_F, W_init);
//    computeWhirlpool(p, centroids_F, W);
//    W_init = W;
/*    W = Eigen::MatrixXd::Zero(W_init.rows(), W_init.cols());
    for (int i = 0; i < W_init.rows(); i++) 
    {
        curMesh->v0.row(i) = projectOntoFace(curMesh->v0.row(i), curMesh->F, curMesh->V, i);
    }
      
    initOptVars(*curMesh, curMesh->v0, curMesh->optVars);
*/
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
    updateView(curMesh, viewer);
}


int main(int argc, char *argv[])
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    if (!igl::readOBJ("../sphere.obj", V, F))
        return -1;
    curMesh = new MeshData(V, F);

    w.handleWeights.resize(F.rows());
    w.handleWeights.setConstant(1.0);
    w.lambdaDreg = 1e-4;
    w.lambdaGeodesic = 1000;
    w.lambdaVD = 1000;
    w.lambdaVW = 1000;
    w.lambdaunit = 10;

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

      viewer.ngui->addWindow(Eigen::Vector2i(-10, 400), "Handles");
      viewer.ngui->addVariable("Norm Vectors", curMesh->vs.normFaceVectors);
      viewer.ngui->addVariable("lambda", lambda);
      viewer.ngui->addVariable("idx0", idx0);
      viewer.ngui->addVariable("px",px0);
      viewer.ngui->addVariable("py",py0);
      viewer.ngui->addVariable("idx1", idx1);
      viewer.ngui->addVariable("px",px1);
      viewer.ngui->addVariable("py",py1);
      viewer.ngui->addVariable("idx0", idx2);
      viewer.ngui->addVariable("px",px2);
      viewer.ngui->addVariable("py",py2);
      viewer.ngui->addVariable("idx0", idx3);
      viewer.ngui->addVariable("px",px3);
      viewer.ngui->addVariable("py",py3);

              // call to generate menu
        viewer.screen->performLayout();
        return false;
    };

    viewer->launch();
}
