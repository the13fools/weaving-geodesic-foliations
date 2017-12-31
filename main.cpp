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

double operator_scale = 1.;
int opt_step = 0;

void logToFile(const Eigen::MatrixXd W, std::string foldername, std::string filename)
{
#ifndef WIN32
    char folderpath[50];
    sprintf(folderpath, "log/%s", foldername.c_str());
    mkdir("log", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(folderpath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    char logpath[50];
    sprintf(logpath, "%s/%s.txt", folderpath, filename.c_str());   
    std::ofstream myfile (logpath);
    for(int i = 0; i < W.rows(); i++)
    {
        if (myfile.is_open())
        {
	    myfile << W.row(i) << "\n";
	}

	else
	{
	    std::cout << "Unable to open file";
	    break;
	}
    } 
    myfile.close();
#endif
}

int descentStep = 0;
void takeGradientDescentStep()
{
    int nfaces = curMesh->F.rows();        

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
    int nfaces = (int)curMesh->F.rows();
    for (int i = 0; i < nfaces; i++)
    {
        curMesh->v0.row(i) = projectOntoFace(curMesh->v0.row(i), curMesh->F, curMesh->V, i);
    }
      
    initOptVars(*curMesh, curMesh->v0, curMesh->optVars);

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
        viewer.ngui->addWindow(Eigen::Vector2i(10, 60), "Weaving");
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

              // call to generate menu
        viewer.screen->performLayout();
        return false;
    };

    viewer->launch();
}
