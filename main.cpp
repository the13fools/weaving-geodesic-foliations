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

double px = 0;
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
    Weights w;
    w.handleWeights.resize(nfaces);
    w.handleWeights.setConstant(1.0);
    w.lambdaDreg = 1;
    w.lambdaGeodesic = 1000;
    w.lambdaVD = 1000;
    w.lambdaVW = 1000;

    for (int loops = 0; loops < desc_loops; loops++) {

        Eigen::MatrixXd Op_Grad(curMesh->F.rows(), 3);
        double t = 1.0;
        Eigen::MatrixXd del_W_V;
        // Not effecient, but will make it feel more correct to update, then show
        for (int i = 0; i < desc_steps; i++)
        {
/*            Op_Grad.setZero();
            computeDelWV(*curMesh, W, W, del_W_V);
            Eigen::MatrixXd dE(curMesh->F.rows(), 3);
            dvEnergy(*curMesh, W, W, dE);
            Op_Grad += dE;
            dwEnergy(*curMesh, W, W, dE);
            Op_Grad += dE;

            double newenergy = 0;
            lineSearch(W, Op_Grad, t, newenergy);
            W -= t*Op_Grad;
  */
	
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
    initOptVars(curMesh->v0, curMesh->optVars);

//    computeDistanceField(p, centroids_F, W_init);
//    computeWhirlpool(p, centroids_F, W_init);
//    computeWhirlpool(p, centroids_F, W);
//    W_init = W;
/*    W = Eigen::MatrixXd::Zero(W_init.rows(), W_init.cols());
    for (int i = 0; i < W_init.rows(); i++) 
    {
	for (int j = 0; j < 3; j++)
	{
	    const Eigen::Vector4i &einfo = E.row(F_edges(i,j));
	    // if face faceidx is on the boundary, return empty matrix
	    if (einfo[2] == -1 || einfo[3] == -1)
	    {
		W.row(i) = W_init.row(i);
		break;
	    }
	    else 
	    {
	       W.row(i) = Eigen::VectorXd::Random(3) * .00001;
	    }
	}

    }
*/
    Weights w;
    int nfaces = (int)curMesh->F.rows();
    w.handleWeights.resize(nfaces);
    w.handleWeights.setConstant(1.0);
    w.lambdaDreg = 1;
    w.lambdaGeodesic = 1000;
    w.lambdaVD = 1000;
    w.lambdaVW = 1000;
    alternatingMinimization(*curMesh, w, curMesh->optVars);

    descentStep = 1;
    updateView(curMesh, viewer);

}


void addNoiseToField() 
{
    double eps = .1;
    Eigen::MatrixXd noise = Eigen::MatrixXd::Random( curMesh->v0.rows(), 3 ) * eps;
    for (int i = 0; i < curMesh->v0.rows(); i++)
    {
        noise(i, 2) = 0.;
    }

    curMesh->v0 += noise;
    initOptVars(curMesh->v0, curMesh->optVars);
  
    descentStep = 1;
    updateView(curMesh, viewer);
}


int main(int argc, char *argv[])
{  
  //   assignFaceVal(F,viz);;
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
  igl::readOBJ("../circ.obj", V, F);
  curMesh = new MeshData(V, F);
  
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
      viewer.ngui->addVariable("Operator Scale",operator_scale);
      viewer.ngui->addVariable("Center X",px);
      viewer.ngui->addVariable("Center Y",py);
      viewer.ngui->addVariable("Descent Steps",desc_steps);
      viewer.ngui->addVariable("Descent Outer Loops",desc_loops);

      // Add a button
      viewer.ngui->addButton("Add Noise to Field", addNoiseToField);
      viewer.ngui->addButton("Recompute Derivative", showVectorField);
      viewer.ngui->addButton("Grad Descent Step", takeGradientDescentStep);

      viewer.ngui->addVariable("Log Folder", folderName);
//      viewer.ngui->addVariable("Shade State", shading_enum_state, true)
//          ->setItems({"Operator Error", "Update Direction", "Update Magnitude"});
//	  ->setCallback([] { updateView(del_W_F, descentStep); });


      // call to generate menu
      viewer.screen->performLayout();
      return false;
  };

  viewer->launch();
}
