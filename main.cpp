#include <igl/viewer/Viewer.h>

#include <math.h>
#include <iostream>
#include <fstream>

// For making directories
#include <sys/stat.h>
// #include <direct.h>

#include "DataLoad.h"
#include "Covariant.h"
#include "FaceBased.h"
#include "StateStructs.h"
#include "RelaxViewer.h"


// stuff that is precomputed about the mesh and doesn't change during optimization
MeshData curMesh;

// stuff that is updated through the course of optimization.  Mostly a convenience to make visualization easier.
// Still good to unpack for most functions to maintain clarity + compile time checks  
OptState optState;


igl::viewer::Viewer *viewer;

double px = 0;
double py = 0;
int desc_steps = 1; // Inner Loop
int desc_loops = 1; // Outer Loop
char fileName[50] = "0-0-10ksteps";
std::string folderName = "logging_location";
std::string fieldName = "field_dt.txt";

double operator_scale = 1.;

void computeDelWV(const MeshData &mesh, 
        const Eigen::MatrixXd &v, 
	    const Eigen::MatrixXd &w,
        Eigen::MatrixXd &result)
{
    result.resize(v.rows(), v.cols());
    int nfaces = mesh.F.rows();
    for (int i = 0; i < nfaces; i++)
    {
        result.row(i) = w.row(i) * (mesh.Ms[i] * v);
    }
}

// energy function
double energy(const MeshData &mesh,
    const Eigen::MatrixXd &v, const Eigen::MatrixXd &w)
{
    Eigen::MatrixXd delwv;
    computeDelWV(mesh, v, w, delwv);

    int nfaces = (int)v.rows();
    double result = 0;
    for(int i=0; i<nfaces; i++)
        result += delwv.row(i).dot(delwv.row(i));
    return result / 2.0;
}

// energy derivative with respect to v
void dvEnergy(const MeshData &mesh,  
    const Eigen::MatrixXd &v,
    const Eigen::MatrixXd &w,
    Eigen::MatrixXd &dE)
{
    int nfaces = v.rows();

    dE.resize(nfaces, 3);
    dE.setZero();

    Eigen::MatrixXd delwv;
    computeDelWV(mesh, v, w, delwv);

    for (int i = 0; i < nfaces; i++)
    {
        Eigen::VectorXd temp = mesh.Ms[i].transpose() * w.row(i).transpose(); 
        dE += temp * delwv.row(i);
    }
}

// energy derivative with respect to w
void dwEnergy(const MeshData &mesh,
    const Eigen::MatrixXd &v,
    const Eigen::MatrixXd &w,
    Eigen::MatrixXd &dE)
{
    int nfaces = v.rows();

    dE.resize(nfaces, 3);
    dE.setZero();

    Eigen::MatrixXd delwv;
    computeDelWV(mesh, v, w, delwv);

    for (int i = 0; i < nfaces; i++)
    {
        Eigen::VectorXd rowTemp = delwv.row(i) * v.transpose();
        dE.row(i) += rowTemp.transpose() * mesh.Ms[i].transpose();
    }
}

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


void computeOperatorGradient_finitedifference(const Eigen::MatrixXd &v, Eigen::MatrixXd &Op_Grad)
{
    int nfaces = v.rows();

    double e = energy(curMesh, v,v); 

    Op_Grad = Eigen::MatrixXd::Zero(nfaces, 3);
    double eps = .0000001;

    for (int i = 0; i < nfaces; i++)
    {
        for (int j = 0; j < 3; j++) 
	{
            Eigen::MatrixXd shifted = v;
	    shifted(i, j) += eps;
	    double diff = energy(curMesh, shifted, shifted);
	    Op_Grad(i, j) = (diff - e) / eps;
	}	
	std::cout << i << "\n";
    }

}


Eigen::MatrixXd Op_Grad_fd;
bool is_fd_set;
void computeFiniteDifference() 
{
    is_fd_set = true;
    computeOperatorGradient_finitedifference(optState.W, Op_Grad_fd);
    logToFile(Op_Grad_fd, folderName, "fd2");
}

void loadFiniteDifference()
{
    Op_Grad_fd = readMatrix("fd.txt");
}

void lineSearch(const Eigen::MatrixXd &v, const Eigen::MatrixXd &gradV, double &t, double &newenergy)
{
    double c1 = 1.5;
    double c2 = 0.5;
    t *= c1;
    double orig = energy(curMesh, v, v);
    std::cout << "original energy " << orig << std::endl;
    while (true)
    {
        Eigen::MatrixXd testv = v - t*gradV;
        newenergy = energy(curMesh, testv, testv);
        std::cout << "new energy, t = " << t << ": " << newenergy << std::endl;
        if (newenergy < orig || t < .00000001)
            return;
        else
            t *= c2;
    }
}


int descentStep = 0;
void takeGradientDescentStep()
{
    for (int loops = 0; loops < desc_loops; loops++) {
        int nfaces = curMesh.F.rows();

        Eigen::MatrixXd Op_Grad(curMesh.F.rows(), 3);
        double t = 1.0;
        Eigen::MatrixXd del_W_V;
        // Not effecient, but will make it feel more correct to update, then show
        for (int i = 0; i < desc_steps; i++)
        {
            Op_Grad.setZero();
            computeDelWV(curMesh, optState.W, optState.W, del_W_V);
            Eigen::MatrixXd dE(curMesh.F.rows(), 3);
            dvEnergy(curMesh, optState.W, optState.W, dE);
            Op_Grad += dE;
            dwEnergy(curMesh, optState.W, optState.W, dE);
            Op_Grad += dE;

            double newenergy = 0;
            lineSearch(optState.W, Op_Grad, t, newenergy);
            optState.W -= t*Op_Grad;
        }
        computeDelWV(curMesh, optState.W, optState.W, del_W_V);
        descentStep++;
        updateView(curMesh, optState, viewer);
    }
}

void loadField()
{
    optState.W = readMatrix(fieldName.c_str());
    Eigen::MatrixXd del_W_V;
    computeDelWV(curMesh, optState.W, optState.W, del_W_V);

    descentStep = 1;
    updateView(curMesh, optState, viewer);
}

void showVectorField()
{
    computeCentroids(curMesh.F,curMesh.V,curMesh.centroids_F);

    Eigen::Vector3d p(px, py,0);
    computeDistanceField(p, curMesh.centroids_F, optState.W);

    computeDistanceField(p, curMesh.centroids_F, curMesh.W_init);
//    computeWhirlpool(p, centroids_F, W_init);
//    computeWhirlpool(p, centroids_F, W);
//    W_init = W;
    
    Eigen::MatrixXd del_W_V;
    computeDelWV(curMesh, optState.W, optState.W, del_W_V);

//    project(W_init);
    logToFile(curMesh.W_init, folderName, "0W");
    descentStep = 1;
    
    updateView(curMesh, optState, viewer);
}


void addNoiseToField() 
{
    double eps = .1;
    Eigen::MatrixXd noise = Eigen::MatrixXd::Random( optState.W.rows(), 3 ) * eps;
    for (int i = 0; i < optState.W.rows(); i++)
    {
        noise(i, 2) = 0.;
    }

    optState.W += noise;
  //  W_init += noise; // This way we see the error from the ground truth

    Eigen::MatrixXd del_W_V;
    computeDelWV(curMesh, optState.W, optState.W, del_W_V);

    descentStep = 1;
    updateView(curMesh, optState, viewer);
}

void testGradients()
{
    int triangleToTest = 100;
    double energyorig = energy(curMesh, optState.W, optState.W);
    for (int i = 0; i < 2; i++)
    {
        Eigen::MatrixXd perturbed = optState.W;
        perturbed(triangleToTest, i) += 1e-6;
        double energynewv = energy(curMesh, perturbed, optState.W);
        double energyneww = energy(curMesh, optState.W, perturbed);
        double findiffv = (energynewv - energyorig) / 1e-6;
        double findiffw = (energyneww - energyorig) / 1e-6;
        Eigen::MatrixXd OpGradv;
        Eigen::MatrixXd OpGradw;
        dvEnergy(curMesh, optState.W, optState.W, OpGradv);
        dwEnergy(curMesh, optState.W, optState.W, OpGradw);
        std::cout << "v: " << findiffv << " vs " << OpGradv(triangleToTest, i) << "    w: " << findiffw << " vs " << OpGradw(triangleToTest, i) << std::endl;                   
    }
}



int main(int argc, char *argv[])
{  
  //   assignFaceVal(F,viz);;

  igl::readOBJ("../circ.obj", curMesh.V, curMesh.F);
  buildEdges(curMesh.F, curMesh.E);
  buildEdgesPerFace(curMesh.F, curMesh.E, curMesh.F_edges);
  computeGradientMatrices(curMesh.F, curMesh.V, curMesh.E, curMesh.F_edges, curMesh.Ms);

  // Plot the mesh  
  viewer = new igl::viewer::Viewer();
  viewer->data.set_mesh(curMesh.V, curMesh.F);
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
      viewer.ngui->addButton("Compute Finite Diff", computeFiniteDifference);
      viewer.ngui->addButton("Test Gradients", testGradients);
      viewer.ngui->addButton("Recompute Derivative", showVectorField);
      viewer.ngui->addButton("Grad Descent Step", takeGradientDescentStep);

      viewer.ngui->addVariable("Log Folder", folderName);
//      viewer.ngui->addVariable("Shade State", shading_enum_state, true)
//          ->setItems({"Operator Error", "Update Direction", "Update Magnitude"});
//	  ->setCallback([] { updateView(del_W_F, descentStep); });


      viewer.ngui->addVariable("Load Field Name", fieldName);
      viewer.ngui->addButton("Load Field", loadField);

      // call to generate menu
      viewer.screen->performLayout();
      return false;
  };

  viewer->launch();
}
