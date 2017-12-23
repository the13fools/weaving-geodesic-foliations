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


Eigen::MatrixXd V;
Eigen::MatrixXi F, E, F_edges;

igl::viewer::Viewer *viewer;

double px = 0;
double py = 0;
int desc_steps = 10;
int desc_loops = 3;
char fileName[50] = "0-0-10ksteps";
std::string folderName = "logging_location";

enum shading_enum {
   OP_ERROR = 0,
   INIT_DIRECTION,
   INIT_MAGNITUDE
};
shading_enum shading_enum_state = INIT_DIRECTION;

void computeCovariantOperator(const Eigen::VectorXd &scalar_F,
    const Eigen::MatrixXi &F,
    const Eigen::MatrixXi &F_edges,
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &E,
    Eigen::MatrixXd &DM_local)
{
    Eigen::VectorXd scalar_E;
    computeEdgeWeights(scalar_F, V, E, scalar_E); 
    
    int nfaces = F.rows();

    DM_local.resize(nfaces, 2);

    for (int i = 0; i < nfaces; i++)
    {
        Eigen::Vector3d u = V.row(F(i, 1)) - V.row(F(i, 0));
        Eigen::Vector3d v = V.row(F(i, 2)) - V.row(F(i, 0));

        double u_weight = 2 * (scalar_E(F_edges(i, 0)) - scalar_E(F_edges(i, 1)));
        double v_weight = 2 * (scalar_E(F_edges(i, 0)) - scalar_E(F_edges(i, 2)));

        DM_local(i, 0) = u_weight;
	DM_local(i, 1) = v_weight; 
    }
}

void evaluateOperator(const Eigen::MatrixXd &DM_local,
	const Eigen::MatrixXd &W_local, 
        Eigen::MatrixXd &del_W_F, int idx)
{
    Eigen::Vector3d component(0, 0, 0);
    component(idx) = 1;

// This is super slow because it does a bunch of irrelevant multiplications...    
//    Eigen::MatrixXd weights = DM_local * W_local.transpose();

    for (int i = 0; i < DM_local.rows(); i++)
    {
 
 	del_W_F.row(i) += component * (DM_local(i, 0) * W_local(i, 0) + DM_local(i,1) * W_local(i, 1));
 
    }
}

void computeCovariantOperatorNew(const Eigen::MatrixXi &F, 
	const Eigen::MatrixXd &V, 
	const Eigen::MatrixXi &E, 
	const Eigen::MatrixXi &F_edges, 
        const Eigen::MatrixXd &v, 
	const Eigen::MatrixXd &w,
        const std::vector<Eigen::SparseMatrix<double> > &Ms,
        Eigen::MatrixXd &result)
{
    result.resize(v.rows(), v.cols());
    int nfaces = F.rows();
    for (int i = 0; i < nfaces; i++)
    {
        result.row(i) = w.row(i) * (Ms[i] * v);
    }
}

Eigen::MatrixXd W; // This can be thought of as 3 ``independent'' scalar fields
Eigen::MatrixXd W_init; // W at init, for visualizing change in descent 

std::vector<Eigen::SparseMatrix<double> > Ms; // the gradient operator; Ms[i] * F gives the gradient of F on triangle i

Eigen::MatrixXd colorField;
Eigen::MatrixXd centroids_F;

Eigen::MatrixXd del_W_F;




void computeOperatorGradient( const std::vector<Eigen::SparseMatrix<double> > &Ms,
	const Eigen::MatrixXd &del_W_F,
	const Eigen::MatrixXd &v,
	Eigen::MatrixXd &Op_Grad)
{
    int nfaces = v.rows();

    Op_Grad = Eigen::MatrixXd::Zero(nfaces, 3);

    for (int i = 0; i < nfaces; i++)
    {
	Eigen::VectorXd temp = Ms[i].transpose() * v.row(i).transpose(); 
	Op_Grad += temp * del_W_F.row(i);

        Eigen::VectorXd rowTemp = del_W_F.row(i) * v.transpose();
	Op_Grad.row(i) += rowTemp * Ms[i].transpose();
    }

}

void logToFile(const Eigen::MatrixXd W, std::string foldername, std::string filename)
{
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
}

double energy_OP = 0.;
void updateView(const Eigen::MatrixXd del_W_F, int step)
{
    int nFaces = F.rows(); 

    logToFile(W, folderName, std::to_string(step)); 

    // Set mesh colors and log operator state
    Eigen::VectorXd Z(nFaces);
    Eigen::Vector3d testDir(1,0,0);
    energy_OP = 0.;
    
    for (int i = 0; i < nFaces; i++)
    {
        switch (shading_enum_state)
	{
            case OP_ERROR:
                Z(i) = log(del_W_F.row(i).norm() + .000005);
		break;
            case INIT_DIRECTION:
		Z(i) = (W-W_init).row(i).normalized().dot(testDir);
		break;
            case INIT_MAGNITUDE:
		Z(i) = log( (W-W_init).row(i).squaredNorm() );
		break;
	}

	energy_OP += del_W_F.row(i).squaredNorm();       
    }
    std::cout << energy_OP << " Operator Energy\n";

    // Average edge length for sizing
    const double avg = igl::avg_edge_length(V,F);
    colorField.resize(nFaces, 3);
    
    //  igl::jet(Z,true,colorField);
    igl::colormap(igl::COLOR_MAP_TYPE_MAGMA,Z, true, colorField);


    // Plot the mesh
    viewer->data.clear();
    viewer->data.set_mesh(V, F);
    viewer->data.set_face_based(true);

    viewer->data.set_colors(colorField);

    Eigen::MatrixXd eps = Eigen::MatrixXd::Constant(nFaces,3,.001);
/*
    Eigen::MatrixXd field;
    field.resize(nFaces, 3);
    for (int i = 0; i < nFaces; i++)
    {
	field.row(i) = W.row(i).normalized();
    }
*/
    const Eigen::RowVector3d green(0.1,0.9,0.2),blue(0.1,0.2,0.8);
  //  viewer->data.add_edges(centroids_f  + del_w_f*avg/2, centroids_f, blue);
    viewer->data.add_edges(centroids_F  + W*avg/2, centroids_F, blue);
    viewer->data.add_edges(centroids_F  + (W - W_init)*avg/2, centroids_F, green);
}


void project(Eigen::MatrixXd &v)
{
	int nfaces = (int)v.rows();
	double fieldMass = 0;
	for (int i = 0; i < nfaces; i++)
	{
	    fieldMass += v.row(i).squaredNorm();
	} 

	v *= (double) sqrt(nfaces / fieldMass);	
}

double energy(const Eigen::MatrixXd &v)
{
	int nfaces = (int)v.rows();
	del_W_F(nfaces, 3);
	del_W_F.setZero();
	computeCovariantOperatorNew(F, V, E, F_edges, v, v, Ms, del_W_F);
	double result = 0;
	for(int i=0; i<nfaces; i++)
		result += del_W_F.row(i).dot(del_W_F.row(i));
	return result;
}

void lineSearch(const Eigen::MatrixXd &v, const Eigen::MatrixXd &gradV, double &t, double &newenergy)
{
    double c1 = 1.5;
    double c2 = 0.5;
    t *= c1;
    double orig = energy(v);
    std::cout << "original energy " << orig << std::endl;
    Eigen::MatrixXd v_cp = v;
    project(v_cp);
    double projected = energy(v_cp);
    std::cout << "projected energy " << projected << std::endl;
    while(true)
    {
    	Eigen::MatrixXd testv = v - t*gradV;    
    	project(testv);
    	newenergy = energy(testv);
	std::cout << "new energy, t = " << t << ": " << newenergy << std::endl;
	if(newenergy < orig || t < .0000000000001)
		return;
	else
		t *= c2;
    }
}



int descentStep = 0;
void takeGradientDescentStep()
{
for(int loops = 0; loops < desc_loops; loops++) {
    int nfaces = F.rows();

    del_W_F.resize(nfaces, 3);
    del_W_F.setZero();
    Eigen::MatrixXd Op_Grad;
    double t = 1.0;

    // Not effecient, but will make it feel more correct to update, then show
    for (int i = 0; i < desc_steps; i++) 
    {
        computeCovariantOperatorNew(F, V, E, F_edges, W, W, Ms, del_W_F);
        computeOperatorGradient(Ms, del_W_F, W, Op_Grad);
	double newenergy = 0;
	lineSearch(W, Op_Grad, t, newenergy);
        W -= t*Op_Grad; 
	project(W);
    }

    computeCovariantOperatorNew(F, V, E, F_edges, W, W, Ms, del_W_F);
       
    descentStep++;
    updateView(del_W_F, descentStep);
}}


void showVectorField()
{
    computeCentroids(F,V,centroids_F);

    Eigen::Vector3d p(px, py,0);
    computeDistanceField(p, centroids_F, W);
    computeDistanceField(p, centroids_F, W_init);
//    computeWhirlpool(p, centroids_F, W);
//    W_init = W;


    computeCovariantOperatorNew(F, V, E, F_edges, W, W, Ms, del_W_F);

    int nFaces = F.rows(); 

    //  Eigen::VectorXd Z = W.col(0); // - del_W_F;// - W_recovered.col(0);
    // Eigen::VectorXd Z = del_W_F.transpose() * del_W_F;// - W_recovered.col(0);

    project(W_init);
    logToFile(W_init, folderName, "0W");
    descentStep = 1;
    updateView(del_W_F, descentStep);   

}



void addNoiseToField() 
{
    double eps = .1;
    Eigen::MatrixXd noise = Eigen::MatrixXd::Random( W.rows(), 3 ) * eps;
    for (int i = 0; i < W.rows(); i++)
    {
        noise(i, 2) = 0.;
    }

    W += noise;
  //  W_init += noise; // This way we see the error from the ground truth

    project(W);

    computeCovariantOperatorNew(F, V, E, F_edges, W, W, Ms, del_W_F);

    descentStep = 1;
    updateView(del_W_F, descentStep);
}



int main(int argc, char *argv[])
{  
  //   assignFaceVal(F,viz);;

  igl::readOBJ("../circ.obj", V, F);
  buildEdges(F, E);
  buildEdgesPerFace(F, E, F_edges);
  computeGradientMatrices(F, V, E, F_edges, Ms);

  // Plot the mesh  
  viewer = new igl::viewer::Viewer();
  viewer->data.set_mesh(V, F);
  viewer->data.set_face_based(true);
  viewer->callback_init = [&](igl::viewer::Viewer& viewer)
  {
      
      viewer.ngui->window()->setVisible(false);
      viewer.ngui->addWindow(Eigen::Vector2i(10, 60), "Weaving"); 
      // Add new group
      viewer.ngui->addGroup("Vector Field Options");

      // Expose a variable
      viewer.ngui->addVariable("Center X",px);
      viewer.ngui->addVariable("Center Y",py);
      viewer.ngui->addVariable("Descent Steps",desc_steps);
      viewer.ngui->addVariable("Descent Outer Loops",desc_loops);

      // Add a button
      viewer.ngui->addButton("Recompute Derivative", showVectorField);
      viewer.ngui->addButton("Add Noise to Field", addNoiseToField);
      viewer.ngui->addButton("Grad Descent Step", takeGradientDescentStep);

      viewer.ngui->addVariable("Log Folder", folderName);
      viewer.ngui->addVariable("Shade State", shading_enum_state, true)
          ->setItems({"Operator Error", "Update Direction", "Update Magnitude"});


      // call to generate menu
      viewer.screen->performLayout();
      return false;
  };

  viewer->launch();
}
