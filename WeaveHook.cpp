#include "WeaveHook.h"
#include "GaussNewton.h"
#include <iostream>

using namespace std;


void WeaveHook::setFaceColors(igl::viewer::Viewer &viewer)
{ 
    int faces = weave->F.rows();
    if ( curFaceEnergies.rows() != faces && shading_state != NONE) { return ; }
 
    igl::ColorMapType viz_color = igl::COLOR_MAP_TYPE_MAGMA;
    
    Eigen::VectorXd Z(faces);
    
    for (int i = 0; i < faces; i++) 
    {
	switch(shading_state)
	{
	    case NONE:
		Z(i) = .7;
		break;
	    case F1_ENERGY:
		std::cout << curFaceEnergies(i,0) << std::endl;
		Z(i) = log(curFaceEnergies(i,0));
	        break;	
	    case F2_ENERGY:
		Z(i) = log(curFaceEnergies(i,1));
	        break;	
	    case F3_ENERGY:
		Z(i) = log(curFaceEnergies(i,2));
	        break;	
	    case TOT_ENERGY:
		Z(i) = log(curFaceEnergies.row(i).norm());
	        break;	
	}	
    }
    switch (shading_state) 
    {
        case NONE: 
            faceColors.setConstant(0.7);
	    break;
	default:
	    igl::colormap(viz_color,Z, true, faceColors);
            break;
    }
    
    viewer.data.set_colors(faceColors);
}

void WeaveHook::renderRenderGeometry(igl::viewer::Viewer &viewer)
{
    viewer.data.clear();
    viewer.data.set_mesh(renderQ, renderF);
    int edges = edgeSegs.rows();
    Eigen::MatrixXd renderPts(2 * edges, 3);
    for (int i = 0; i < edges; i++)
    {
        Eigen::Vector3d vec = edgeVecs.row(i);
        if (normalizeVectors)
        {
            if(vec.norm() != 0.0)
                vec *= baseLength / vec.norm() * sqrt(3.0) / 6.0;
        }
        renderPts.row(2 * i) = edgePts.row(i);
        renderPts.row(2 * i + 1) = edgePts.row(i) + vectorScale*vec.transpose();
    }
    viewer.data.set_edges(renderPts, edgeSegs, edgeColors);      
    setFaceColors(viewer);   
}

bool WeaveHook::simulateOneStep()
{
    //GNtestFiniteDifferences(*weave, params);
    //exit(-1);

    oneStep(*weave, params);
    faceEnergies(*weave, params, tempFaceEnergies);
    return false;
}
