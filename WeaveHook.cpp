#include "WeaveHook.h"
#include "GaussNewton.h"
#include <iostream>
#include "Permutations.h"

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

void WeaveHook::drawTraceCenterlines(igl::viewer::Viewer &viewer)
{
    if (isDeleteLastTrace)
    {
        trace->popLastCurve();
        isDeleteLastTrace = false;
    }
    if (isDrawTrace)
    {
        trace->traceCurve(*weave, trace_state, traceIdx, traceSign, traceFaceId, traceSteps);
        isDrawTrace = false;
    }
    if (isSaveTrace)
    {
        trace->logRibbonsToFile( "rods", "example" );
        isSaveTrace = false;	
    }

    Eigen::RowVector3d red(0.9, .1, .1), green(.1, .9, .1);
    for (int i = 0; i < trace->curves.size(); i++)
    {
        int rows = trace->curves[i].rows();
        Eigen::MatrixXd s1 = trace->curves[i].block(0, 0, rows - 1, 3);
        Eigen::MatrixXd s2 = trace->curves[i].block(1, 0, rows - 1, 3);
         
	switch (trace->modes[i])
	{
	    case GEODESIC:
                viewer.data.add_edges(s1, s2, red);
		break;
	    case FIELD:
		viewer.data.add_edges(s1, s2, green);
                if( showBending )
		{
                    Eigen::MatrixXd bend_colors = Eigen::MatrixXd::Zero(s2.rows(),3);
	       	    igl::ColorMapType viz_color = igl::COLOR_MAP_TYPE_JET;
		    igl::colormap(viz_color,trace->bending[i].tail(s2.rows()), true, bend_colors);
		    viewer.data.add_points( s2, bend_colors );
		}
		break;
	}
    }
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
                vec *= baseLength / vec.norm() * sqrt(3.0) / 6.0 * 0.75;
        }
        renderPts.row(2 * i) = edgePts.row(i) - vectorScale*vec.transpose();
        renderPts.row(2 * i + 1) = edgePts.row(i) + vectorScale*vec.transpose();
    }
    if ( !hideVectors )
    {
	viewer.data.set_edges(renderPts, edgeSegs, edgeColors);      
    }
    setFaceColors(viewer);  
    drawTraceCenterlines(viewer); 
}

bool WeaveHook::simulateOneStep()
{
    //GNtestFiniteDifferences(*weave, params);
    //exit(-1);
    //reassignPermutations();
    oneStep(*weave, params);
    faceEnergies(*weave, params, tempFaceEnergies);
    return false;
}

void WeaveHook::reassignPermutations()
{
    int flipped = ::reassignPermutations(*weave);
    std::cout << flipped << " permutations changed" << std::endl;
}

void WeaveHook::normalizeFields()
{
    weave->normalizeFields();
    updateRenderGeometry();
}
