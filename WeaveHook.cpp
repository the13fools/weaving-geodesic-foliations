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
            faceColors = clicked; // faceColors.setConstant(0.7);
            showCutVertexSelection(viewer);     
	    break;
	default:
	    igl::colormap(viz_color,Z, true, faceColors);
            break;
    }
    
    viewer.data.set_colors(faceColors);
}


void WeaveHook::showCutVertexSelection(igl::viewer::Viewer &viewer)
{

    Eigen::MatrixXd selectedVert = Eigen::MatrixXd::Zero(3, 3);
    int counter = 0;
    for (int i = 0; i < vertexSelect.rows(); i++)
    {
        if ( vertexSelect(i) > -1 )
	{
	    if (counter < 2) 
	    { 
		counter++;
		std::cout << i << "\n";
		selectedVert.row(counter) = weave->V.row( weave->F( i, vertexSelect(i) ) );
	    }
    	}

    }	
//    std::cout << counter << "\n"; 

    Eigen::RowVector3d teal(.1, .9, .9);
    viewer.data.add_points( selectedVert, teal ); 
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
	for (int i = 0; i < 3; i++)
	{
        //    trace->traceCurve(*weave, trace_state, i, 1, traceFaceId, traceSteps);
        //    trace->traceCurve(*weave, trace_state, i, -1, traceFaceId, traceSteps);
	}

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

void WeaveHook::updateSingularVerts(igl::viewer::Viewer &viewer)
{

    nonIdentityEdges = Eigen::MatrixXd::Zero(weave->E.size(), 3);
    for (int i = 0; i < weave->Ps.size(); i++)
    {
        bool id = true;
	for (int j = 0; j < 3; j++)
	{
            if (weave->Ps[i](j,j) != 1)
	    {
                id = false;
	    }
	}
	if (!id)
	{
            nonIdentityEdges.row(i) = ( weave->V.row(weave->edgeVerts(i, 0)) + 
		                        weave->V.row(weave->edgeVerts(i, 1)) ) * .5;
	    
	}
    }	

    Eigen::RowVector3d green(.1, .9, .1);
    Eigen::RowVector3d blue(.1, .1, .9);
    viewer.data.add_points( singularVerts_topo, green ); 
    viewer.data.add_points( nonIdentityEdges, blue ); 
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
    updateSingularVerts(viewer);
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
    std::vector<int> topsingularities;
    std::vector<std::pair<int, int> > geosingularities;
    findSingularVertices(*weave, topsingularities, geosingularities);
    std::cout << "now " << topsingularities.size() << " topological and " << geosingularities.size() << " geometric singularities" << std::endl;


    singularVerts_geo = Eigen::MatrixXd::Zero(geosingularities.size(), 3);
    for (int i = 0; i < geosingularities.size(); i++)
    {
//        singularVerts_geo.row(i) = weave->V.row(geosingularities[i]);
    }
    singularVerts_topo = Eigen::MatrixXd::Zero(topsingularities.size(), 3);
    for (int i = 0; i < topsingularities.size(); i++)
    {
        singularVerts_topo.row(i) = weave->V.row(topsingularities[i]);
    }

    nonIdentityEdges = Eigen::MatrixXd::Zero(weave->E.size(), 3);
    for (int i = 0; i < weave->Ps.size(); i++)
    {
        bool id = true;
	for (int j = 0; j < 3; j++)
	{
            if (weave->Ps[i](j,j) != 1)
	    {
                id = false;
	    }
	}
	if (!id)
	{
            nonIdentityEdges.row(i) = ( weave->V.row(weave->edgeVerts(i, 0)) + 
		                        weave->V.row(weave->edgeVerts(i, 1)) ) * .5;
	    
	}
    }
    params.edgeWeights = Eigen::VectorXd::Constant(weave->nEdges(), 1);
}

void WeaveHook::normalizeFields()
{
    weave->normalizeFields();
    updateRenderGeometry();
}

void WeaveHook::serializeVectorField()
{
    weave->serialize(vectorFieldName);
}

void WeaveHook::exportVectorField()
{
    weave->serialize_forexport(vectorFieldName);
}

void WeaveHook::deserializeVectorField()
{
    weave->deserialize(vectorFieldName);
    updateRenderGeometry();
}
