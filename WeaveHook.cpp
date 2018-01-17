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
            for (int i = 0; i < (int)selectedVertices.size(); i++)
                faceColors.row(selectedVertices[i].first) << 1, 0, 0;
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
    Eigen::RowVector3d teal(.1, .9, .9);
    int nsel = renderSelectedVertices.size();
    Eigen::MatrixXd P(nsel, 3);
    Eigen::MatrixXd C(nsel, 3);
    for (int i = 0; i < (int)renderSelectedVertices.size(); i++)
    {
        P.row(i) = renderSelectedVertices[i];
        C.row(i) = teal;
    }
    viewer.data.add_points(P, C);
    
}

void WeaveHook::drawCuts(igl::viewer::Viewer &viewer)
{
    Eigen::RowVector3d blue(0.9, .1, .9);
    Eigen::MatrixXd C(cutPos.rows(), 3);
    for(int i=0; i<3; i++)
        C.col(i).setConstant(blue[i]);
    viewer.data.add_points(cutPos, C);
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

/*    nonIdentityEdges = Eigen::MatrixXd::Zero(weave->E.size(), 3);
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
    }	*/

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
    drawCuts(viewer);
}

bool WeaveHook::simulateOneStep()
{
    params.edgeWeights.resize(weave->nEdges());
    params.edgeWeights.setConstant(1.0);
    for (int i = 0; i < (int)weave->cuts.size(); i++)
    {
        for (int j = 0; j < (int)weave->cuts[i].path.size(); j++)
        {
            params.edgeWeights[weave->cuts[i].path[j].first] = 0.0;
        }
    }
    oneStep(*weave, params);
    faceEnergies(*weave, params, tempFaceEnergies);
    return false;
}

void WeaveHook::reassignPermutations()
{
    int flipped = reassignCutPermutations(*weave);
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
            if (weave->Ps[i](j, j) != 1)
            {
                id = false;
            }
        }
        if (!id)
        {
            nonIdentityEdges.row(i) = (weave->V.row(weave->edgeVerts(i, 0)) +
                weave->V.row(weave->edgeVerts(i, 1))) * .5;

        }
    }
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

void WeaveHook::resetCutSelection()
{
    selectedVertices.clear();
    renderSelectedVertices.clear();
}

void WeaveHook::addCut()
{
    //    cuts.push_back();
    //    shortestPath
    int idx1 = -1;
    int idx2 = -1;

    if (selectedVertices.size() < 2)
        return;

    idx1 = weave->F(selectedVertices[0].first, selectedVertices[0].second);
    idx2 = weave->F(selectedVertices[1].first, selectedVertices[1].second);

    if (idx1 == idx2)
        return;

    Cut c;    
    weave->shortestPath(idx1, idx2, c.path);
    weave->cuts.push_back(c);    

    updateRenderGeometry();
}

void WeaveHook::removePrevCut()
{
    if (weave->cuts.size() > 0)
        weave->cuts.pop_back();
    updateRenderGeometry();
} 


