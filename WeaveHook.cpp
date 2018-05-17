#include "WeaveHook.h"
#include "GaussNewton.h"
#include <iostream>
#include "Permutations.h"
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include "ImGuiDouble.h"
#include "Surface.h"
#include "CoverMesh.h"

using namespace std;

void WeaveHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
{
    ImGui::InputText("Mesh", meshName);
    ImGui::Combo("GUI Mode", (int *)&gui_mode, cover ? "Weave\0Cover\0\0" : "Weave\0\0");

    if (gui_mode == GUIMode_Enum::WEAVE)
    {
        if (ImGui::CollapsingHeader("Visualization (Weave)", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::InputDoubleScientific("Vector Scale", &vectorScale);
            ImGui::Checkbox("Normalize Vectors", &normalizeVectors);
            ImGui::Checkbox("Hide Vectors", &hideVectors);
            ImGui::Combo("Shading", (int *)&weave_shading_state, "None\0F1 Energy\0F2 Energy\0F3 Energy\0Total Energy\0Connection\0\0");
            if (ImGui::Button("Normalize Fields", ImVec2(-1, 0)))
                normalizeFields();
            ImGui::Checkbox("Fix Fields", &weave->fixFields);
        }
        if (ImGui::CollapsingHeader("Solver Parameters", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::InputDoubleScientific("Compatilibity Lambda", &params.lambdacompat);
            ImGui::InputDoubleScientific("Tikhonov Reg", &params.lambdareg);
            ImGui::InputDoubleScientific("V curl reg", &params.curlreg);
            if (ImGui::Button("Reassign Permutations", ImVec2(-1, 0)))
                reassignPermutations();
            if (ImGui::Button("Remove Singularities", ImVec2(-1, 0)))
                removeSingularities();
            if (ImGui::Button("Create Cover", ImVec2(-1, 0)))
                augmentField();            
        }
        if (ImGui::CollapsingHeader("Save/Load Field", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::InputText("Filename", vectorFieldName);
            if (ImGui::Button("Save Field", ImVec2(-1, 0)))
                serializeVectorField();
            if (ImGui::Button("Load Field", ImVec2(-1, 0)))
                deserializeVectorField();
            if (ImGui::Button("Load Field (Old Format)", ImVec2(-1, 0)))
                deserializeVectorFieldOld();
        }
        if (ImGui::CollapsingHeader("Cuts", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::Button("Reset Cut Select", ImVec2(-1, 0)))
                resetCutSelection();
            if (ImGui::Button("Add Cut", ImVec2(-1, 0)))
                addCut();
            if (ImGui::Button("Remove Prev Cut", ImVec2(-1, 0)))
                removePrevCut();
        }

        if (ImGui::CollapsingHeader("Traces", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::Button("Draw Trace", ImVec2(-1, 0)))
            {
                isDrawTrace = true;
            }
            if (ImGui::Button("Delete Last Trace", ImVec2(-1, 0)))
            {
                isDeleteLastTrace = true;
            }
            if (ImGui::Button("Save Traces", ImVec2(-1, 0)))
            {
                isSaveTrace = true;
            }
        }

        menu.callback_draw_custom_window = [&]()
        {
            // Define next window position + size
            ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(200, 200), ImGuiSetCond_FirstUseEver);
            ImGui::Begin(
                "Handles", nullptr,
                ImGuiWindowFlags_NoSavedSettings
            );

            ImGui::InputInt("Face Location", &handleLocation, 0, 0);
            ImGui::InputDoubleScientific("P0", &handleParams[0]);
            ImGui::InputDoubleScientific("P1", &handleParams[1]);
            ImGui::InputDoubleScientific("P2", &handleParams[2]);
            ImGui::InputDoubleScientific("P3", &handleParams[3]);
            ImGui::InputDoubleScientific("P4", &handleParams[4]);
            ImGui::InputDoubleScientific("P5", &handleParams[5]);

            ImGui::End();

            ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 210), ImGuiSetCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(200, 300), ImGuiSetCond_FirstUseEver);

            ImGui::Begin(
                "Manipulate", nullptr,
                ImGuiWindowFlags_NoSavedSettings
            );

            if (ImGui::CollapsingHeader("Tracing Controls", ImGuiTreeNodeFlags_DefaultOpen))
            {
                ImGui::InputInt("Trace Face", &traceFaceId, 0, 0);
                ImGui::InputInt("Trace Steps", &traceSteps, 0, 0);
                ImGui::InputInt("Trace Field", &traceIdx, 0, 0);
                ImGui::InputInt("Trace Sign", &traceSign, 0, 0);
                ImGui::Combo("Trace Mode", (int *)&trace_state, "Geodesic\0Field\0\0");
                ImGui::Checkbox("Show Bending", &showBending);
                ImGui::InputText("Trace File", traceFile);
                if (ImGui::Button("Save Traces", ImVec2(-1, 0)))
                    saveTraces();
                if (ImGui::Button("Load Traces", ImVec2(-1, 0)))
                    loadTraces();
                if (ImGui::Button("Load Sampled Traces", ImVec2(-1, 0)))
                    loadSampledTraces();
            }
            ImGui::End();
        };
    }
    else if (gui_mode == GUIMode_Enum::COVER)
    {
        if (ImGui::CollapsingHeader("Visualization (Cover)", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::InputDoubleScientific("Vector Scale", &vectorScale);
            ImGui::Checkbox("Hide Vectors", &hideVectors);
            ImGui::Combo("Shading", (int *)&cover_shading_state, "None\0Theta Value\0Connection\0\0");
        }

        if (ImGui::CollapsingHeader("Cover Controls", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::InputDoubleScientific("Initial Scales", &scalesInit);
            if (ImGui::Button("Compute Function Value", ImVec2(-1, 0)))
                computeFunc();
            ImGui::InputInt("Num Isolines", &numISOLines);
            if (ImGui::Button("Draw Isolines", ImVec2(-1, 0)))
                drawISOLines();
        }
    }
}

bool WeaveHook::mouseClicked(igl::opengl::glfw::Viewer &viewer, int button)
{
    if (gui_mode != GUIMode_Enum::WEAVE)
        return false;
    int fid;
    Eigen::Vector3f bc;
    // Cast a ray in the view direction starting from the mouse position
    double x = viewer.current_mouse_x;
    double y = viewer.core.viewport(3) - viewer.current_mouse_y;
    if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.core.model,
        viewer.core.proj, viewer.core.viewport, this->weave->fs->data().V, this->weave->fs->data().F, fid, bc))
    {
        std::cout << fid << " - clicked on vertex #\n"; 
        bool found = false;
        for (int i = 0; i < (int)selectedVertices.size(); i++)
        {
            if(selectedVertices[i].first == fid)
            { 
                found = true;
                if (selectedVertices[i].second < 2)
                    selectedVertices[i].second++;
                else
                    selectedVertices.erase(selectedVertices.begin() + i);
            }
        }
        if(!found && selectedVertices.size() < 2)
        {
            std::pair<int, int> newsel(fid, 0);
            selectedVertices.push_back(newsel);
        }
        renderSelectedVertices.clear();
        for (int i = 0; i < (int)selectedVertices.size(); i++)
        {
            renderSelectedVertices.push_back(weave->fs->data().V.row(weave->fs->data().F(selectedVertices[i].first, selectedVertices[i].second)));
        }
        return true;
    }
    return false;
}

void WeaveHook::initSimulation()
{
    if (weave)
        delete weave;
    if (cover)
        delete cover;
    weave = new Weave(meshName, 3);    
    cover = NULL;
    gui_mode = GUIMode_Enum::WEAVE;
    Handle h;
    h.face = 0;
    h.dir << 1, 0;
    h.field = 2;
    weave->addHandle(h);
    h.face = 0;
    h.dir << 0, 1;
    h.field = 1;
    weave->addHandle(h);
    h.face = 0;
    h.dir << 1, -1;
    h.field = 0;
    weave->addHandle(h);
    curFaceEnergies = Eigen::MatrixXd::Zero(3, 3);
    selectedVertices.clear();
    renderSelectedVertices.clear();
    params.edgeWeights = Eigen::VectorXd::Constant(weave->fs->nEdges(), 1);    

    weave->fixFields = false;    
}

void WeaveHook::setFaceColorsCover(igl::opengl::glfw::Viewer &viewer)
{
    int faces = cover->fs->data().F.rows();
    // if ( curFaceEnergies.rows() != faces && shading_state != NONE) { return ; }
    // cout << "fuck" << endl;

    igl::ColorMapType viz_color = igl::COLOR_MAP_TYPE_MAGMA;

    Eigen::VectorXd Z(faces);
    Z.setConstant(0.7);

    Eigen::VectorXd FVAL(cover->fs->nVerts());

    if (cover_shading_state == FUN_VAL)
    {
        for (int i = 0; i < cover->fs->nVerts(); i++)
        {            
            FVAL(i) = cover->theta[i];            
        }
        viewer.data().set_face_based(false);
    }
    else
    {
        viewer.data().set_face_based(true);
    }

    if (cover_shading_state == CS_CONNECTION_ENERGY)
    {
        cover->fs->connectionEnergy(Z);
    }

    Eigen::MatrixXd faceColors(cover->fs->nFaces(), 3);

    switch (cover_shading_state) 
    {   
    case CS_NONE:
        faceColors.setConstant(0.7);
        break;

    case FUN_VAL:
        igl::colormap(viz_color, FVAL, true, faceColors);        
        if (drawLine)
        {
            for (int i = 0; i < paths.size(); i ++)
            {
                const Eigen::RowVector3d red(0.9,.1,.1), blue(0.1,.1,.9);
                Eigen::MatrixXd line_starts = paths[i].block(0, 0, paths[i].rows() - 1, 3);
                Eigen::MatrixXd line_ends  = paths[i].block(1, 0, paths[i].rows() - 1, 3);
                viewer.data().add_edges( line_starts, line_ends, red);
            }
        }
        break;
    default:
        igl::colormap(viz_color,Z, true, faceColors);
        break;
    }


    viewer.data().set_colors(faceColors);
}

void WeaveHook::setFaceColorsWeave(igl::opengl::glfw::Viewer &viewer)
{
    int faces = weave->fs->data().F.rows();
    // if ( curFaceEnergies.rows() != faces && shading_state != NONE) { return ; }
    // cout << "fuck" << endl;

    igl::ColorMapType viz_color = igl::COLOR_MAP_TYPE_MAGMA;

    Eigen::VectorXd Z(faces);    

    for (int i = 0; i < faces; i++)
    {
        switch (weave_shading_state)
        {
        case F1_ENERGY:
            Z(i) = log(curFaceEnergies(i, 0));
            break;
        case F2_ENERGY:
            Z(i) = log(curFaceEnergies(i, 1));
            break;
        case F3_ENERGY:
            Z(i) = log(curFaceEnergies(i, 2));
            break;
        case TOT_ENERGY:
            Z(i) = log(curFaceEnergies.row(i).norm());
            break;
        case WS_NONE:
        default:
            Z(i) = .7;
            break;
        }
    }

    if (weave_shading_state == WS_CONNECTION_ENERGY)
    {
        weave->fs->connectionEnergy(Z);
    }

    viewer.data().set_face_based(true);
   
    Eigen::MatrixXd faceColors(weave->fs->nFaces(), 3);

    switch (weave_shading_state)
    {
    case WS_NONE:
        faceColors.setConstant(0.7);
        for (int i = 0; i < (int)selectedVertices.size(); i++)
            faceColors.row(selectedVertices[i].first) << 1, 0, 0;
        showCutVertexSelection(viewer);
        break;
    default:
        igl::colormap(viz_color, Z, true, faceColors);
        break;
    }

    viewer.data().set_colors(faceColors);
}


void WeaveHook::showCutVertexSelection(igl::opengl::glfw::Viewer &viewer)
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
    viewer.data().add_points(P, C);
    
}

void WeaveHook::drawCuts(igl::opengl::glfw::Viewer &viewer)
{
    Eigen::RowVector3d blue(0.9, .1, .9);
    Eigen::MatrixXd C(cutPos1.rows(), 3);
    for(int i=0; i<3; i++)
        C.col(i).setConstant(blue[i]);
    viewer.data().add_edges(cutPos1, cutPos2, C);
    
}

void WeaveHook::drawTraceCenterlines(igl::opengl::glfw::Viewer &viewer)
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
            trace->traceCurve(*weave, trace_state, i, 1, traceFaceId, traceSteps);
            trace->traceCurve(*weave, trace_state, i, -1, traceFaceId, traceSteps);
        }

     //   trace->traceCurve(*weave, trace_state, traceIdx, traceSign, traceFaceId, traceSteps);
        isDrawTrace = false;
    }
    if (isSaveTrace)
    {
        trace->logRibbonsToFile( "rods", "example", *weave );
        isSaveTrace = false;    
    }

    Eigen::RowVector3d red(0.9, .1, .1), green(.1, .9, .1);
    for (int i = 0; i < trace->curves.size(); i++)
    {
        int rows = trace->curves[i].rows();
        Eigen::MatrixXd s1 = trace->curves[i].block(0, 0, rows - 1, 3);
        Eigen::MatrixXd s2 = trace->curves[i].block(1, 0, rows - 1, 3);
        Eigen::MatrixXd norms = trace->normals[i].block(0, 0, rows - 1, 3);
         
        switch (trace->modes[i])
        {
            case GEODESIC:
                viewer.data().add_edges(s1, s2, red);
                break;
            case FIELD:
                viewer.data().add_edges(s1, s2, green);
             ///   cout << s1.rows() << " " << trace->normals[i].rows() << "\n";
//z                viewer.data.add_edges(s1, s1 + norms, red);
                if( showBending )
                {
                    Eigen::MatrixXd bend_colors = Eigen::MatrixXd::Zero(s2.rows(),3);
                    igl::ColorMapType viz_color = igl::COLOR_MAP_TYPE_JET;
                    igl::colormap(viz_color,trace->bending[i].tail(s2.rows()), true, bend_colors);
                    viewer.data().add_points( s2, bend_colors );
                }
                break;
        }
    }
}

void WeaveHook::updateSingularVerts(igl::opengl::glfw::Viewer &viewer)
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
    }   */

    Eigen::RowVector3d green(.1, .9, .1);
    Eigen::RowVector3d blue(.1, .1, .9);
    viewer.data().add_points( singularVerts_topo, green ); 
    viewer.data().add_points( nonIdentityEdges, blue ); 
}

void WeaveHook::renderRenderGeometry(igl::opengl::glfw::Viewer &viewer)
{
    viewer.data().clear();

    if (gui_mode == GUIMode_Enum::WEAVE)
    {
        viewer.data().set_mesh(renderQWeave, renderFWeave);
        int edges = edgeSegsWeave.rows();
        Eigen::MatrixXd renderPts(2 * edges, 3);
        for (int i = 0; i < edges; i++)
        {
            Eigen::Vector3d vec = edgeVecsWeave.row(i);
            if (normalizeVectors)
            {
                if (vec.norm() != 0.0)
                    vec *= baseLength / vec.norm() * sqrt(3.0) / 6.0 * 0.75;
            }
            renderPts.row(2 * i) = edgePtsWeave.row(i) - vectorScale*vec.transpose();
            renderPts.row(2 * i + 1) = edgePtsWeave.row(i) + vectorScale*vec.transpose();
        }
        if (!hideVectors)
        {
            viewer.data().set_edges(renderPts, edgeSegsWeave, edgeColorsWeave);
        }
        setFaceColorsWeave(viewer);
        drawTraceCenterlines(viewer);
        updateSingularVerts(viewer);
        drawCuts(viewer);
    }
    else if (gui_mode == GUIMode_Enum::COVER)
    {
        viewer.data().set_mesh(renderQCover, renderFCover);
        int edges = edgeSegsCover.rows();
        Eigen::MatrixXd renderPts(2 * edges, 3);
        for (int i = 0; i < edges; i++)
        {
            Eigen::Vector3d vec = edgeVecsCover.row(i);
            if (normalizeVectors)
            {
                if (vec.norm() != 0.0)
                    vec *= baseLength / vec.norm() * sqrt(3.0) / 6.0 * 0.75;
            }
            renderPts.row(2 * i) = edgePtsCover.row(i) - vectorScale*vec.transpose();
            renderPts.row(2 * i + 1) = edgePtsCover.row(i) + vectorScale*vec.transpose();
        }
        if (!hideVectors)
        {
            viewer.data().set_edges(renderPts, edgeSegsCover, edgeColorsCover);
        }
        setFaceColorsCover(viewer);
    }
}

bool WeaveHook::simulateOneStep()
{
    params.edgeWeights.resize(weave->fs->nEdges());
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
        singularVerts_topo.row(i) = weave->fs->data().V.row(topsingularities[i]);
    }

    nonIdentityEdges = Eigen::MatrixXd::Zero(weave->fs->data().E.size(), 3);
    for (int i = 0; i < weave->fs->Ps_.size(); i++)
    {
        bool id = true;
        for (int j = 0; j < 3; j++)
        {
            if (weave->fs->Ps(i)(j, j) != 1)
            {
                id = false;
            }
        }
        if (!id)
        {
            nonIdentityEdges.row(i) = (weave->fs->data().V.row(weave->fs->data().edgeVerts(i, 0)) +
                weave->fs->data().V.row(weave->fs->data().edgeVerts(i, 1))) * .5;

        }
    }
}

void WeaveHook::normalizeFields()
{
    weave->fs->normalizeFields();
    updateRenderGeometry();
}

void WeaveHook::serializeVectorField()
{
    std::ofstream ofs(vectorFieldName, ios::binary);
    weave->serialize(ofs);
}

void WeaveHook::augmentField()
{
    if (cover) delete cover;
    cover = weave->createCover();    
    updateRenderGeometry();
}

void WeaveHook::computeFunc()
{
    if(cover)
        cover->computeFunc(scalesInit);
}

void WeaveHook::drawISOLines()
{
    paths.clear();
    drawLine = true;
    //TODO
    /*weave->drawISOLines(numISOLines);
    for (int p_cnt = 0; p_cnt < weave->isoLines.size(); p_cnt ++)
    {
        std::vector<Eigen::Vector3d> curPath = weave->isoLines[p_cnt];
        Eigen::MatrixXd path(curPath.size(), 3);
        if (curPath.size() < 10){
            // cout << "skip the " << p_cnt << " iso line" << endl;
            continue;
        }
        for (int i = 0; i < curPath.size(); i ++)
            path.row(i) = curPath[i];
        paths.push_back(path);
        // cout << "Generating the " << p_cnt << " iso line" << endl;
        // Eigen::MatrixXd line_starts = path.block(0, 0, curPath.size() - 1, 3);
        // Eigen::MatrixXd line_ends  = path.block(1, 0, curPath.size() - 1, 3);
        // viewer.data.add_edges( line_starts, line_ends, red);
    }
    trace->loadGeneratedCurves(weave->isoLines, weave->isoNormal);*/
}

void WeaveHook::deserializeVectorField()
{
    std::ifstream ifs(vectorFieldName, ios::binary);
    weave->deserialize(ifs);
    updateRenderGeometry();
}

void WeaveHook::deserializeVectorFieldOld()
{
    std::ifstream ifs(vectorFieldName);
    weave->deserializeOldRelaxFile(ifs);
    updateRenderGeometry();
}

void WeaveHook::resetCutSelection()
{
    selectedVertices.clear();
    renderSelectedVertices.clear();
}

void WeaveHook::removeSingularities()
{
    std::vector<int> topsingularities;
    std::vector<std::pair<int, int> > geosingularities;
    findSingularVertices(*weave, topsingularities, geosingularities);

    std::vector<int> todelete = topsingularities;
//    for (int i = 0; i < geosingularities.size(); i++)
//        todelete.push_back(geosingularities[i].first);

    weave->removePointsFromMesh(todelete);
    updateRenderGeometry();
}

void WeaveHook::addCut()
{
    int idx1 = -1;
    int idx2 = -1;

    std::cout << selectedVertices.size() << std::endl;

    if (selectedVertices.size() < 2)
        return;

    idx1 = weave->fs->data().F(selectedVertices[0].first, selectedVertices[0].second);
    idx2 = weave->fs->data().F(selectedVertices[1].first, selectedVertices[1].second);

    if (idx1 == idx2)
        return;

    Cut c;    
    weave->fs->shortestPath(idx1, idx2, c.path);
    weave->cuts.push_back(c);    
    /*for (int i = 0; i < c.path.size(); i++ )
    {
        std::cout << weave->fs->data().edgeVerts(c.path[i].first, c.path[i].second) << "\n"; // " " <<  weave->edgeVerts.row(c.path[i].first);
    }*/

    updateRenderGeometry();
}

void WeaveHook::removePrevCut()
{
    if (weave->cuts.size() > 0)
        weave->cuts.pop_back();
    updateRenderGeometry();
} 

void WeaveHook::saveTraces()
{
    trace->save(traceFile);
}

void WeaveHook::loadTraces()
{
    trace->load(traceFile);
}

void WeaveHook::loadSampledTraces()
{
    trace->loadSampledCurves(traceFile);
}

void WeaveHook::updateRenderGeometry()
{
    renderQWeave = weave->fs->data().V;
    renderFWeave = weave->fs->data().F;
    weave->createVisualizationEdges(edgePtsWeave, edgeVecsWeave, edgeSegsWeave, edgeColorsWeave);
    weave->createVisualizationCuts(cutPos1, cutPos2);
    baseLength = weave->fs->data().averageEdgeLength;
    curFaceEnergies = tempFaceEnergies;


    weave->handles[0].face = handleLocation;
    weave->handles[0].dir(0) = handleParams(0);
    weave->handles[0].dir(1) = handleParams(1);
    weave->handles[1].face = handleLocation;
    weave->handles[1].dir(0) = handleParams(2);
    weave->handles[1].dir(1) = handleParams(3);
    weave->handles[2].face = handleLocation;
    weave->handles[2].dir(0) = handleParams(4);
    weave->handles[2].dir(1) = handleParams(5);

    if (cover)
    {
        renderQCover = cover->fs->data().V;
        renderFCover = cover->fs->data().F;
        cover->createVisualizationEdges(edgePtsCover, edgeVecsCover, edgeSegsCover, edgeColorsCover);
    }
    else
    {
        renderQCover.resize(0, 3);
        renderFCover.resize(0, 3);
        edgePtsCover.resize(0, 3);
        edgeVecsCover.resize(0, 3);
        edgeSegsCover.resize(0, 2);
        edgeColorsCover.resize(0, 3);
    }
}