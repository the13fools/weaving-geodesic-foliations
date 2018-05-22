#include "WeaveHook.h"
#include "GaussNewton.h"
#include <iostream>
#include "Permutations.h"
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include "ImGuiDouble.h"
#include "Surface.h"
#include "CoverMesh.h"

#include <igl/decimate.h>
#include <igl/upsample.h>

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
            ImGui::Checkbox("Wireframe", &wireframe);
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
            if (ImGui::Button("Reassign Permutations", ImVec2(-1, 0)))
                reassignPermutations();
        }
        if (ImGui::CollapsingHeader("Misc", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::InputInt("Target # faces", &targetResolution, 0, 0);
            if (ImGui::Button("Resample Mesh", ImVec2(-1, 0)))
                resample();
            ImGui::InputInt("Num Isolines", &numISOLines);
            if (ImGui::Button("Draw Isolines", ImVec2(-1, 0)))
                drawISOLines();
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
            ImGui::SetNextWindowSize(ImVec2(200, 500), ImGuiSetCond_FirstUseEver);

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
            }
            if (ImGui::CollapsingHeader("Traces", ImGuiTreeNodeFlags_DefaultOpen))
            {
                ImGui::Checkbox("Show Traces", &showTraces);

                if (ImGui::Button("Draw Trace", ImVec2(-1, 0)))
                {
                    computeTrace();
                }
                if (ImGui::Button("Delete Last Trace", ImVec2(-1, 0)))
                {
                    deleteLastTrace();
                }
                if (ImGui::Button("Clear All Traces", ImVec2(-1, 0)))
                {
                    clearTraces();
                }                
            }
            if (ImGui::CollapsingHeader("Rods", ImGuiTreeNodeFlags_DefaultOpen))
            {
                ImGui::Checkbox("Show Rod Segments", &showRatTraces);
                ImGui::InputDoubleScientific("Extend Traces By", &extendTrace);
                ImGui::InputDoubleScientific("Segment Lenght", &segLen);
                ImGui::InputDoubleScientific("Max Curvature", &maxCurvature);
                ImGui::InputDoubleScientific("Min Rod Length", &minRodLen);
                if (ImGui::Button("Generate From Traces", ImVec2(-1, 0)))
                {
                    rationalizeTraces();
                }
                if (ImGui::Button("Save To Rod File", ImVec2(-1, 0)))
                {
//                    saveRods();
                }
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
            ImGui::Combo("Shading", (int *)&cover_shading_state, "None\0S Value\0Theta Value\0Connection\0\0");
            ImGui::Checkbox("Show Cuts", &showCoverCuts);
        }

        if (ImGui::CollapsingHeader("Cover Controls", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::InputDoubleScientific("Regularization", &initSReg);
            if (ImGui::Button("Initialize S", ImVec2(-1,0)))
                initializeS();
            ImGui::InputDoubleScientific("Global Rescaling", &globalSScale);
            if (ImGui::Button("Compute Function Value", ImVec2(-1, 0)))
                computeFunc();
            ImGui::InputInt("Num Isolines", &numISOLines);
            if (ImGui::Button("Draw Isolines", ImVec2(-1, 0)))
                drawISOLines();
        }

        menu.callback_draw_custom_window = NULL;
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

void WeaveHook::clear()
{
    if (cover)
        delete cover;
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
    
    singularVerts_topo.resize(0,3);
    singularVerts_geo.resize(0,3);
    nonIdentityEdges.resize(0,3);
    cutPos1Weave.resize(0,3);
    cutPos2Weave.resize(0,3);

    weave->fixFields = false;    
    pathstarts.resize(0,3);
    pathends.resize(0,3);

    traces.clear();
}

void WeaveHook::initSimulation()
{
    if (weave)
        delete weave;
    weave = new Weave(meshName, 3);    
    clear();    
}

void WeaveHook::resample()
{
    Eigen::MatrixXd Vcurr = weave->fs->data().V;
    Eigen::MatrixXi Fcurr  = weave->fs->data().F;
    while ( Fcurr.rows() < targetResolution * 2)
    {
        Eigen::MatrixXd Vtmp = Vcurr;
        Eigen::MatrixXi Ftmp = Fcurr;
        igl::upsample(Vtmp, Ftmp, Vcurr, Fcurr);
    }
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::VectorXi J;
    
    igl::decimate(Vcurr, Fcurr, targetResolution, V, F, J);
    igl::writeOBJ("resampled.obj",V,F);
    
    delete weave;
    
    weave = new Weave(V, F, 3);    
    clear();  

    // Hacky... 
    updateRenderGeometry();
}

void WeaveHook::setFaceColorsCover(igl::opengl::glfw::Viewer &viewer)
{
    int faces = cover->fs->data().F.rows();
    // if ( curFaceEnergies.rows() != faces && shading_state != NONE) { return ; }
    // cout << "fuck" << endl;

    igl::ColorMapType viz_color = igl::COLOR_MAP_TYPE_MAGMA;

    Eigen::VectorXd Z(faces);
    Z.setConstant(0.7);

    int nsplitverts = cover->splitMesh().nVerts();
    Eigen::VectorXd FVAL(nsplitverts);

    if (cover_shading_state == FUN_VAL)
    {
        for (int i = 0; i < nsplitverts; i++)
        {            
            FVAL(i) = cover->theta[cover->visMeshToCoverMesh(i)];            
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
    
    if (cover_shading_state == CS_S_VAL)
    {
        Z = cover->s;
    }

    Eigen::MatrixXd faceColors(cover->fs->nFaces(), 3);

    switch (cover_shading_state) 
    {   
    case CS_NONE:
        faceColors.setConstant(0.7);
        break;

    case FUN_VAL:
        igl::colormap(viz_color, FVAL, true, faceColors);        
        break;
    default:
        igl::colormap(viz_color,Z, true, faceColors);
        break;
    }
    
    const Eigen::RowVector3d green(.1,.9,.1);
    viewer.data().add_edges( pathstarts, pathends, green);
    
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
    if (gui_mode == GUIMode_Enum::WEAVE)
    {
        Eigen::RowVector3d blue(0.9, .1, .9);
        Eigen::MatrixXd C(cutPos1Weave.rows(), 3);
        for (int i = 0; i < 3; i++)
            C.col(i).setConstant(blue[i]);
        viewer.data().add_edges(cutPos1Weave, cutPos2Weave, C);
    }
    else if (gui_mode == GUIMode_Enum::COVER)
    {
        Eigen::RowVector3d blue(0.9, .1, .9);
        Eigen::MatrixXd C(cutPos1Cover.rows(), 3);
        for (int i = 0; i < 3; i++)
            C.col(i).setConstant(blue[i]);
        viewer.data().add_edges(cutPos1Cover, cutPos2Cover, cutColorsCover);
    }
}

void WeaveHook::deleteLastTrace()
{
    traces.popLastCurve();
    updateRenderGeometry();
}

void WeaveHook::clearTraces()
{
    traces.clear();
    updateRenderGeometry();
}

void WeaveHook::rationalizeTraces()
{
    traces.rationalizeTraces(maxCurvature, extendTrace, segLen, minRodLen);
    updateRenderGeometry();
}

void WeaveHook::computeTrace()
{
    for (int i = 0; i < 3; i++)
    {
        traces.traceCurve(*weave->fs, trace_state, i, 1, traceFaceId, traceSteps);
        traces.traceCurve(*weave->fs, trace_state, i, -1, traceFaceId, traceSteps);
    }
    updateRenderGeometry();
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
        if(showTraces)
            viewer.data().add_edges( tracestarts, traceends, tracecolors );
        Eigen::RowVector3d orange(0.9, 0.9, 0.1);
        Eigen::RowVector3d red(0.9, 0.1, 0.1);
        if (showRatTraces)
        {
            viewer.data().add_edges(rattracestarts, rattraceends, orange);
            viewer.data().add_points(ratcollisions, red);
        }

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
            if (vec.norm() != 0.0)
                vec *= cover->renderScale() * baseLength / vec.norm() * sqrt(3.0) / 6.0 * 0.75;
            renderPts.row(2 * i) = edgePtsCover.row(i) - vectorScale*vec.transpose();
            renderPts.row(2 * i + 1) = edgePtsCover.row(i) + vectorScale*vec.transpose();
        }
        if (!hideVectors)
        {
            viewer.data().set_edges(renderPts, edgeSegsCover, edgeColorsCover);
        }
        setFaceColorsCover(viewer);
        if(showCoverCuts)
            drawCuts(viewer);
    }

    viewer.data().show_faces = !wireframe;
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

    std::vector<Eigen::Vector3d> centers;
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
            Eigen::Vector3d midpt = (weave->fs->data().V.row(weave->fs->data().edgeVerts(i, 0)) + weave->fs->data().V.row(weave->fs->data().edgeVerts(i, 1))) * .5;
            centers.push_back(midpt);

        }
    }
    int ncenters = centers.size();
    nonIdentityEdges.resize(ncenters, 3);
    for(int i=0; i<ncenters; i++)
        nonIdentityEdges.row(i) = centers[i];
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
    if (cover)
    {
        traces.purgeTraces(cover->fs);
        delete cover;
    }
    cover = weave->createCover();    
    updateRenderGeometry();
    gui_mode = GUIMode_Enum::COVER;
}

void WeaveHook::computeFunc()
{
    if(cover)
        cover->computeFunc(globalSScale);
    updateRenderGeometry();
}

void WeaveHook::drawISOLines()
{
    if(cover)
    {
        traces.purgeTraces(cover->fs);
        std::vector<Trace> newtraces;
        cover->recomputeIsolines(numISOLines, newtraces);
        for (auto it : newtraces)
            traces.addTrace(it);
    }

    updateRenderGeometry();
}

void WeaveHook::deserializeVectorField()
{
    clear();
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
    for (int i = 0; i < geosingularities.size(); i++)
        todelete.push_back(geosingularities[i].first);

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

void WeaveHook::updateRenderGeometry()
{
    renderQWeave = weave->fs->data().V;
    renderFWeave = weave->fs->data().F;
    weave->createVisualizationEdges(edgePtsWeave, edgeVecsWeave, edgeSegsWeave, edgeColorsWeave);
    weave->createVisualizationCuts(cutPos1Weave, cutPos2Weave);
    baseLength = weave->fs->data().averageEdgeLength;
    curFaceEnergies = tempFaceEnergies;

    if (weave->handles.size() < 3)
        weave->handles.resize(3);
    weave->handles[0].face = handleLocation;
    weave->handles[0].dir(0) = handleParams(0);
    weave->handles[0].dir(1) = handleParams(1);
    weave->handles[1].face = handleLocation;
    weave->handles[1].dir(0) = handleParams(2);
    weave->handles[1].dir(1) = handleParams(3);
    weave->handles[2].face = handleLocation;
    weave->handles[2].dir(0) = handleParams(4);
    weave->handles[2].dir(1) = handleParams(5);

    int tracesegs = 0;
    for (int i = 0; i < traces.nTraces(); i++)
    {
        tracesegs += traces.trace(i).segs.size();
    }
    tracestarts.resize(tracesegs, 3);
    traceends.resize(tracesegs, 3);
    tracecolors.resize(tracesegs, 3);
    Eigen::RowVector3d red(0.9, .1, .1), green(.1, .9, .1);

    int tridx = 0;
    for (int i = 0; i < traces.nTraces(); i++)
    {
        int nsegs = traces.trace(i).segs.size();
        std::vector<Eigen::Vector3d> verts;
        std::vector<Eigen::Vector3d> normals;
        traces.renderTrace(i, verts, normals);
        for (int j = 0; j < nsegs; j++)
        {
            tracestarts.row(tridx) = verts[j].transpose() + 0.0001*normals[j].transpose();
            traceends.row(tridx) = verts[j + 1].transpose() + 0.0001*normals[j].transpose();

            switch (traces.trace(i).type_)
            {
            case GEODESIC:
                tracecolors.row(tridx) = red;
                break;
            case FIELD:
                tracecolors.row(tridx) = green;
                break;
            }

            tridx++;
        }
    }

    int rattracesegs = 0;
    for (int i = 0; i < traces.nRationalizedTraces(); i++)
    {
        rattracesegs += traces.rationalizedTrace(i).pts.rows() - 1;
    }
    rattracestarts.resize(rattracesegs, 3);
    rattraceends.resize(rattracesegs, 3);

    tridx = 0;
    for (int i = 0; i < traces.nRationalizedTraces(); i++)
    {
        int nsegs = traces.rationalizedTrace(i).pts.rows() - 1;
        for (int j = 0; j < nsegs; j++)
        {
            Eigen::Vector3d normal = traces.rationalizedTrace(i).normals.row(j);
            rattracestarts.row(tridx) = traces.rationalizedTrace(i).pts.row(j) + 0.001*normal.transpose();
            rattraceends.row(tridx) = traces.rationalizedTrace(i).pts.row(j+1) + 0.001*normal.transpose();
            tridx++;
        }
    }

    ratcollisions.resize(2*traces.nCollisions(), 3);
    for (int i = 0; i < traces.nCollisions(); i++)
    {
        Eigen::Vector3d pt0, pt1;
        traces.collisionPoint(i, pt0, pt1);
        ratcollisions.row(2 * i) = pt0.transpose();
        ratcollisions.row(2 * i + 1) = pt1.transpose();
    }

    if (cover)
    {
        cover->createVisualization(renderQCover, renderFCover, edgePtsCover, edgeVecsCover, edgeSegsCover, edgeColorsCover,
            cutPos1Cover, cutPos2Cover, cutColorsCover);

        int totsegs = 0;
        int ntraces = traces.nTraces();
        for (int i = 0; i < ntraces; i++)
        {
            if (traces.trace(i).parent_ == cover->fs)
                totsegs += traces.trace(i).segs.size();
        }
        pathstarts.resize(totsegs, 3);
        pathends.resize(totsegs, 3);
        int idx = 0;
        for (int j = 0; j < ntraces; j++)
        {
            if (traces.trace(j).parent_ != cover->fs)
                continue;

            Eigen::MatrixXd pathstart, pathend;
            cover->drawTraceOnSplitMesh(traces.trace(j), pathstart, pathend);
            int nsegs = traces.trace(j).segs.size();
            for (int i = 0; i < nsegs; i++)
            {
                pathstarts.row(idx) = pathstart.row(i);
                pathends.row(idx) = pathend.row(i);
                idx++;
            }
        }
    }
    else
    {
        renderQCover.resize(0, 3);
        renderFCover.resize(0, 3);
        edgePtsCover.resize(0, 3);
        edgeVecsCover.resize(0, 3);
        edgeSegsCover.resize(0, 2);
        edgeColorsCover.resize(0, 3);
        cutPos1Cover.resize(0, 3);
        cutPos2Cover.resize(0, 3);
        cutColorsCover.resize(0, 3);
        pathstarts.resize(0, 3);
        pathends.resize(0, 3);
    }
}

void WeaveHook::initializeS()
{
    if(cover) cover->initializeS(initSReg);
    updateRenderGeometry();
}
