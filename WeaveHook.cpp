#include "WeaveHook.h"
#include "GaussNewton.h"
#include "LinearSolver.h"
#include <iostream>
#include "Permutations.h"
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include "Surface.h"
#include "CoverMesh.h"
#include "CurlLocalIntegration.h"
#include "SpectralLocalIntegration.h"
#include "MIGlobalIntegration.h"
#include "GNGlobalIntegration.h"
#include <igl/decimate.h>
#include <igl/upsample.h>
#include <igl/hsv_to_rgb.h>
#include <igl/local_basis.h>
#include <random>
#include <chrono>

using namespace std;

void WeaveHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
{
    ImGui::InputText("Mesh", meshName);
    ImGui::Combo("GUI Mode", (int *)&gui_mode, cover ? "Weave\0Cover\0\0" : "Weave\0\0");
    int nf = (weave ? weave->fs->nFields() : 0);
    ImGui::Text("Number of fields: %d", nf);
    if (rosyN)
    {
        ImGui::Text("Is %d-RoSy", rosyN);
        if(ImGui::Button("Split into Fields", ImVec2(-1,0)))
            splitFromRoSy();
    }
    else
    {
        if (nf == 1)
        {
            ImGui::InputInt("Symmetry Degree:", &desiredRoSyN, 0, 0);
            if (ImGui::Button("Convert to RoSy", ImVec2(-1, 0)))
                convertToRoSy();
        }
        else
        {
            ImGui::Text("Not RoSy");
        }
    }

    ImGui::Combo("Solver Mode", (int *)&solver_mode, "Curl Free\0Dirchlet (Knoppel '13)\0\0");
    ImGui::Checkbox("Soft Handle Constraint", &params.softHandleConstraint);
    ImGui::Checkbox("Disable Curl Constraint", &params.disableCurlConstraint);
    if (ImGui::Button("One Step", ImVec2(-1, 0)))
    {
        simulateOneStep();
        updateRenderGeometry();
    }

    if (gui_mode == GUIMode_Enum::WEAVE)
    {
        if (ImGui::CollapsingHeader("Visualization (Weave)", ImGuiTreeNodeFlags_DefaultOpen))
        {
            bool needsrender = false;
            if (rosyN)
            {
                needsrender |= ImGui::Combo("Show", (int *)&rosyVisMode, "Nothing\0RoSy\0Rep. vector\0\0");
            }
            else
            {
                needsrender |= ImGui::Combo("Show", (int *)&vectorVisMode, "Nothing\0V\0V and delta\0V + delta\0\0");
            }
            needsrender |= ImGui::InputDouble("Vector Scale", &vectorScale);
            needsrender |= ImGui::Checkbox("Normalize Vectors", &normalizeVectors);
            ImGui::Checkbox("Wireframe", &wireframe);
            ImGui::InputDouble("Handle Scale", &params.handleScale);
            ImGui::Combo("Shading", (int *)&weave_shading_state, "None\0F1 Energy\0F2 Energy\0F3 Energy\0Total Energy\0Connection\0\0");
            if (ImGui::Button("Normalize Fields", ImVec2(-1, 0)))
                normalizeFields();
            ImGui::Checkbox("Fix Fields", &weave->fixFields);

            if (needsrender)
                updateRenderGeometry();
        }
        if (ImGui::CollapsingHeader("Solver Parameters", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::InputDouble("Compatilibity Lambda", &params.lambdacompat);
            ImGui::InputDouble("Tikhonov Reg", &params.lambdareg);
            ImGui::InputDouble("Curl Viz Face threshold", &params.curlreg);

            ImGui::InputDouble("vizVectorCurl", &params.vizVectorCurl);
            ImGui::InputDouble("vizCorrectionCurl", &params.vizCorrectionCurl);
            ImGui::Checkbox("Normalize Viz Vecs", &params.vizNormalizeVecs);
            ImGui::Checkbox("Show Curl Sign", &params.vizShowCurlSign);

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
            if (ImGui::Button("Load Field (Paul)", ImVec2(-1,0)))
                deserializePaulField();
        }
        if (ImGui::CollapsingHeader("Cuts", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::Button("Reset Cut Select", ImVec2(-1, 0)))
                resetCutSelection();
            if (ImGui::Button("Add Cut", ImVec2(-1, 0)))
                addCut();
            if (ImGui::Button("Remove Prev Cut", ImVec2(-1, 0)))
                removePrevCut();
            if (ImGui::Button("Clear All Cuts", ImVec2(-1, 0)))
                clearCuts();
            if (ImGui::Button("Reassign Permutations", ImVec2(-1, 0)))
                reassignPermutations();
        }
        if (ImGui::CollapsingHeader("Misc", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::InputInt("Target # faces", &targetResolution, 0, 0);
            ImGui::InputInt("Num Fields", &fieldCount, 0, 0);
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

            ImGui::InputInt("Handle Face", &handleLocation[0], 0, 0);
            ImGui::InputInt("Handle Field", &handleLocation[1], 0, 0);
            ImGui::InputDouble("P0", &handleParams[0]);
            ImGui::InputDouble("P1", &handleParams[1]);
            ImGui::InputDouble("P2", &handleParams[2]);

            if (ImGui::Button("Add Handle", ImVec2(-1, 0)))
                addHandle();
            if (ImGui::Button("Remove Handle", ImVec2(-1, 0)))
                removeHandle();

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
                if (ImGui::Button("Draw Random Traces", ImVec2(-1, 0)))
                {
                    computeRandomTraces(numRandomTraces);
                }
                ImGui::InputInt("Number of Them", &numRandomTraces);
                
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
                ImGui::InputDouble("Extend Traces By", &extendTrace);
                ImGui::InputDouble("Segment Lenght", &segLen);
                ImGui::InputDouble("Max Curvature", &maxCurvature);
                ImGui::InputDouble("Min Rod Length", &minRodLen);
                if (ImGui::Button("Generate From Traces", ImVec2(-1, 0)))
                {
                    rationalizeTraces();
                }
                if (ImGui::Button("Smooth", ImVec2(-1,0)))
                {
                    traces.smoothRationalizedTraces(0.1);
                    updateRenderGeometry();
                }
                ImGui::InputText("Rod Filename", rodFilename);
                if (ImGui::Button("Save To Rod File", ImVec2(-1, 0)))
                {
                    saveRods();
                }
            }
            ImGui::End();
        };
    }
    else if (gui_mode == GUIMode_Enum::COVER)
    {
        if (ImGui::CollapsingHeader("Visualization (Cover)", ImGuiTreeNodeFlags_DefaultOpen))
        {
            bool needsrender = false;
            needsrender |= ImGui::InputDouble("Vector Scale", &vectorScale);
            needsrender |= ImGui::Checkbox("Hide Vectors", &hideCoverVectors);
            ImGui::Combo("Shading", (int *)&cover_shading_state, "None\0S Value\0Theta Value\0Connection\0Theta Grad Error\0");
            ImGui::Checkbox("Show Cuts", &showCoverCuts);

            if (needsrender)
                updateRenderGeometry();
        }

        if (ImGui::CollapsingHeader("Cover Controls", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::Combo("Local Method", (int *)&local_field_integration_method, "Nothing\0Curl Correction\0Our Spectral\0\0");
            if (local_field_integration_method == LFI_SPECTRAL || local_field_integration_method == LFI_CURLCORRECT)
            {
                ImGui::InputDouble("Regularization", &initSReg);
            }
            ImGui::Combo("Global Method", (int *)&global_field_integration_method, "Our Gauss-Newton\0Mixed Integer\0\0");
            ImGui::InputDouble("Global Rescaling", &globalSScale);
            
            if(global_field_integration_method == GFI_GN)
            {
                ImGui::InputInt("Alternations", &globalAlternations);                
                ImGui::InputInt("Power Iterations", &globalPowerIters);
            }
            else if (global_field_integration_method == GFI_MI)
            {
                ImGui::InputDouble("Anisotropy", &bommesAniso);
                ImGui::InputDouble("Regularization", &globalThetaReg);
            }
            if (ImGui::Button("Compute Function Value", ImVec2(-1, 0)))
                computeFunc();
            ImGui::InputInt("Num Isolines", &numISOLines);
            if (ImGui::Button("Round Antipodal Covers", ImVec2(-1, 0)))
                roundCovers();
            if (ImGui::Button("Draw Isolines", ImVec2(-1, 0)))
                drawISOLines();
            ImGui::InputText("Export Prefix", exportPrefix);
            if (ImGui::Button("Export for Rendering", ImVec2(-1,0)))
                exportForRendering();
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
    if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view,
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

// TODO add handle reset functionality...

void WeaveHook::clear()
{
    if (cover)
        delete cover;
    cover = NULL;
    gui_mode = GUIMode_Enum::WEAVE;
    ls.clearHandles();
    
    for (int i = 0; i < fieldCount; i++)
    {
        Handle h;
        h.face = 0;
        Eigen::Vector3d handleVec(sin( ( 2 * 3.1415 * i) / fieldCount), cos( ( 2 * 3.1415 * i) / fieldCount), 0);
        
        Eigen::Matrix<double, 3, 2> B = weave->fs->data().Bs[h.face];
        Eigen::Matrix<double, 2, 3> toBarys = (B.transpose()*B).inverse() * B.transpose();

        h.dir = toBarys * handleVec;
        h.field = i;
        ls.addHandle(h);
    }
    weave->handles = ls.handles;

    curFaceEnergies = Eigen::MatrixXd::Zero(3, 3);
    selectedVertices.clear();
    renderSelectedVertices.clear();
    params.edgeWeights = Eigen::VectorXd::Constant(weave->fs->nEdges(), 1);    
    
    singularVerts_topo.resize(0,3);
    singularVerts_geo.resize(0,3);
    nonIdentity1Weave.resize(0,3);
    nonIdentity2Weave.resize(0,3);
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
    weave = new Weave(meshName, fieldCount);    
    rosyN = 0;
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
    
    weave = new Weave(V, F, fieldCount); 
    std::cout << fieldCount << std::endl;   
    clear();  

    // Hacky... 
    updateRenderGeometry();
}

void WeaveHook::addHandle()
{
    Handle h;
    h.face = handleLocation[0];
    h.face = h.face < weave->fs->nFaces() ? h.face : weave->fs->nFaces() - 1;

    Eigen::Vector3d handleVec(handleParams[0], handleParams[1], handleParams[2]);
    handleVec.normalize();
    
    Eigen::Matrix<double, 3, 2> B = weave->fs->data().Bs[h.face];
    Eigen::Matrix<double, 2, 3> toBarys = (B.transpose()*B).inverse() * B.transpose();

    h.dir = toBarys * handleVec;
    h.field = handleLocation[1];
    h.field = h.field < weave->fs->nFields() ? h.field : weave->fs->nFields() - 1;
    std::cout << "Just added a handle to face: " << h.face << " field: " << h.field << " projected vector " << (weave->fs->data().Bs[h.face] * h.dir).transpose() << std::endl;

    bool toAdd = true;
    for (int i = 0; i < ls.handles.size(); i++)
    {
        if (ls.handles[i].face == h.face && ls.handles[i].field == h.field)
        {
            ls.handles[i].dir = h.dir;
            toAdd = false;
        }
    }
    if (toAdd)
        ls.addHandle(h);

    weave->handles = ls.handles;
    updateRenderGeometry();
}

void WeaveHook::removeHandle()
{
    if (ls.handles.size() > 0)
        ls.handles.pop_back();
    weave->handles = ls.handles;
    std::cout << " There are now " << ls.handles.size() << " handles " << std::endl;
    weave->handles = ls.handles;
    updateRenderGeometry();
}

void WeaveHook::setFaceColorsCover(igl::opengl::glfw::Viewer &viewer)
{
    int faces = cover->fs->data().F.rows();
    
    igl::ColorMapType viz_color = igl::COLOR_MAP_TYPE_MAGMA;

    int nsplitverts = cover->splitMesh().nVerts();

    Eigen::MatrixXd faceColors(cover->fs->nFaces(), 3);
    Eigen::MatrixXd vertColors(nsplitverts, 3);

    switch (cover_shading_state) 
    {   
        case FUN_VAL:
        {
            const double PI = 3.1415926535898;
            for (int i = 0; i < nsplitverts; i++)
            {            
                double theta = cover->theta[cover->visMeshToCoverMesh(i)];            
                // as HSV: [ 360 * (theta+pi) / 2 pi, 1.0, 0.5 ]
                double H = 360.0 * (theta + PI) / (2.0 * PI);
                igl::hsv_to_rgb(H, 1.0, 0.5, vertColors(i,0), vertColors(i,1), vertColors(i,2));
            }   
            break;
        }
        case CS_S_VAL:
        {
            igl::colormap(viz_color, cover->scales, true, faceColors);        
            break;
        }
        case CS_CONNECTION_ENERGY:
        {
            Eigen::VectorXd Z(faces);
            cover->fs->connectionEnergy(Z, 0., params);
            igl::colormap(viz_color, Z, true, faceColors);   
            break;     
        }
        case CS_GRAD_DEVIATION:
        {
            Eigen::VectorXd Z(faces);
            cover->gradThetaDeviation(Z);
            igl::colormap(viz_color, Z, false, faceColors);
            break;
        }
        case CS_NONE:
        default:
        {
            faceColors.setConstant(0.7);
            break;
        }
    }
    
    // fade deleted faces
    for(int i=0; i<cover->fs->nFaces(); i++)
    {
        if(cover->fs->isFaceDeleted(i))
        {
            faceColors(i,0) = 0.5 + 0.5 * faceColors(i,0);
            faceColors(i,1) *= 0.5;
            faceColors(i,2) *= 0.5;
        }
    }
    
    const Eigen::RowVector3d green(.1,.9,.1);
    viewer.data().add_edges( pathstarts, pathends, green);
    
    if(cover_shading_state == FUN_VAL)
    {
        viewer.data().set_colors(vertColors);
        viewer.data().set_face_based(false);
    }
    else
    {
        viewer.data().set_colors(faceColors);
        viewer.data().set_face_based(true);
    }
}

void WeaveHook::setFaceColorsWeave(igl::opengl::glfw::Viewer &viewer)
{
    int faces = weave->fs->data().F.rows();
    // if ( curFaceEnergies.rows() != faces && shading_state != NONE) { return ; }
    // cout << "fuck" << endl;

    igl::ColorMapType viz_color = igl::COLOR_MAP_TYPE_PARULA;

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
        weave->fs->connectionEnergy(Z, params.curlreg, params); // TODO make real var
    }

    viewer.data().set_face_based(true);
   
    Eigen::MatrixXd faceColors(weave->fs->nFaces(), 3);

    switch (weave_shading_state)
    {
    case WS_NONE:
        faceColors.setConstant(0.7);
        for (int i = 0; i < (int)selectedVertices.size(); i++)
            faceColors.row(selectedVertices[i].first) << 0.5, 0.5, 0;
        showCutVertexSelection(viewer);
        break;
    default:
        igl::colormap(viz_color, Z, true, faceColors);
        break;
    }

    // fade deleted faces
    for(int i=0; i<weave->fs->nFaces(); i++)
    {
        if(weave->fs->isFaceDeleted(i))
        {
            faceColors(i,0) = 0.5 + 0.5 * faceColors(i,0);
            faceColors(i,1) *= 0.5;
            faceColors(i,2) *= 0.5;
        }
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
        Eigen::RowVector3d blue(.1, .1, 0.9);
        Eigen::RowVector3d purple(0.9, .1, .9);
        Eigen::MatrixXd C1(cutPos1Weave.rows(), 3);
        for (int i = 0; i < 3; i++)
            C1.col(i).setConstant(blue[i]);
        viewer.data().add_edges(cutPos1Weave, cutPos2Weave, C1);
        Eigen::MatrixXd C2(nonIdentity1Weave.rows(), 3);
        for(int i=0; i<3; i++)
            C2.col(i).setConstant(purple[i]);
        viewer.data().add_edges(nonIdentity1Weave, nonIdentity2Weave, C2);
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

// geodesic thing
void WeaveHook::computeTrace()
{
    // for (int i = 0; i < 3; i++)
    // {
        traces.traceCurve(*weave->fs, trace_state, 0, 1, traceFaceId, traceSteps);
    //     traces.traceCurve(*weave->fs, trace_state, i, -1, traceFaceId, traceSteps);
    // }
    updateRenderGeometry();
}

void WeaveHook::computeRandomTraces(int numtraces)
{
    int nfaces = weave->fs->nFaces();
    int nfields = weave->fs->nFields();
    
    std::default_random_engine generator;
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_int_distribution<int> facedistribution(0,nfaces-1);
    std::uniform_int_distribution<int> fielddistribution(0,nfields-1);
    std::uniform_int_distribution<int> dirdistribution(0, 1);
    
    for (int i = 0; i < numtraces; i++)
    {
        int face = facedistribution(generator);
        int field = fielddistribution(generator);
        int dir = 2*dirdistribution(generator) - 1;
        traces.traceCurve(*weave->fs, trace_state, field, dir, face, traceSteps);
    }
    updateRenderGeometry();
}

void WeaveHook::updateSingularVerts(igl::opengl::glfw::Viewer &viewer)
{
    Eigen::RowVector3d green(.1, .9, .1);
    Eigen::RowVector3d blue(.1, .1, .9);
    viewer.data().add_points( singularVerts_topo, green ); 
}

void WeaveHook::renderRenderGeometry(igl::opengl::glfw::Viewer &viewer)
{
    viewer.data().clear();

    if (gui_mode == GUIMode_Enum::WEAVE)
    {
        viewer.data().set_mesh(renderQWeave, renderFWeave);
        viewer.data().set_edges(edgePtsWeave, edgeSegsWeave, edgeColorsWeave);
        setFaceColorsWeave(viewer);
        if(showTraces)
            viewer.data().add_edges( tracestarts, traceends, tracecolors );
        Eigen::RowVector3d orange(0.9, 0.9, 0.1);
        Eigen::RowVector3d red(0.9, 0.1, 0.1);
        if (showRatTraces)
        {
            viewer.data().add_edges(rattracestarts, rattraceends, orange);
      //      viewer.data().add_points(ratcollisions, red);
        }

        updateSingularVerts(viewer);
        drawCuts(viewer);
    }
    else if (gui_mode == GUIMode_Enum::COVER)
    {
        viewer.data().set_mesh(renderQCover, renderFCover);        
        viewer.data().set_edges(edgePtsCover, edgeSegsCover, edgeColorsCover);
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
    params.rosyN = rosyN; // make this the same...

    int nfaces = weave->fs->data().F.rows();
    int nfields = weave->fs->nFields();

    if ( solver_mode == Solver_Enum::CURLFREE )
    {
        Eigen::VectorXd primal = weave->fs->vectorFields.segment(0, 2*nfaces*nfields);
        Eigen::VectorXd dual = weave->fs->vectorFields.segment(2*nfaces*nfields, 2*nfaces*nfields);

        const int numDesignIters = 10;
        ls.takeSomeSteps(*weave, params, primal, dual, rosyN != 0, numDesignIters);

        weave->fs->vectorFields.segment(0, 2*nfaces*nfields) = primal;
        weave->fs->vectorFields.segment(2*nfaces*nfields, 2*nfaces*nfields) = dual;


 //       std::cout << primal.transpose()<< std::endl << primal.norm() << std::endl<< std::endl;
  //      std::cout << dual.transpose() << std::endl <<dual.norm() << std::endl<< std::endl;
        std::cout << "primal norm " << primal.norm() << " dual norm " <<dual.norm() <<  std::endl;
    }
    else 
    {
        Eigen::VectorXd curField = weave->fs->vectorFields.segment(0, 2*nfaces*nfields);
        weave->fs->vectorFields.setZero();
        weave->fs->vectorFields.segment(0, 2*nfaces*nfields) = curField;
        oneStep(*weave, params);
        faceEnergies(*weave, params, tempFaceEnergies);
    }
    Eigen::VectorXd temp;
    std::cout << "Total Geodesic Energy" << weave->fs->getGeodesicEnergy(temp, params) << std::endl;

    std::cout << "ran a step" << std::endl;
    return false;
}

void WeaveHook::reassignPermutations()
{
   // int flipped = reassignCutPermutations(*weave);
    int flipped = reassignAllPermutations(*weave);

    std::cout << flipped << " permutations changed" << std::endl;
    
    std::vector<std::pair<int, int> > topsingularities;
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
        singularVerts_topo.row(i) = weave->fs->data().V.row(topsingularities[i].first);
    }

    std::vector<int> nonIdentityEdges;
    int m = weave->fs->nFields();
    for (int i = 0; i < weave->fs->Ps_.size(); i++)
    {
        bool id = true;
        
        for (int j = 0; j < m; j++)
        {
            if (weave->fs->Ps(i)(j, j) != 1)
            {
                id = false;
            }
        }
        if (!id)
        {
            nonIdentityEdges.push_back(i);  // TODO: Fix viz bug!
        }
    }
    int ncuts = nonIdentityEdges.size();
    nonIdentity1Weave.resize(ncuts, 3);
    nonIdentity2Weave.resize(ncuts, 3);
    for(int i=0; i<ncuts; i++)
    {
        Eigen::Vector3d normal(0,0,0);
        int f0 = weave->fs->data().E(nonIdentityEdges[i], 0);
        if(f0 != -1)
            normal += weave->fs->faceNormal(f0);
        int f1 = weave->fs->data().E(nonIdentityEdges[i], 1);
        if(f1 != -1)
            normal += weave->fs->faceNormal(f1);
        Eigen::Vector3d offset = 0.001*normal/ normal.norm();
        nonIdentity1Weave.row(i) = weave->fs->data().V.row(weave->fs->data().edgeVerts(nonIdentityEdges[i], 0)) + offset.transpose();
        nonIdentity2Weave.row(i) = weave->fs->data().V.row(weave->fs->data().edgeVerts(nonIdentityEdges[i], 1)) + offset.transpose();
    }
}

void WeaveHook::normalizeFields()
{
    weave->fs->normalizeFields();
    updateRenderGeometry();
}

void WeaveHook::augmentField()
{
    if (cover)
    {
        traces.purgeTraces(cover->fs);
        delete cover;
        cover = NULL;
    }

    weave->fs->undeleteAllFaces();
    
    if (rosyN)
    {
        Weave *splitWeave = weave->splitFromRosy(rosyN);
        reassignAllPermutations(*splitWeave);

        std::vector<std::pair<int, int> > topsingularities;
        std::vector<std::pair<int, int> > geosingularities;
        findSingularVertices(*splitWeave, topsingularities, geosingularities);

        std::vector<std::pair<int, int> > todelete = topsingularities;
        for (int i = 0; i < geosingularities.size(); i++)
            todelete.push_back(geosingularities[i]);

        cover = splitWeave->createCover(todelete);
        delete splitWeave;
    }
    else
    {
        std::vector<std::pair<int, int> > topsingularities;
        std::vector<std::pair<int, int> > geosingularities;
        findSingularVertices(*weave, topsingularities, geosingularities);

        std::vector<std::pair<int, int> > todelete = topsingularities;
        for (int i = 0; i < geosingularities.size(); i++)
            todelete.push_back(geosingularities[i]);

        cover = weave->createCover(todelete);
    }
    updateRenderGeometry();
    gui_mode = GUIMode_Enum::COVER;
}

void WeaveHook::computeFunc()
{
    if (cover)
    {
        LocalFieldIntegration *method;
        if (local_field_integration_method == LFI_TRIVIAL)
            method = new TrivialLocalIntegration();
        else if (local_field_integration_method == LFI_CURLCORRECT)
            method = new CurlLocalIntegration(initSReg);
        else if(local_field_integration_method == LFI_SPECTRAL)
            method = new SpectralLocalIntegration(initSReg);
        else
        {
            assert(!"Unknown local integration method");
            return;
        }
        GlobalFieldIntegration *gmethod;
        if (global_field_integration_method == GFI_GN)
            gmethod = new GNGlobalIntegration(globalAlternations, globalPowerIters);
        else if(global_field_integration_method == GFI_MI)
            gmethod = new MIGlobalIntegration(bommesAniso, globalThetaReg);

        cover->integrateField(method, gmethod, globalSScale);
        delete method;
        delete gmethod;
    }
    updateRenderGeometry();
}

void WeaveHook::roundCovers()
{
    if (cover && numISOLines > 0)
    {
        cover->roundAntipodalCovers(numISOLines);
    }
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

static const int magic = 0x4242;

void WeaveHook::serializeVectorField()
{
    int currentRLXVersion = 2;
    std::ofstream ofs(vectorFieldName, ios::binary);
    ofs.write((char *)&magic, sizeof(int));
    ofs.write((char *)&currentRLXVersion, sizeof(int));
    ofs.write((char *)&rosyN, sizeof(int));
    weave->serialize(ofs);
}

void WeaveHook::deserializeVectorField()
{
    std::ifstream ifs(vectorFieldName, ios::binary);
    if (!ifs)
    {
        std::cerr << "Cannot open vector field file: " << vectorFieldName << std::endl;
        return;
    }
    clear();    
    int header=0;
    ifs.read((char *)&header, sizeof(int));
    if (header != magic)
    {
           // old version
        ifs.clear();
        ifs.seekg(0);
        rosyN = 0;
        weave->deserialize(ifs);
    }
    else
    {
        int saveversion = 0;
        ifs.read((char *)&saveversion, sizeof(int));
        if (saveversion == 1)
        {
            char isrosyc=0;
            ifs.read(&isrosyc, 1);
            rosyN = (isrosyc ? 3 : 0);
        }
        else
        {
            ifs.read((char *)&rosyN, sizeof(int));
        }
        weave->deserialize(ifs);
    }
    updateRenderGeometry();
}

void WeaveHook::deserializeVectorFieldOld()
{
    std::ifstream ifs(vectorFieldName);
    weave->deserializeOldRelaxFile(ifs);
    rosyN = 0;
    updateRenderGeometry();
}

void WeaveHook::deserializePaulField()
{
    std::ifstream ifs(vectorFieldName);
    weave->deserializePaulFile(ifs);
    rosyN = 0;
    updateRenderGeometry();
}

void WeaveHook::resetCutSelection()
{
    selectedVertices.clear();
    renderSelectedVertices.clear();
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
    if (rosyN)
    {
        weave->createVisualizationEdges(edgePtsWeave, edgeSegsWeave, edgeColorsWeave,
            rosyVisMode, normalizeVectors, vectorScale, rosyN);
    }
    else
    {
        weave->createVisualizationEdges(edgePtsWeave, edgeSegsWeave, edgeColorsWeave,
            vectorVisMode, normalizeVectors, vectorScale);
    }
   
    weave->createVisualizationCuts(cutPos1Weave, cutPos2Weave);
    curFaceEnergies = tempFaceEnergies;

    // TODO: refactor, just used for rendering 
    weave->handles = ls.handles;

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
        cover->createVisualization(renderQCover, renderFCover, edgePtsCover, edgeSegsCover, edgeColorsCover,
            cutPos1Cover, cutPos2Cover, cutColorsCover,
            hideCoverVectors, vectorScale);

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
        edgeSegsCover.resize(0, 2);
        edgeColorsCover.resize(0, 3);
        cutPos1Cover.resize(0, 3);
        cutPos2Cover.resize(0, 3);
        cutColorsCover.resize(0, 3);
        pathstarts.resize(0, 3);
        pathends.resize(0, 3);
    }
}

void WeaveHook::saveRods()
{
    traces.exportRodFile(rodFilename.c_str(), weave->fs->nFaces());
}

void WeaveHook::exportForRendering()
{    
    std::string meshName = exportPrefix + std::string("_mesh.obj");
    igl::writeOBJ(meshName.c_str(), weave->fs->data().V, weave->fs->data().F);
    std::string fieldName = exportPrefix + std::string("_field.csv");
    std::ofstream vfs(fieldName.c_str());
    int nfaces = weave->fs->nFaces();
    int nfields = weave->fs->nFields();
    int nverts = weave->fs->nVerts();
    for(int i=0; i<nfaces; i++)
    {
        Eigen::Vector3d centroid(0,0,0);
        for(int j=0; j<3; j++)
            centroid += weave->fs->data().V.row( weave->fs->data().F(i,j) ).transpose();
        centroid /= 3.0;
        for(int j=0; j<nfields; j++)
        {
            Eigen::Vector3d vf = weave->fs->data().Bs[i] * weave->fs->v(i, j);
            if (vf.norm() != 0.0)
                vf *= weave->fs->data().averageEdgeLength / vf.norm() * sqrt(3.0) / 6.0;
            vfs << centroid[0]-vf[0] << ", " << centroid[1]-vf[1] << ", " << centroid[2]-vf[2] << ", " << centroid[0] + vf[0] << ", " << centroid[1] + vf[1] << ", " << centroid[2] + vf[2] << std::endl;                        
        }
    }
    for(int i=0; i<2*nfields; i++)
    {
        /*std::stringstream ss;
        ss << exportPrefix << "_s_" << i << ".csv";
        std::ofstream sfs(ss.str().c_str());
        for(int j=0; j<nfaces; j++)
            sfs << cover->s[i*nfaces + j] << ",\t 0,\t0" << std::endl;*/
            
        std::stringstream ss2;
        ss2 << exportPrefix << "_theta_" << i << ".csv";
        std::ofstream thetafs(ss2.str().c_str());
        for(int j=0; j<nverts; j++)
        {
            thetafs << cover->theta[cover->visMeshToCoverMesh(i*nverts+j)] << ",\t 0,\t0" << std::endl;
        }
    }

    std::string coverMeshName = exportPrefix + std::string("_covermesh.obj");
    igl::writeOBJ(coverMeshName.c_str(), cover->splitMesh().data().V, cover->splitMesh().data().F);
    for(int i=0; i<2*nfields; i++)
    {       
        std::stringstream ssfb;
        ssfb << exportPrefix << "_facebased_" << i << ".csv";
        std::ofstream fbfs(ssfb.str().c_str());
        Eigen::VectorXd connection(nfaces); 
        Eigen::VectorXd deviation(nfaces);
        cover->fs->connectionEnergy(connection, 0., params);// make sure calculation per face not per cover...
        cover->gradThetaDeviation(deviation);
        for(int j=0; j<nfaces; j++)
        {
            int idx = i*nfaces + j;
            fbfs << cover->scales(idx) << ",\t" << connection(idx) << ",\t" << deviation(idx) << std::endl;
        }
    }

    std::stringstream ss3;
    ss3 << exportPrefix << "_geoeng" << ".csv";
    std::ofstream geoengfs(ss3.str().c_str());
    Eigen::VectorXd energy(nfaces);
    weave->fs->connectionEnergy(energy, 0, params);
    for(int i=0; i<nfaces; i++)
    {
        geoengfs << energy(i) << ",\t 0,\t0" << std::endl;
    }

    std::string cutsname = exportPrefix + std::string("_cuts.csv");
    std::ofstream cfs(cutsname.c_str());
    int nsegs = nonIdentity1Weave.rows();
    for(int i=0; i<nsegs; i++)
    {
        cfs << nonIdentity1Weave(i,0) << ", " << nonIdentity1Weave(i,1) << ", " << nonIdentity1Weave(i,2) << ", " << nonIdentity2Weave(i, 0) << ", " << nonIdentity2Weave(i,1) << ", " << nonIdentity2Weave(i,2) << std::endl;
    }
    std::string singname_topo = exportPrefix + std::string("_toposing.csv");
    std::ofstream singfs_topo(singname_topo.c_str());
    int nsing = singularVerts_topo.rows();
    for(int i=0; i<nsing; i++)
    {
        singfs_topo << singularVerts_topo(i,0) << ", " << singularVerts_topo(i,1) << ", " << singularVerts_topo(i,2) << std::endl;
    }
    
    std::string singname_geom = exportPrefix + std::string("_geomsing.csv");
    std::ofstream singfs(singname_geom.c_str());
    nsing = singularVerts_geo.rows();
    for(int i=0; i<nsing; i++)
    {
        singfs << singularVerts_geo(i,0) << ", " << singularVerts_geo(i,1) << ", " << singularVerts_geo(i,2) << std::endl;
    }
    
    std::string tracename = exportPrefix + std::string("_traces.csv");
    traces.exportForRendering(tracename.c_str());
}

void WeaveHook::convertToRoSy()
{
    if (!weave || rosyN > 0 || desiredRoSyN < 1 || weave->fs->nFields() != 1)
        return;

    weave->convertToRoSy(desiredRoSyN);
    ls.clearHandles();
    weave->handles = ls.handles;
    rosyN = desiredRoSyN;
    updateRenderGeometry();
}

void WeaveHook::splitFromRoSy()
{
    if(!weave || rosyN == 0)
        return;
    
    clear();
    Weave *splitWeave = weave->splitFromRosy(rosyN);
    delete weave;
    weave = splitWeave;
    rosyN = 0;
    reassignAllPermutations(*weave);    
    int m = weave->fs->nFields();
    for (int i = 0; i < weave->fs->Ps_.size(); i++)
    {
        bool id = true;
        for (int j = 0; j < m; j++)
        {
            if (weave->fs->Ps(i)(j, j) != 1)
            {
                id = false;
            }
        }
        if (!id)
        {
            Cut c;
            std::pair<int, int> cutedge(i, 1);
            c.path.push_back(cutedge);
            weave->cuts.push_back(c);
        }
    }

    ls.clearHandles();
    for (int i = 0; i < m; i++)
    {
        Handle h;
        
        h.face = 0;
        h.dir = weave->fs->v(0, i);
        h.field = i;

        ls.addHandle(h);
    }
    weave->handles = ls.handles;

    updateRenderGeometry();
}

void WeaveHook::clearCuts()
{
    weave->cuts.clear();
    updateRenderGeometry();
}
