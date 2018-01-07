#include "WeaveHook.h"
#include "GaussNewton.h"
#include <iostream>

using namespace std;

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
    viewer.data.set_colors(faceColors);
}

bool WeaveHook::simulateOneStep()
{
    //GNtestFiniteDifferences(*weave, params);
    //exit(-1);

    oneStep(*weave, params);
    return false;
}