#include <igl/viewer/Viewer.h>
#include <thread>
#include "WeaveHook.h"

static PhysicsHook *hook = NULL;

void toggleSimulation()
{
    if (!hook)
        return;

    if (hook->isPaused())
        hook->run();
    else
        hook->pause();
}

void resetSimulation()
{
    if (!hook)
        return;

    hook->reset();
}

bool drawCallback(igl::viewer::Viewer &viewer)
{
    if (!hook)
        return false;

    hook->render(viewer);
    return false;
}

bool keyCallback(igl::viewer::Viewer& viewer, unsigned int key, int modifiers)
{
    if (key == ' ')
    {
        toggleSimulation();
        return true;
    }
    return false;
}

void drawTrace()
{
    hook->isDrawTrace = true;
}

void deleteLastTrace()
{
    hook->isDeleteLastTrace = true;
}

void saveTraces()
{
    hook->isSaveTrace = true;
}

bool initGUI(igl::viewer::Viewer &viewer)
{
    viewer.ngui->window()->setVisible(false);
    viewer.ngui->addWindow(Eigen::Vector2i(10, 10), "Weaving");
    viewer.ngui->addButton("Run/Pause Sim", toggleSimulation);
    viewer.ngui->addButton("Reset Sim", resetSimulation);
    hook->initGUI(viewer);
    viewer.ngui->addButton("Draw Trace", drawTrace);
    viewer.ngui->addButton("Delete Last Trace", deleteLastTrace);
    viewer.ngui->addButton("Save Traces to Rod", saveTraces);
    viewer.screen->performLayout();
    return false;
}

int main(int argc, char *argv[])
{
  igl::viewer::Viewer viewer;

  hook = new WeaveHook();
  hook->reset();

  viewer.data.set_face_based(true);
  viewer.core.is_animating = true;
  viewer.callback_key_pressed = keyCallback;
  viewer.callback_pre_draw = drawCallback;
  viewer.callback_init = initGUI;
  viewer.launch();
}
