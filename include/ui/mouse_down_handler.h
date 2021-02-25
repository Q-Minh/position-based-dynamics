#ifndef PBD_MOUSE_DOWN_HANDLER_H
#define PBD_MOUSE_DOWN_HANDLER_H

#include "pbd/deformable_mesh.h"
#include "picking_state.h"

#include <GLFW/glfw3.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/unproject_onto_mesh.h>

namespace ui {

struct mouse_down_handler_t
{
    std::function<bool()> is_model_ready;
    picking_state_t* picking_state;
    pbd::deformable_mesh_t* model;

    mouse_down_handler_t(
        std::function<bool()> is_model_ready,
        picking_state_t* picking_state,
        pbd::deformable_mesh_t* model)
        : is_model_ready(is_model_ready), picking_state(picking_state), model(model)
    {
    }

    bool operator()(igl::opengl::glfw::Viewer& viewer, int button, int modifier);
};

} // namespace ui

#endif // PBD_MOUSE_DOWN_HANDLER_H