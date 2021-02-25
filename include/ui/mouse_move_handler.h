#ifndef PBD_MOUSE_MOVE_HANDLER_H
#define PBD_MOUSE_MOVE_HANDLER_H

#include "pbd/deformable_mesh.h"
#include "picking_state.h"

#include <igl/opengl/glfw/Viewer.h>
#include <igl/unproject.h>

namespace ui {

struct mouse_move_handler_t
{
    std::function<bool()> is_model_ready;
    picking_state_t* picking_state;
    pbd::deformable_mesh_t* model;
    Eigen::MatrixX3d* fext;

    mouse_move_handler_t(
        std::function<bool()> is_model_ready,
        picking_state_t* picking_state,
        pbd::deformable_mesh_t* model,
        Eigen::MatrixX3d* fext)
        : is_model_ready(is_model_ready), picking_state(picking_state), model(model), fext(fext)
    {
    }

    bool operator()(igl::opengl::glfw::Viewer& viewer, int button, int modifier);
};

} // namespace ui

#endif // PBD_MOUSE_MOVE_HANDLER_H