#include "ui/mouse_down_handler.h"

namespace ui {

bool mouse_down_handler_t::operator()(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
{
    if (!is_model_ready())
        return false;

    using button_type = igl::opengl::glfw::Viewer::MouseButton;
    if (static_cast<button_type>(button) != button_type::Left)
        return false;

    bool const process_pick = modifier == GLFW_MOD_CONTROL || modifier == GLFW_MOD_SHIFT;

    int fid;
    double const x = static_cast<double>(viewer.current_mouse_x);
    double const y = viewer.core().viewport(3) - static_cast<double>(viewer.current_mouse_y);
    picking_state->mouse_x = viewer.current_mouse_x;
    picking_state->mouse_y = viewer.current_mouse_y;

    Eigen::Vector3f bc{};

    bool const hit = igl::unproject_onto_mesh(
        Eigen::Vector2f(x, y),
        viewer.core().view,
        viewer.core().proj,
        viewer.core().viewport,
        model->positions(),
        model->faces(),
        fid,
        bc);

    if (!hit)
        return false;

    Eigen::Vector3i const face{
        model->faces()(fid, 0),
        model->faces()(fid, 1),
        model->faces()(fid, 2)};
    unsigned int closest_vertex = face(0);

    if (bc(1) > bc(0) && bc(1) > bc(2))
    {
        closest_vertex = face(1);
    }
    else if (bc(2) > bc(0) && bc(2) > bc(1))
    {
        closest_vertex = face(2);
    }

    if (modifier == GLFW_MOD_CONTROL)
    {
        picking_state->is_picking = true;
        picking_state->vertex     = closest_vertex;
    }
    if (modifier == GLFW_MOD_SHIFT)
    {
        model->toggle_fixed(closest_vertex, physics_params->mass_per_particle);
    }

    return process_pick;
}

} // namespace ui
