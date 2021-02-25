#include "ui/mouse_move_handler.h"

namespace ui {

bool mouse_move_handler_t::operator()(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
{
    if (!is_model_ready())
        return false;

    if (!picking_state->is_picking)
        return false;

    double const x1 = static_cast<double>(picking_state->mouse_x);
    double const y1 = viewer.core().viewport(3) - static_cast<double>(picking_state->mouse_y);

    double const x2 = static_cast<double>(viewer.current_mouse_x);
    double const y2 = viewer.core().viewport(3) - static_cast<double>(viewer.current_mouse_y);

    Eigen::Vector3d const p1 = igl::unproject(
                                   Eigen::Vector3f(x1, y1, .5f),
                                   viewer.core().view,
                                   viewer.core().proj,
                                   viewer.core().viewport)
                                   .cast<double>();

    Eigen::Vector3d const p2 = igl::unproject(
                                   Eigen::Vector3f(x2, y2, .5f),
                                   viewer.core().view,
                                   viewer.core().proj,
                                   viewer.core().viewport)
                                   .cast<double>();

    Eigen::Vector3d const direction = (p2 - p1).normalized();

    fext->row(picking_state->vertex) =
        direction.transpose() * static_cast<double>(picking_state->force);

    viewer.data().add_points(
        model->positions().row(picking_state->vertex),
        Eigen::RowVector3d(1., 0., 0.));

    picking_state->mouse_x = viewer.current_mouse_x;
    picking_state->mouse_y = viewer.current_mouse_y;

    return true;
}

} // namespace ui
