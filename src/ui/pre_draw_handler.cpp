#include "ui/pre_draw_handler.h"

namespace ui {

bool pre_draw_handler_t::operator()(igl::opengl::glfw::Viewer& viewer)
{
    if (!is_model_ready())
        return false;

    if (viewer.core().is_animating)
    {
        fext->col(1).array() -= physics_params->is_gravity_active ? 9.81 : 0.;
        pbd::solve(
            *model,
            *fext,
            physics_params->dt,
            static_cast<std::uint32_t>(physics_params->solver_iterations),
            static_cast<std::uint32_t>(physics_params->solver_substeps));
        fext->setZero();
        viewer.data().clear();
        viewer.data().set_mesh(model->positions(), model->faces());
    }

    for (auto i = 0u; i < model->positions().rows(); ++i)
    {
        if (!model->is_fixed(i))
            continue;

        viewer.data().add_points(model->positions().row(i), Eigen::RowVector3d{1., 0., 0.});
    }

    return false; // do not return from drawing loop
}

} // namespace ui
