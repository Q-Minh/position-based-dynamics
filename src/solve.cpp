#include "solve.h"

namespace pbd {

void solve(
    deformable_mesh_t& model,
    Eigen::MatrixX3d const& fext,
    double dt,
    std::uint32_t iterations)
{
    auto& v = model.velocity();
    auto& x = model.positions();

    auto const& m            = model.mass();
    Eigen::VectorXd const w  = 1. / m.array();
    Eigen::MatrixX3d const a = fext.array().colwise() * w.array();

    // explicit euler step
    auto vexplicit    = v + dt * a;
    Eigen::MatrixXd p = x + dt * vexplicit;

    // generate collision constraints here ...

    // sequential gauss seidel type solve
    for (auto n = 0u; n < iterations; ++n)
    {
        for (auto const& constraint : model.constraints())
            constraint->project(p, m);
    }

    // update solution
    for (auto i = 0u; i < x.rows(); ++i)
    {
        if (model.is_fixed(i))
            continue;

        v.row(i) = (p.row(i) - x.row(i)) / dt;
        x.row(i) = p.row(i);
    }

    // friction or other non-conservative forces here ...
}

} // namespace pbd
