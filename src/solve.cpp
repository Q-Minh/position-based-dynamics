#include "solve.h"

namespace pbd {

void solve(
    deformable_mesh_t& model,
    Eigen::MatrixX3d const& fext,
    double timestep,
    std::uint32_t iterations,
    std::uint32_t substeps)
{
    auto const num_iterations = iterations / substeps;
    double dt                 = timestep / static_cast<double>(substeps);

    for (auto s = 0u; s < substeps; ++s)
    {
        auto& v = model.velocity();
        auto& x = model.positions();

        auto const& m            = model.mass();
        Eigen::MatrixX3d const a = fext.array().colwise() / m.array();

        // explicit euler step
        auto vexplicit    = v + dt * a;
        Eigen::MatrixXd p = x + dt * vexplicit;

        // generate collision constraints here ...

        // sequential gauss seidel type solve
        for (auto n = 0u; n < num_iterations; ++n)
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
}

} // namespace pbd
