#include "solver/solve.h"

namespace solver {

void solve(
    xpbd::deformable_mesh_t& model,
    Eigen::MatrixX3d const& fext,
    double timestep,
    std::uint32_t iterations,
    std::uint32_t substeps)
{
    auto const num_iterations = iterations / substeps;
    double dt                 = timestep / static_cast<double>(substeps);
    auto const& constraints   = model.constraints();
    auto const J              = constraints.size();
    std::vector<double> lagrange_multipliers(J, 0.);

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
        std::fill(lagrange_multipliers.begin(), lagrange_multipliers.end(), 0.0);
        for (auto n = 0u; n < num_iterations; ++n)
        {
            for (auto j = 0u; j < J; ++j)
            {
                auto const& constraint = constraints[j];
                constraint->project(p, m, lagrange_multipliers[j], dt);
            }
        }

        // update solution
        for (auto i = 0u; i < x.rows(); ++i)
        {
            v.row(i) = (p.row(i) - x.row(i)) / dt;
            x.row(i) = p.row(i);
        }

        // friction or other non-conservative forces here ...
    }
}

} // namespace solver
