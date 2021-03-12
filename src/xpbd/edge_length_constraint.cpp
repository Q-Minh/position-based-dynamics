#include "xpbd/edge_length_constraint.h"

namespace xpbd {

edge_length_constraint_t::scalar_type
edge_length_constraint_t::evaluate(positions_type const& p, masses_type const& M) const
{
    auto const v0 = indices().at(0);
    auto const v1 = indices().at(1);
    auto const p0 = p.row(v0);
    auto const p1 = p.row(v1);

    return (p0 - p1).norm() - d_;
}

void edge_length_constraint_t::project(
    positions_type& p,
    masses_type const& M,
    scalar_type& lagrange,
    scalar_type const dt) const
{
    auto const& indices = this->indices();
    auto const v0       = indices.at(0);
    auto const v1       = indices.at(1);
    auto const p0       = p.row(v0);
    auto const p1       = p.row(v1);
    auto const w0       = 1. / M(v0);
    auto const w1       = 1. / M(v1);
    auto const n        = (p0 - p1).normalized();
    auto const C        = evaluate(p, M);

    // <n,n> = 1 and <-n,-n> = 1
    scalar_type const weighted_sum_of_gradients = w0 + w1;
    scalar_type const alpha_tilde               = alpha_ / (dt * dt);
    scalar_type const delta_lagrange =
        -(C + alpha_tilde * lagrange) / (weighted_sum_of_gradients + alpha_tilde);

    lagrange += delta_lagrange;
    p.row(v0) += w0 * n * delta_lagrange;
    p.row(v1) += w1 * -n * delta_lagrange;
}

} // namespace xpbd
