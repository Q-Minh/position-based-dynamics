#include "edge_length_constraint.h"

namespace pbd {

edge_length_constraint_t::scalar_type
edge_length_constraint_t::evaluate(positions_type const& p, masses_type const& M) const
{
    auto const v0 = indices().at(0);
    auto const v1 = indices().at(1);
    auto const p0 = p.row(v0);
    auto const p1 = p.row(v1);

    return (p0 - p1).norm() - d_;
}

void edge_length_constraint_t::project(positions_type& p, masses_type const& M) const
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

    p.row(v0) += w0 / (w0 + w1) * C * -n;
    p.row(v1) += w1 / (w0 + w1) * C * n;
}

} // namespace pbd
