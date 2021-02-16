#include "edge_length_constraint.h"

namespace pbd {

edge_length_constraint_t::scalar_type
edge_length_constraint_t::evaluate(vertices_type const& V, masses_type const& M) const
{
    auto const v0 = indices().at(0);
    auto const v1 = indices().at(1);
    auto const p0 = V.row(v0);
    auto const p1 = V.row(v1);

    return (p1 - p0).norm() - d_;
}

void edge_length_constraint_t::project(vertices_type& V, masses_type const& M) const
{
    auto const& indices = this->indices();
    auto const v0 = indices.at(0);
    auto const v1 = indices.at(1);
    auto const p0 = V.row(v0).transpose();
    auto const p1 = V.row(v1).transpose();
    auto const w0 = 1. / M(v0);
    auto const w1 = 1. / M(v1);
    auto const n  = (p1 - p0).normalized();
    auto const C  = evaluate(V, M);
    
    V.row(v0) = w0 / (w0 + w1) * C * -n.transpose();
    V.row(v1) = w1 / (w0 + w1) * C *  n.transpose();
}

} // namespace pbd
