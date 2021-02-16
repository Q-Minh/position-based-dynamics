#include "deformable_mesh.h"

#include "edge_length_constraint.h"

#include <unordered_map>

namespace pbd {

void deformable_mesh_t::constrain_edge_lengths()
{
    auto const& positions = this->positions();
    auto const& faces     = this->faces();
    auto const& elements  = this->elements();

    using edges_type = std::unordered_multimap<int, std::pair<int, bool>>;
    edges_type edges;
    edges.reserve(6u * positions.size());
    for (auto i = 0u; i < this->elements().rows(); ++i)
    {
        auto const element = this->elements().row(i);
        for (int j = 0; j < element.cols(); ++j)
        {
            auto const e0 = faces(i, j);
            auto const e1 = faces(i, (j + 1) % 3);

            edges.insert({e0, {e1, true}});
        }
    }

    for (auto& edge : edges)
    {
        auto const& e0         = edge.first;
        auto const& e1         = edge.second.first;
        bool const& is_counted = edge.second.second;

        if (!is_counted)
            continue;

        auto const p0   = positions.row(e0);
        auto const p1   = positions.row(e1);
        auto const d    = (p1 - p0).norm();
        auto constraint = std::make_unique<edge_length_constraint_t>(
            std::initializer_list<std::uint32_t>{
                static_cast<std::uint32_t>(e0),
                static_cast<std::uint32_t>(e1)},
            d);

        this->constraints().push_back(std::move(constraint));

        auto range        = edges.equal_range(e1);
        auto reverse_edge = std::find_if(
            range.first,
            range.second,
            [=](typename edges_type::value_type const& key_value) {
                auto const& reverse_edge_endpoint = key_value.second.first;
                return reverse_edge_endpoint == e0;
            });

        if (reverse_edge == range.second)
            continue;

        bool& is_reverse_edge_counted = reverse_edge->second.second;
        is_reverse_edge_counted       = false;
    }
}

} // namespace pbd