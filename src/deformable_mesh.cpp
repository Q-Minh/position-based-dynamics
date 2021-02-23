#include "deformable_mesh.h"

#include "edge_length_constraint.h"

#include <array>
#include <igl/barycenter.h>
#include <igl/boundary_facets.h>
#include <igl/copyleft/tetgen/cdt.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/edges.h>
#include <igl/winding_number.h>
#include <iostream>
#include <unordered_map>

namespace pbd {

void deformable_mesh_t::tetrahedralize(Eigen::MatrixXd const& V, Eigen::MatrixXi const& F)
{
    Eigen::MatrixXd TV;
    Eigen::MatrixXi TT, TF;
    igl::copyleft::tetgen::CDTParam cdt_params;
    cdt_params.flags = "cpYR";
    if (igl::copyleft::tetgen::cdt(V, F, cdt_params, TV, TT, TF))
        return;

    TT = TT.rowwise().reverse().eval();
    TF = TF.rowwise().reverse().eval();

    Eigen::MatrixXd BC;
    igl::barycenter(TV, TT, BC);

    Eigen::VectorXd W;
    igl::winding_number(V, F, BC, W);

    Eigen::MatrixXi IT((W.array() > 0.5).count(), 4);
    std::size_t k = 0u;
    for (auto t = 0; t < TT.rows(); ++t)
    {
        if (W(t) <= 0.5)
            continue;

        IT.row(k++) = TT.row(t);
    }

    Eigen::MatrixXi G;
    igl::boundary_facets(IT, G);
    G = G.rowwise().reverse().eval();

    this->p_ = TV;
    this->E_ = IT;
    this->F_ = G;
    this->constraints_.clear();
}

void deformable_mesh_t::constrain_edge_lengths()
{
    auto const& positions = this->positions();
    auto const& faces     = this->faces();
    auto const& elements  = this->elements();

    using edges_type = std::unordered_multimap<int, std::pair<int, bool>>;
    edges_type edges;
    Eigen::MatrixXi E;
    igl::edges(elements, E);

    for (auto i = 0u; i < E.rows(); ++i)
    {
        auto const edge = E.row(i);
        auto const e0   = edge(0);
        auto const e1   = edge(1);

        auto const d = (positions.row(e0) - positions.row(e1)).norm();

        auto constraint = std::make_unique<edge_length_constraint_t>(
            std::initializer_list<std::uint32_t>{
                static_cast<std::uint32_t>(e0),
                static_cast<std::uint32_t>(e1)},
            d);

        this->constraints().push_back(std::move(constraint));
    }
}

} // namespace pbd