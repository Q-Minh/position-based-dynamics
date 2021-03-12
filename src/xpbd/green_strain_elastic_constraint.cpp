#include "xpbd/green_strain_elastic_constraint.h"

#include <Eigen/Dense>
#include <eigen/SVD>

namespace xpbd {

green_strain_elastic_constraint_t::green_strain_elastic_constraint_t(
    std::initializer_list<index_type> indices,
    positions_type const& p,
    scalar_type young_modulus,
    scalar_type poisson_ratio)
    : base_type(indices), V0_{0.}, DmInv_{}, mu_{}, lambda_{}
{
    assert(indices.size() == 4u);

    auto const v1 = this->indices().at(0);
    auto const v2 = this->indices().at(1);
    auto const v3 = this->indices().at(2);
    auto const v4 = this->indices().at(3);

    auto const p1 = p.row(v1);
    auto const p2 = p.row(v2);
    auto const p3 = p.row(v3);
    auto const p4 = p.row(v4);

    Eigen::Matrix3d Dm;
    Dm.col(0) = (p1 - p4).transpose();
    Dm.col(1) = (p2 - p4).transpose();
    Dm.col(2) = (p3 - p4).transpose();

    V0_     = std::abs((1. / 6.) * Dm.determinant());
    DmInv_  = Dm.inverse();
    mu_     = (young_modulus) / (2. * (1 + poisson_ratio));
    lambda_ = (young_modulus * poisson_ratio) / ((1 + poisson_ratio) * (1 - 2 * poisson_ratio));
}

green_strain_elastic_constraint_t::scalar_type
green_strain_elastic_constraint_t::signed_volume(positions_type const& V) const
{
    Eigen::RowVector3d const p1 = V.row(indices()[0]);
    Eigen::RowVector3d const p2 = V.row(indices()[1]);
    Eigen::RowVector3d const p3 = V.row(indices()[2]);
    Eigen::RowVector3d const p4 = V.row(indices()[3]);

    auto const vol = (1. / 6.) * (p2 - p1).cross(p3 - p1).dot(p4 - p1);
    // Eigen::Matrix3d Ds;
    // Ds.col(0)      = (p1 - p4).transpose();
    // Ds.col(1)      = (p2 - p4).transpose();
    // Ds.col(2)      = (p3 - p4).transpose();
    // auto const vol = (1. / 6.) * Ds.determinant();
    return vol;
}

void green_strain_elastic_constraint_t::project_inverted(positions_type& p, masses_type const& m)
    const
{
}

void green_strain_elastic_constraint_t::project(positions_type& p, masses_type const& m) const
{
    auto const v1 = this->indices().at(0);
    auto const v2 = this->indices().at(1);
    auto const v3 = this->indices().at(2);
    auto const v4 = this->indices().at(3);

    auto const p1 = p.row(v1);
    auto const p2 = p.row(v2);
    auto const p3 = p.row(v3);
    auto const p4 = p.row(v4);

    auto const w1 = m(v1);
    auto const w2 = m(v2);
    auto const w3 = m(v3);
    auto const w4 = m(v4);

    // auto const Vsigned         = signed_volume(p);
    // auto const V               = std::abs(Vsigned);
    // bool const is_tet_inverted = Vsigned <= 0.0;

    Eigen::Matrix3d Ds;
    Ds.col(0) = (p1 - p4).transpose();
    Ds.col(1) = (p2 - p4).transpose();
    Ds.col(2) = (p3 - p4).transpose();

    Eigen::Matrix3d const F = Ds * DmInv_;
    Eigen::Matrix3d const I = Eigen::Matrix3d::Identity();

    scalar_type constexpr epsilon = 1e-20;

    Eigen::Matrix3d Piola;
    scalar_type psi{};

    /*if (is_tet_inverted)
    {
        Eigen::JacobiSVD<Eigen::Matrix3d> UFhatV(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Vector3d const Fsigma = UFhatV.singularValues();
        Eigen::Matrix3d Fhat;
        Fhat.setZero();
        Fhat(0, 0) = Fsigma(0);
        Fhat(1, 1) = Fsigma(1);
        Fhat(2, 2) = Fsigma(2);

        Eigen::Matrix3d U       = UFhatV.matrixU();
        Eigen::Matrix3d const V = UFhatV.matrixV();

        auto smallest_element_idx = 0;
        if (Fsigma(0) < Fsigma(1) && Fsigma(0) < Fsigma(2))
            smallest_element_idx = 0;
        if (Fsigma(1) < Fsigma(0) && Fsigma(1) < Fsigma(2))
            smallest_element_idx = 1;
        if (Fsigma(2) < Fsigma(0) && Fsigma(2) < Fsigma(1))
            smallest_element_idx = 2;

        Fhat(smallest_element_idx, smallest_element_idx) =
            -Fhat(smallest_element_idx, smallest_element_idx);
        U.col(smallest_element_idx) = -U.col(smallest_element_idx);

        // stress reaches maximum at 58% compression
        scalar_type constexpr min_singular_value = 0.577;
        Fhat(0, 0)                               = std::min(Fhat(0, 0), min_singular_value);
        Fhat(1, 1)                               = std::min(Fhat(1, 1), min_singular_value);
        Fhat(2, 2)                               = std::min(Fhat(2, 2), min_singular_value);

        Eigen::Matrix3d const Ehat     = 0.5 * (Fhat.transpose() * Fhat - I);
        scalar_type const trace        = Ehat.trace();
        Eigen::Matrix3d const Piolahat = Fhat * ((2. * mu_ * Ehat) + (lambda_ * trace * I));

        Eigen::Matrix3d const E = U * Ehat * V.transpose();
        psi = mu_ * (E.array() * E.array()).sum() + 0.5 * lambda_ * E.trace() * E.trace();

        Piola = U * Piolahat * V.transpose();
    }
    else*/
    //{
        Eigen::Matrix3d const E = 0.5 * (F.transpose() * F - I);

        scalar_type const trace = E.trace();
        psi   = mu_ * (E.array() * E.array()).sum() + 0.5 * lambda_ * trace * trace;
        Piola = F * ((2. * mu_ * E) + (lambda_ * E.trace() * I));
    //}

    // We can omit the negative sign of the H matrix,
    // because the position based corrections are already
    // made in the negative gradient direction
    Eigen::Matrix3d const H  = - V0_ * Piola * DmInv_.transpose();
    Eigen::Vector3d const f1 = H.col(0);
    Eigen::Vector3d const f2 = H.col(1);
    Eigen::Vector3d const f3 = H.col(2);
    Eigen::Vector3d const f4 = -(f1 + f2 + f3);

    // clang-format off
     auto const lambda =
        w1 * f1.squaredNorm() +
        w2 * f2.squaredNorm() +
        w3 * f3.squaredNorm() +
        w4 * f4.squaredNorm();
    // clang-format on

    if (lambda < epsilon)
        return;

    scalar_type const C = V0_ * psi;
    auto const s        = C / lambda;

    p.row(v1) += s * w1 * f1;
    p.row(v2) += s * w2 * f2;
    p.row(v3) += s * w3 * f3;
    p.row(v4) += s * w4 * f4;
}

} // namespace xpbd