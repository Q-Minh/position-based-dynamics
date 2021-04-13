#include "xpbd/green_strain_elastic_constraint.h"

#include <Eigen/Dense>
#include <eigen/SVD>

namespace xpbd {

green_strain_elastic_constraint_t::green_strain_elastic_constraint_t(
    std::initializer_list<index_type> indices,
    positions_type const& p,
    scalar_type young_modulus,
    scalar_type poisson_ratio,
    scalar_type const alpha)
    : base_type(indices, alpha), V0_{0.}, DmInv_{}, mu_{}, lambda_{}
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

    V0_     = (1. / 6.) * Dm.determinant();
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

    Eigen::Matrix3d Ds;
    Ds.col(0)      = (p1 - p4).transpose();
    Ds.col(1)      = (p2 - p4).transpose();
    Ds.col(2)      = (p3 - p4).transpose();
    auto const vol = (1. / 6.) * Ds.determinant();
    return vol;
}

void green_strain_elastic_constraint_t::project(
    positions_type& p,
    masses_type const& m,
    scalar_type& lagrange,
    scalar_type const dt) const
{
    auto const v1 = this->indices().at(0);
    auto const v2 = this->indices().at(1);
    auto const v3 = this->indices().at(2);
    auto const v4 = this->indices().at(3);

    auto const p1 = p.row(v1);
    auto const p2 = p.row(v2);
    auto const p3 = p.row(v3);
    auto const p4 = p.row(v4);

    auto const w1 = 1. / m(v1);
    auto const w2 = 1. / m(v2);
    auto const w3 = 1. / m(v3);
    auto const w4 = 1. / m(v4);

    auto const Vsigned        = signed_volume(p);
    bool const is_V_positive  = Vsigned >= 0.;
    bool const is_V0_positive = V0_ >= 0.;
    bool const is_tet_inverted =
        (is_V_positive && !is_V0_positive) || (!is_V_positive && is_V0_positive);
    scalar_type constexpr epsilon = 1e-20;

    Eigen::Matrix3d Ds;
    Ds.col(0) = (p1 - p4).transpose();
    Ds.col(1) = (p2 - p4).transpose();
    Ds.col(2) = (p3 - p4).transpose();

    Eigen::Matrix3d const F = Ds * DmInv_;
    Eigen::Matrix3d const I = Eigen::Matrix3d::Identity();

    // TODO: Implement correct inversion handling described in
    // Irving, Geoffrey, Joseph Teran, and Ronald Fedkiw. "Invertible finite elements for robust
    // simulation of large deformation." Proceedings of the 2004 ACM SIGGRAPH/Eurographics symposium
    // on Computer animation. 2004.
    // scalar_type psi{};
    // Eigen::Matrix3d Piola;
    Eigen::JacobiSVD<Eigen::Matrix3d> UFhatV(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Vector3d const Fsigma = UFhatV.singularValues();
    Eigen::Matrix3d Fhat;
    Fhat.setZero();
    Fhat(0, 0) = Fsigma(0);
    Fhat(1, 1) = Fsigma(1);
    Fhat(2, 2) = Fsigma(2);

    Eigen::Matrix3d U       = UFhatV.matrixU();
    Eigen::Matrix3d const V = UFhatV.matrixV();

    if (is_tet_inverted)
    {
        Fhat(2, 2) = -Fhat(2, 2);
        U.col(2)   = -U.col(2);
    }

    // stress reaches maximum at 58% compression
    scalar_type constexpr min_singular_value = 0.577;
    Fhat(0, 0)                               = std::max(Fhat(0, 0), min_singular_value);
    Fhat(1, 1)                               = std::max(Fhat(1, 1), min_singular_value);
    Fhat(2, 2)                               = std::max(Fhat(2, 2), min_singular_value);

    Eigen::Matrix3d const Ehat     = 0.5 * (Fhat.transpose() * Fhat - I);
    scalar_type const EhatTrace    = Ehat.trace();
    Eigen::Matrix3d const Piolahat = Fhat * ((2. * mu_ * Ehat) + (lambda_ * EhatTrace * I));

    Eigen::Matrix3d const E  = U * Ehat * V.transpose();
    scalar_type const Etrace = E.trace();
    scalar_type const psi = mu_ * (E.array() * E.array()).sum() + 0.5 * lambda_ * Etrace * Etrace;

    Eigen::Matrix3d const Piola = U * Piolahat * V.transpose();

    // H is the negative gradient of the elastic potential
    scalar_type const V0     = std::abs(V0_);
    Eigen::Matrix3d const H  = -V0 * Piola * DmInv_.transpose();
    Eigen::Vector3d const f1 = H.col(0);
    Eigen::Vector3d const f2 = H.col(1);
    Eigen::Vector3d const f3 = H.col(2);
    Eigen::Vector3d const f4 = -(f1 + f2 + f3);

    // clang-format off
     auto const weighted_sum_of_gradients =
        w1 * f1.squaredNorm() +
        w2 * f2.squaredNorm() +
        w3 * f3.squaredNorm() +
        w4 * f4.squaredNorm();
    // clang-format on

    if (weighted_sum_of_gradients < epsilon)
        return;

    scalar_type const C           = V0 * psi;
    scalar_type const alpha_tilde = alpha_ / (dt * dt);
    scalar_type const delta_lagrange =
        -(C + alpha_tilde * lagrange) / (weighted_sum_of_gradients + alpha_tilde);

    lagrange += delta_lagrange;
    // because f = - grad(potential), then grad(potential) = -f and thus grad(C) = -f
    p.row(v1) += w1 * -f1 * delta_lagrange;
    p.row(v2) += w2 * -f2 * delta_lagrange;
    p.row(v3) += w3 * -f3 * delta_lagrange;
    p.row(v4) += w4 * -f4 * delta_lagrange;
}

green_strain_elastic_constraint_t::scalar_type
green_strain_elastic_constraint_t::evaluate(positions_type const& p, masses_type const& m) const
{
    auto const v1 = this->indices().at(0);
    auto const v2 = this->indices().at(1);
    auto const v3 = this->indices().at(2);
    auto const v4 = this->indices().at(3);

    auto const p1 = p.row(v1);
    auto const p2 = p.row(v2);
    auto const p3 = p.row(v3);
    auto const p4 = p.row(v4);

    auto const w1 = 1. / m(v1);
    auto const w2 = 1. / m(v2);
    auto const w3 = 1. / m(v3);
    auto const w4 = 1. / m(v4);

    auto const Vsigned        = signed_volume(p);
    bool const is_V_positive  = Vsigned >= 0.;
    bool const is_V0_positive = V0_ >= 0.;
    bool const is_tet_inverted =
        (is_V_positive && !is_V0_positive) || (!is_V_positive && is_V0_positive);
    scalar_type constexpr epsilon = 1e-20;

    Eigen::Matrix3d Ds;
    Ds.col(0) = (p1 - p4).transpose();
    Ds.col(1) = (p2 - p4).transpose();
    Ds.col(2) = (p3 - p4).transpose();

    Eigen::Matrix3d const F = Ds * DmInv_;
    Eigen::Matrix3d const I = Eigen::Matrix3d::Identity();

    // TODO: Implement correct inversion handling described in
    // Irving, Geoffrey, Joseph Teran, and Ronald Fedkiw. "Invertible finite elements for robust
    // simulation of large deformation." Proceedings of the 2004 ACM SIGGRAPH/Eurographics symposium
    // on Computer animation. 2004.
    // scalar_type psi{};
    // Eigen::Matrix3d Piola;
    Eigen::JacobiSVD<Eigen::Matrix3d> UFhatV(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Vector3d const Fsigma = UFhatV.singularValues();
    Eigen::Matrix3d Fhat;
    Fhat.setZero();
    Fhat(0, 0) = Fsigma(0);
    Fhat(1, 1) = Fsigma(1);
    Fhat(2, 2) = Fsigma(2);

    Eigen::Matrix3d U       = UFhatV.matrixU();
    Eigen::Matrix3d const V = UFhatV.matrixV();

    if (is_tet_inverted)
    {
        Fhat(2, 2) = -Fhat(2, 2);
        U.col(2)   = -U.col(2);
    }

    // stress reaches maximum at 58% compression
    scalar_type constexpr min_singular_value = 0.577;
    Fhat(0, 0)                               = std::max(Fhat(0, 0), min_singular_value);
    Fhat(1, 1)                               = std::max(Fhat(1, 1), min_singular_value);
    Fhat(2, 2)                               = std::max(Fhat(2, 2), min_singular_value);

    Eigen::Matrix3d const Ehat     = 0.5 * (Fhat.transpose() * Fhat - I);
    scalar_type const EhatTrace    = Ehat.trace();
    Eigen::Matrix3d const Piolahat = Fhat * ((2. * mu_ * Ehat) + (lambda_ * EhatTrace * I));

    Eigen::Matrix3d const E  = U * Ehat * V.transpose();
    scalar_type const Etrace = E.trace();
    scalar_type const psi = mu_ * (E.array() * E.array()).sum() + 0.5 * lambda_ * Etrace * Etrace;

    scalar_type const V0 = std::abs(V0_);
    scalar_type const C = V0 * psi;
    return C;
}

} // namespace xpbd