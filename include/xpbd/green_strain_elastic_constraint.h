#ifndef PBD_GREEN_STRAIN_ELASTIC_CONSTRAINT_H
#define PBD_GREEN_STRAIN_ELASTIC_CONSTRAINT_H

#include "xpbd/constraint.h"

namespace xpbd {

class green_strain_elastic_constraint_t : public constraint_t
{
  public:
    using self_type      = green_strain_elastic_constraint_t;
    using base_type      = constraint_t;
    using index_type     = std::uint32_t;
    using scalar_type    = typename constraint_t::scalar_type;
    using masses_type    = Eigen::VectorXd;
    using positions_type = typename base_type::positions_type;
    using gradient_type  = typename base_type::gradient_type;

  public:
    green_strain_elastic_constraint_t(
        std::initializer_list<index_type> indices,
        positions_type const& p,
        scalar_type young_modulus,
        scalar_type poisson_ratio);

    virtual void project(positions_type& p, masses_type const& m) const override;

  protected:
    scalar_type signed_volume(positions_type const& V) const;
    void project_inverted(positions_type& p, masses_type const& m) const;

  private:
    scalar_type V0_;
    Eigen::Matrix3d DmInv_;
    scalar_type mu_;
    scalar_type lambda_;
};

} // namespace xpbd

#endif // PBD_GREEN_STRAIN_ELASTIC_CONSTRAINT_H