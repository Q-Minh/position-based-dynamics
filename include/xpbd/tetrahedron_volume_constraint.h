#ifndef PBD_TETRAHEDRON_VOLUME_CONSTRAINT_H
#define PBD_TETRAHEDRON_VOLUME_CONSTRAINT_H

#include "constraint.h"

namespace xpbd {

class tetrahedron_volume_constraint_t : public constraint_t
{
  public:
    using self_type      = tetrahedron_volume_constraint_t;
    using base_type      = constraint_t;
    using index_type     = std::uint32_t;
    using scalar_type    = double;
    using masses_type    = Eigen::VectorXd;
    using positions_type = typename base_type::positions_type;
    using gradient_type  = typename base_type::gradient_type;

  public:
    tetrahedron_volume_constraint_t(
        std::initializer_list<index_type> indices,
        positions_type const& p,
        scalar_type const alpha = 0.0);

    scalar_type evaluate(positions_type const& V, masses_type const& M) const;
    virtual void
    project(positions_type& V, masses_type const& M, scalar_type& lagrange, scalar_type const dt)
        const override;

  protected:
    scalar_type volume(positions_type const& V) const;

  private:
    scalar_type V0_;
};

} // namespace xpbd

#endif // PBD_TETRAHEDRON_VOLUME_CONSTRAINT_H
