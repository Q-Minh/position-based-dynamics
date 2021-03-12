#ifndef PBD_CONSTRAINT_H
#define PBD_CONSTRAINT_H

#include <Eigen/Core>
#include <vector>

namespace xpbd {

class constraint_t
{
  public:
    using index_type     = std::uint32_t;
    using masses_type    = Eigen::VectorXd;
    using positions_type = Eigen::MatrixXd;
    using gradient_type  = Eigen::Vector3d;
    using scalar_type    = double;

  public:
    constraint_t(std::initializer_list<index_type> indices) : indices_(indices), alpha_(0.0) {}
    constraint_t(std::initializer_list<index_type> indices, scalar_type const alpha)
        : indices_(indices), alpha_(alpha)
    {
    }

    virtual void
    project(positions_type& V, masses_type const& M, scalar_type& lagrange, scalar_type const dt)
        const = 0;
    std::vector<index_type> const& indices() const { return indices_; }

  protected:
    scalar_type alpha_;

  private:
    std::vector<index_type> indices_;
};

} // namespace xpbd

#endif // PBD_CONSTRAINT_H