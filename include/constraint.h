#ifndef PBD_CONSTRAINT_H
#define PBD_CONSTRAINT_H

#include <Eigen/Core>
#include <vector>

namespace pbd {

class constraint_t
{
  public:
    using index_type    = std::uint32_t;
    using masses_type   = Eigen::VectorXd;
    using vertices_type = Eigen::MatrixX3d;
    using gradient_type = Eigen::Vector3d;
    using scalar_type   = double;

  public:
    constraint_t(std::initializer_list<index_type> indices) : indices_(indices) {}

    virtual scalar_type evaluate(vertices_type const& V, masses_type const& M) const = 0;
    virtual void project(vertices_type& V, masses_type const& M) const               = 0;
    virtual std::vector<index_type> const& indices() const { return indices_; }

  private:
    std::vector<index_type> indices_;
};

} // namespace pbd

#endif // PBD_CONSTRAINT_H
