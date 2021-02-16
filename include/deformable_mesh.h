#ifndef PBD_DEFORMABLE_MESH_H
#define PBD_DEFORMABLE_MESH_H

#include <Eigen/Core>
#include <constraint.h>

namespace pbd {

class deformable_mesh_t
{
    using positions_type   = Eigen::MatrixXd;
    using masses_type      = Eigen::VectorXd;
    using velocities_type  = Eigen::MatrixX3d;
    using faces_type       = Eigen::MatrixXi;
    using elements_type    = Eigen::MatrixXi;
    using constraints_type = std::vector<std::unique_ptr<constraint_t>>;

  public:
    deformable_mesh_t() = default;

    deformable_mesh_t(
        positions_type const& positions,
        faces_type const& faces,
        elements_type const& elements,
        masses_type const& masses)
        : p_(positions),
          F_(faces),
          E_(elements),
          m_(masses),
          v_(positions.rows(), positions.cols()),
          constraints_{},
          fixed_(positions.rows(), false)
    {
    }

    deformable_mesh_t(
        positions_type const& positions,
        faces_type const& faces,
        elements_type const& elements)
        : p_(positions),
          F_(faces),
          E_(elements),
          m_(positions.rows()),
          v_(positions.rows(), positions.cols()),
          constraints_{},
          fixed_(positions.rows(), false)
    {
        m_.setOnes();
        v_.setZero();
    }

    bool is_fixed(int i) const { return fixed_[i]; }
    void fix(int i) { fixed_[i] = true; }
    void unfix(int i) { fixed_[i] = false; }
    void toggle_fixed(int i) { fixed_[i] = !fixed_[i]; }

    positions_type const& positions() const { return p_; }
    faces_type const& faces() const { return F_; }
    elements_type const& elements() const { return E_; }
    constraints_type const& constraints() const { return constraints_; }
    velocities_type const& velocity() const { return v_; }
    masses_type const& mass() const { return m_; }
    std::vector<bool> const& fixed() const { return fixed_; }

    positions_type& positions() { return p_; }
    faces_type& faces() { return F_; }
    elements_type& elements() { return E_; }
    constraints_type& constraints() { return constraints_; }
    velocities_type& velocity() { return v_; }
    masses_type& mass() { return m_; }
    std::vector<bool>& fixed() { return fixed_; }

    void constrain_edge_lengths();

  private:
    positions_type p_;             ///< Positions
    faces_type F_;                 ///< Faces
    elements_type E_;              ///< Elements
    masses_type m_;                ///< Per-vertex mass
    velocities_type v_;            ///< Per-vertex velocity
    constraints_type constraints_; ///< PBD constraints
    std::vector<bool> fixed_;      ///< Flags fixed positions
};

} // namespace pbd

#endif // PBD_DEFORMABLE_MESH_H