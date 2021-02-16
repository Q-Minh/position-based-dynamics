#ifndef PBD_DEFORMABLE_MESH_H
#define PBD_DEFORMABLE_MESH_H

#include <Eigen/Core>
#include <constraint.h>

namespace pbd {

class deformable_mesh_t
{
    using vertices_type    = Eigen::MatrixX3d;
    using faces_type       = Eigen::MatrixX3i;
    using elements_type    = Eigen::MatrixX4i;
    using constraints_type = std::vector<constraint_t>;

  public:
    vertices_type const& vertices() const { return V_; }
    faces_type const& faces() const { return F_; }
    elements_type const& elements() const { return E_; }
    constraints_type const& constraints() const { return constraints_; }

    vertices_type& vertices() { return V_; }
    faces_type& faces() { return F_; }
    elements_type& elements() { return E_; }
    constraints_type& constraints() { return constraints_; }

  private:
    vertices_type V_;              ///< Vertices
    faces_type F_;                 ///< Faces
    elements_type E_;              ///< Elements
    constraints_type constraints_; ///< PBD constraints
};

} // namespace pbd

#endif // PBD_DEFORMABLE_MESH_H