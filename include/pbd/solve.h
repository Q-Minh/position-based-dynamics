#ifndef PBD_SOLVE_H
#define PBD_SOLVE_H

#include "deformable_mesh.h"

namespace pbd {

void solve(
    deformable_mesh_t& model,
    Eigen::MatrixX3d const& fext,
    double dt                = 0.01,
    std::uint32_t iterations = 10, 
    std::uint32_t substeps = 10);

} // namespace pbd

#endif // PBD_SOLVE_H
