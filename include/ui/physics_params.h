#ifndef PBD_UI_PHYSICS_PARAMS_H
#define PBD_UI_PHYSICS_PARAMS_H

namespace ui {

struct physics_params_t
{
    bool is_gravity_active = false;
    float dt               = 0.166667;
    int solver_iterations  = 10;
    int solver_substeps    = 10;

    // fem
    float young_modulus = 1.f;
    float poisson_ratio = 0.3f;
};

} // namespace ui

#endif // PBD_UI_PHYSICS_PARAMS_H