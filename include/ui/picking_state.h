#ifndef PBD_UI_PICKING_STATE_H
#define PBD_UI_PICKING_STATE_H

namespace ui {

struct picking_state_t
{
    bool is_picking = false;
    int vertex      = 0;
    float force     = 400.f;
    int mouse_x, mouse_y;
};

} // namespace ui

#endif // PBD_UI_PICKING_STATE_H