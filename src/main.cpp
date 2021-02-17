#include "deformable_mesh.h"
#include "edge_length_constraint.h"
#include "solve.h"

#include <GLFW/glfw3.h>
#include <chrono>
#include <filesystem>
#include <igl/file_dialog_open.h>
#include <igl/marching_tets.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>
#include <igl/unproject.h>
#include <igl/unproject_onto_mesh.h>

int main(int argc, char** argv)
{
    // Simulation state
    pbd::deformable_mesh_t model{};
    Eigen::MatrixX3d fext;

    auto const is_model_ready = [&]() {
        return model.positions().rows() > 0;
    };

    Eigen::MatrixXd Vbox(8, 3);
    Eigen::MatrixXi Ebox(12, 2);

    // clang-format off
	Vbox <<
		0., 0., 0.,
		1., 0., 0.,
		1., .01, 0.,
		0., .01, 0.,
		0., 0., 1.,
		1., 0., 1.,
		1., .01, 1.,
		0., .01, 1.;

	Ebox <<
		0, 1,
		1, 2,
		2, 3,
		3, 0,
		4, 5,
		5, 6,
		6, 7,
		7, 4,
		0, 4,
		1, 5,
		2, 6,
		7, 3;
    // clang-format on

    igl::opengl::glfw::Viewer viewer;
    viewer.data().point_size = 10.f;
    viewer.data().show_lines = false;

    auto const draw_floor_points = [&]() {
        viewer.data().add_points(Vbox, Eigen::RowVector3d(1, 0, 0));
    };
    auto const draw_floor_edges = [&]() {
        for (auto i = 0u; i < Ebox.rows(); ++i)
        {
            viewer.data().add_edges(
                Vbox.row(Ebox(i, 0)),
                Vbox.row(Ebox(i, 1)),
                Eigen::RowVector3d(1, 0, 0));
        }
    };

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    struct picking_state_t
    {
        bool is_picking = false;
        int vertex      = 0;
        float force     = 400.f;
        int mouse_x, mouse_y;
    } picking_state;

    viewer.callback_mouse_down =
        [&](igl::opengl::glfw::Viewer& viewer, int button, int modifier) -> bool {
        if (!is_model_ready())
            return false;

        using button_type = igl::opengl::glfw::Viewer::MouseButton;
        if (static_cast<button_type>(button) != button_type::Left)
            return false;

        bool const process_pick = modifier == GLFW_MOD_CONTROL || modifier == GLFW_MOD_SHIFT;

        int fid;
        double const x = static_cast<double>(viewer.current_mouse_x);
        double const y = viewer.core().viewport(3) - static_cast<double>(viewer.current_mouse_y);
        picking_state.mouse_x = viewer.current_mouse_x;
        picking_state.mouse_y = viewer.current_mouse_y;

        Eigen::Vector3f bc{};

        bool const hit = igl::unproject_onto_mesh(
            Eigen::Vector2f(x, y),
            viewer.core().view,
            viewer.core().proj,
            viewer.core().viewport,
            model.positions(),
            model.faces(),
            fid,
            bc);

        if (!hit)
            return false;

        Eigen::Vector3i const face{
            model.faces()(fid, 0),
            model.faces()(fid, 1),
            model.faces()(fid, 2)};
        unsigned int closest_vertex = face(0);

        if (bc(1) > bc(0) && bc(1) > bc(2))
        {
            closest_vertex = face(1);
        }
        else if (bc(2) > bc(0) && bc(2) > bc(1))
        {
            closest_vertex = face(2);
        }

        if (modifier == GLFW_MOD_CONTROL)
        {
            picking_state.is_picking = true;
            picking_state.vertex     = closest_vertex;
        }
        if (modifier == GLFW_MOD_SHIFT)
        {
            model.toggle_fixed(closest_vertex);
        }

        return process_pick;
    };

    viewer.callback_mouse_move =
        [&](igl::opengl::glfw::Viewer& viewer, int button, int modifier) -> bool {
        if (!is_model_ready())
            return false;

        if (!picking_state.is_picking)
            return false;

        double const x1 = static_cast<double>(picking_state.mouse_x);
        double const y1 = viewer.core().viewport(3) - static_cast<double>(picking_state.mouse_y);

        double const x2 = static_cast<double>(viewer.current_mouse_x);
        double const y2 = viewer.core().viewport(3) - static_cast<double>(viewer.current_mouse_y);

        Eigen::Vector3d const p1 = igl::unproject(
                                       Eigen::Vector3f(x1, y1, .5f),
                                       viewer.core().view,
                                       viewer.core().proj,
                                       viewer.core().viewport)
                                       .cast<double>();
        ;
        Eigen::Vector3d const p2 = igl::unproject(
                                       Eigen::Vector3f(x2, y2, .5f),
                                       viewer.core().view,
                                       viewer.core().proj,
                                       viewer.core().viewport)
                                       .cast<double>();

        Eigen::Vector3d const direction = (p2 - p1).normalized();

        fext.row(picking_state.vertex) =
            direction.transpose() * static_cast<double>(picking_state.force);

        viewer.data().add_points(
            model.positions().row(picking_state.vertex),
            Eigen::RowVector3d(1., 0., 0.));

        picking_state.mouse_x = viewer.current_mouse_x;
        picking_state.mouse_y = viewer.current_mouse_y;

        return true;
    };

    viewer.callback_mouse_up =
        [&](igl::opengl::glfw::Viewer& viewer, int button, int modifier) -> bool {
        if (picking_state.is_picking)
            picking_state.is_picking = false;

        return false;
    };

    float dt = 0.166667;
    menu.callback_draw_viewer_window =
        [&]() {
            ImGui::SetNextWindowSize(ImVec2(300.0f, 480.0f), ImGuiSetCond_FirstUseEver);
            ImGui::Begin("Position Based Dynamics");

            float w = ImGui::GetContentRegionAvailWidth();
            float p = ImGui::GetStyle().FramePadding.x;

            if (ImGui::Button("Load##Mesh", ImVec2((w - p) / 2.f, 0)))
            {
                std::string const filename = igl::file_dialog_open();
                std::filesystem::path const mesh{filename};
                if (std::filesystem::exists(mesh) && std::filesystem::is_regular_file(mesh))
                {
                    Eigen::MatrixXd V;
                    Eigen::MatrixXi F;
                    igl::readOBJ(mesh.string(), V, F);

                    Eigen::RowVector3d vbox_mean = Vbox.colwise().mean();
                    Eigen::RowVector3d v_mean    = V.colwise().mean();
                    Vbox.rowwise() -= vbox_mean;
                    V.rowwise() -= v_mean;
                    V.rowwise() += Eigen::RowVector3d(0., .15, 0.);

                    model = pbd::deformable_mesh_t{
                    V,
                    F,
                    F /* When we'll work with tet meshes, we will give tets as argument here */};

                    model.constrain_edge_lengths();

                    fext.resizeLike(model.positions());
                    fext.setZero();

                    viewer.data().set_mesh(model.positions(), model.faces());
                    viewer.core().align_camera_center(model.positions());
                }
            }
            ImGui::Checkbox("Simulate", &viewer.core().is_animating);
            ImGui::InputFloat("Timestep", &dt, 0.01f, 0.1f, "%.4f");
            if (ImGui::CollapsingHeader("Picking", ImGuiTreeNodeFlags_DefaultOpen))
            {
                ImGui::InputFloat(
                    "Dragging force",
                    &picking_state.force,
                    1.f,
                    10.f,
                    "%.1f");
            }
            ImGui::End();
        };

    auto const draw_fixed_points = [&]() {
        for (auto i = 0u; i < model.positions().rows(); ++i)
        {
            if (!model.is_fixed(i))
                continue;

            viewer.data().add_points(model.positions().row(i), Eigen::RowVector3d{1., 0., 0.});
        }
    };

    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer& viewer) -> bool {
        if (!is_model_ready())
            return false;

        if (viewer.core().is_animating)
        {
            fext.col(1).array() -= 9.81;
            pbd::solve(model, fext, dt);
            fext.setZero();
        }

        viewer.data().clear();
        viewer.data().set_mesh(model.positions(), model.faces());
        draw_fixed_points();
        //draw_floor_points();
        //draw_floor_edges();

        return false; // do not return from drawing loop
    };

    viewer.core().is_animating = false;
    viewer.launch();

    return 0;
}