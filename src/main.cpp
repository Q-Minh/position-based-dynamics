#include <GLFW/glfw3.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>
#include <igl/unproject.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/unproject_ray.h>

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "Usage: ./pbd <path to obj mesh>\n";
        return 1;
    }

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::readOBJ(argv[1], V, F);

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

    Eigen::RowVector3d vbox_mean = Vbox.colwise().mean();
    Eigen::RowVector3d v_mean    = V.colwise().mean();
    Vbox.rowwise() -= vbox_mean;
    V.rowwise() -= v_mean;
    V.rowwise() += Eigen::RowVector3d(0., .15, 0.);

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.data().point_size = 10.f;
    viewer.core().align_camera_center(V);
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

    draw_floor_points();
    draw_floor_edges();

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    struct picking_state_t
    {
        bool is_picking = false;
        int vertex      = 0;
        float speed     = 0.001f;
        int mouse_x, mouse_y;
    } picking_state;

    viewer.callback_mouse_down =
        [&](igl::opengl::glfw::Viewer& viewer, int button, int modifier) -> bool {
        using button_type = igl::opengl::glfw::Viewer::MouseButton;
        if (static_cast<button_type>(button) != button_type::Left)
            return false;

        if (modifier != GLFW_MOD_CONTROL)
            return false;

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
            V,
            F,
            fid,
            bc);

        if (!hit)
            return false;

        Eigen::Vector3i const face{F(fid, 0), F(fid, 1), F(fid, 2)};
        unsigned int closest_vertex = face(0);

        if (bc(1) > bc(0) && bc(1) > bc(2))
        {
            closest_vertex = face(1);
        }
        else if (bc(2) > bc(0) && bc(2) > bc(1))
        {
            closest_vertex = face(2);
        }

        viewer.data().clear_points();
        draw_floor_points();
        viewer.data().add_points(V.row(closest_vertex), Eigen::RowVector3d(1., 0., 0.));
        picking_state.is_picking = true;
        picking_state.vertex     = closest_vertex;

        return true;
    };

    viewer.callback_mouse_move =
        [&](igl::opengl::glfw::Viewer& viewer, int button, int modifier) -> bool {
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

        V.row(picking_state.vertex) += direction * static_cast<double>(picking_state.speed);
        viewer.data().set_mesh(V, F);

        viewer.data().clear_points();
        draw_floor_points();
        viewer.data().add_points(V.row(picking_state.vertex), Eigen::RowVector3d(1., 0., 0.));

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

    menu.callback_draw_custom_window = [&]() {
        ImGui::Begin("Position Based Dynamics");
        if (ImGui::CollapsingHeader("Picking", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::InputFloat("Speed", &picking_state.speed, 0.0001f);
        }
        ImGui::End();
    };

    viewer.launch();

    return 0;
}