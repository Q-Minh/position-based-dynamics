#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOBJ.h>

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "Usage: ./pbd <path to ply mesh>\n";
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
		1., .1, 0.,
		0., .1, 0.,
		0., 0., 1.,
		1., 0., 1.,
		1., .1, 1.,
		0., .1, 1.;

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
    viewer.data().add_points(Vbox, Eigen::RowVector3d(1, 0, 0));
    viewer.data().point_size = 10.f;
    viewer.core().align_camera_center(Vbox);
    viewer.data().show_lines = false;
    for (auto i = 0u; i < Ebox.rows(); ++i)
    {
        viewer.data().add_edges(
            Vbox.row(Ebox(i, 0)),
            Vbox.row(Ebox(i, 1)),
            Eigen::RowVector3d(1, 0, 0));
    }

    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    menu.callback_draw_custom_window = [&]() {
        ImGui::Begin("Position Based Dynamics");
        ImGui::End();
    };

    viewer.launch();

    return 0;
}