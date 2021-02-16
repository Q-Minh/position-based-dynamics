#ifndef PBD_GET_SIMPLE_CLOTH_MODEL_H
#define PBD_GET_SIMPLE_CLOTH_MODEL_H

#include <Eigen/Core>

namespace geometry {

std::pair<Eigen::MatrixXd, Eigen::MatrixXi> get_simple_cloth_model(int rows, int cols)
{
    std::vector<Eigen::RowVector3d> cloth_positions;
    std::vector<Eigen::RowVector3i> cloth_faces;

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            auto const xoffset = static_cast<double>(i);
            auto const yoffset = static_cast<double>(j);
            cloth_positions.push_back(Eigen::RowVector3d{xoffset, yoffset, 0.});

            if (i == rows - 1u || j == cols - 1u)
                continue;

            auto const upper_left_corner  = i * cols + j;
            auto const upper_right_corner = i * cols + (j + 1);
            auto const lower_left_corner  = (i + 1) * cols + j;
            auto const lower_right_corner = (i + 1) * cols + (j + 1);

            cloth_faces.push_back(
                Eigen::RowVector3i{lower_left_corner, upper_right_corner, upper_left_corner});
            cloth_faces.push_back(
                Eigen::RowVector3i{lower_left_corner, lower_right_corner, upper_right_corner});
        }
    }

    Eigen::MatrixXd VCloth(cloth_positions.size(), 3);
    Eigen::MatrixXi FCloth(cloth_faces.size(), 3);
    for (auto i = 0u; i < cloth_positions.size(); ++i)
        VCloth.row(i) = cloth_positions[i];
    for (auto i = 0u; i < cloth_faces.size(); ++i)
        FCloth.row(i) = cloth_faces[i];
}

} // namespace geometry

#endif // PBD_GET_SIMPLE_CLOTH_MODEL_H
