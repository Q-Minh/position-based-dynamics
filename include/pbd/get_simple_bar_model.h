#ifndef PBD_GET_SIMPLE_BAR_MODEL_H
#define PBD_GET_SIMPLE_BAR_MODEL_H

#include <Eigen/Core>
#include <igl/boundary_facets.h>
#include <tuple>

namespace geometry {

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXi>
get_simple_bar_model(std::size_t width, std::size_t height, std::size_t depth)
{
    Eigen::MatrixXd V;
    V.resize(width * height * depth, 3);

    for (auto i = 0; i < width; ++i)
    {
        for (auto j = 0; j < height; ++j)
        {
            for (auto k = 0; k < depth; ++k)
            {
                auto const row = i * height * depth + j * depth + k;
                V.row(row)     = Eigen::Vector3d(
                    static_cast<double>(i),
                    static_cast<double>(j),
                    static_cast<double>(k));
            }
        }
    }

    auto const tet_count = (width - 1u) * (height - 1u) * (depth - 1u) * 5u;
    Eigen::MatrixXi T;
    T.resize(tet_count, 4);
    for (std::size_t i = 0; i < width - 1u; i++)
    {
        for (std::size_t j = 0; j < height - 1u; j++)
        {
            for (std::size_t k = 0; k < depth - 1u; k++)
            {
                //     7*-----*6
                //     /|    /|
                //    / |   / |
                //  4*-----*5 |
                //   | 3*--|--*2
                //   | /   | /
                //   |/    |/
                //  0*-----*1

                // clang-format off
                int p0 = i        * height * depth + j        * depth + k       ;
                int p1 = (i + 1u) * height * depth + j        * depth + k       ;
                int p2 = (i + 1u) * height * depth + (j + 1u) * depth + k       ;
                int p3 = i        * height * depth + (j + 1u) * depth + k       ;
                int p4 = i        * height * depth + j        * depth + (k + 1u);
                int p5 = (i + 1u) * height * depth + j        * depth + (k + 1u);
                int p6 = (i + 1u) * height * depth + (j + 1u) * depth + (k + 1u);
                int p7 = i        * height * depth + (j + 1u) * depth + (k + 1u);
                // clang-format on

                auto const row = (i * (height - 1u) * (depth - 1u) + j * (depth - 1u) + k) * 5u;
                if ((i + j + k) % 2 == 1)
                {
                    T.row(row + 0u) = Eigen::RowVector4i{p1, p0, p5, p2};
                    T.row(row + 1u) = Eigen::RowVector4i{p5, p2, p7, p6};
                    T.row(row + 2u) = Eigen::RowVector4i{p7, p0, p5, p4};
                    T.row(row + 3u) = Eigen::RowVector4i{p2, p0, p7, p3};
                    T.row(row + 4u) = Eigen::RowVector4i{p5, p0, p7, p2};
                }
                else
                {
                    T.row(row + 0u) = Eigen::RowVector4i{p3, p1, p4, p0};
                    T.row(row + 1u) = Eigen::RowVector4i{p6, p1, p3, p2};
                    T.row(row + 2u) = Eigen::RowVector4i{p4, p1, p6, p5};
                    T.row(row + 3u) = Eigen::RowVector4i{p6, p3, p4, p7};
                    T.row(row + 4u) = Eigen::RowVector4i{p3, p1, p6, p4};
                }
            }
        }
    }

    Eigen::MatrixXi F;
    igl::boundary_facets(T, F);
    T = T.rowwise().reverse().eval();
    F = F.rowwise().reverse().eval();
    return {V, T, F};
}

} // namespace geometry

#endif // PBD_GET_SIMPLE_BAR_MODEL_H
