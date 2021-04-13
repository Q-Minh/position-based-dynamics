cmake_minimum_required(VERSION 3.14)
project(pbd VERSION 0.0.1 LANGUAGES CXX)

set_property(GLOBAL PROPERTY USE_FOLDER ON)

include(FetchContent)

set(LIBIGL_USE_STATIC_LIBRARY     ON  CACHE STRING   "Use libigl as static library" )
set(LIBIGL_WITH_CGAL              OFF CACHE INTERNAL "Use CGAL"                     )
set(LIBIGL_WITH_COMISO            OFF CACHE INTERNAL "Use CoMiso"                   )
set(LIBIGL_WITH_CORK              OFF CACHE INTERNAL "Use Cork"                     )
set(LIBIGL_WITH_EMBREE            OFF CACHE INTERNAL "Use Embree"                   )
set(LIBIGL_WITH_MATLAB            OFF CACHE INTERNAL "Use Matlab"                   )
set(LIBIGL_WITH_MOSEK             OFF CACHE INTERNAL "Use MOSEK"                    )
set(LIBIGL_WITH_PNG               OFF CACHE INTERNAL "Use PNG"                      )
set(LIBIGL_WITH_TETGEN            ON  CACHE INTERNAL "Use Tetgen"                   )
set(LIBIGL_WITH_TRIANGLE          OFF CACHE INTERNAL "Use Triangle"                 )
set(LIBIGL_WITH_PREDICATES        OFF CACHE INTERNAL "Use exact predicates"         )
set(LIBIGL_WITH_XML               OFF CACHE INTERNAL "Use XML"                      )
set(LIBIGL_WITHOUT_COPYLEFT       OFF CACHE INTERNAL "Disable Copyleft libraries"   )
set(LIBIGL_WITH_OPENGL            ON  CACHE INTERNAL "download libigl with opengl"  )
set(LIBIGL_WITH_OPENGL_GLFW       ON  CACHE INTERNAL "download libigl with glfw"    )
set(LIBIGL_WITH_OPENGL_GLFW_IMGUI ON  CACHE INTERNAL "download libigl with imgui"   )
FetchContent_Declare(
  _libigl
  GIT_REPOSITORY https://github.com/libigl/libigl
  GIT_TAG        v2.2.0
)
FetchContent_MakeAvailable(_libigl)

FetchContent_Declare(matplotplusplus
    GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus
    GIT_TAG origin/master
)
FetchContent_GetProperties(matplotplusplus)
if(NOT matplotplusplus_POPULATED)
    FetchContent_Populate(matplotplusplus)
    add_subdirectory(${matplotplusplus_SOURCE_DIR} ${matplotplusplus_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

add_executable(pbd)
set_target_properties(pbd PROPERTIES FOLDER position-based-dynamics)
target_compile_features(pbd PRIVATE cxx_std_17)

target_include_directories(pbd
PRIVATE
    include
)

target_sources(pbd
PRIVATE
    # source files

    ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp

    # xpbd
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xpbd/deformable_mesh.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xpbd/edge_length_constraint.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xpbd/green_strain_elastic_constraint.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xpbd/neohookean_elasticity_constraint.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xpbd/tetrahedron_volume_constraint.cpp

    # solver
    ${CMAKE_CURRENT_SOURCE_DIR}/src/solver/solve.cpp

    # ui
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ui/mouse_down_handler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ui/mouse_move_handler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ui/pre_draw_handler.cpp

    # header files

    # geometry
    ${CMAKE_CURRENT_SOURCE_DIR}/include/geometry/get_simple_bar_model.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/geometry/get_simple_cloth_model.h

    # xpbd
    ${CMAKE_CURRENT_SOURCE_DIR}/include/xpbd/constraint.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/xpbd/deformable_mesh.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/xpbd/edge_length_constraint.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/xpbd/green_strain_elastic_constraint.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/xpbd/neohookean_elasticity_constraint.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/xpbd/tetrahedron_volume_constraint.h

    # solver
    ${CMAKE_CURRENT_SOURCE_DIR}/include/solver/solve.h

    # ui
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ui/mouse_down_handler.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ui/mouse_move_handler.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ui/physics_params.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ui/picking_state.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ui/pre_draw_handler.h
)

target_link_libraries(pbd 
PRIVATE 
    igl::core 
    igl::tetgen
    igl::opengl_glfw_imgui
)

add_executable(pbd-plot)
set_target_properties(pbd-plot PROPERTIES FOLDER position-based-dynamics)
target_compile_features(pbd-plot PRIVATE cxx_std_17)

target_include_directories(pbd-plot
PRIVATE
    include
)

target_sources(pbd-plot
PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/plot.cpp

    # xpbd
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xpbd/deformable_mesh.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xpbd/edge_length_constraint.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xpbd/green_strain_elastic_constraint.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xpbd/neohookean_elasticity_constraint.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/xpbd/tetrahedron_volume_constraint.cpp

    # solver
    ${CMAKE_CURRENT_SOURCE_DIR}/src/solver/solve.cpp
)

target_link_libraries(pbd-plot PRIVATE matplot igl::core igl::tetgen)
