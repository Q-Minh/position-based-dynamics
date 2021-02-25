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

# find_package(range-v3 CONFIG REQUIRED)

add_executable(pbd)
set_target_properties(pbd PROPERTIES FOLDER position-based-dynamics)
target_compile_features(pbd PRIVATE cxx_std_17)

target_include_directories(pbd
PRIVATE
    include
)

target_sources(pbd
PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/src/pbd/deformable_mesh.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/pbd/edge_length_constraint.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/pbd/green_strain_elastic_constraint.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/pbd/solve.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/pbd/tetrahedron_volume_constraint.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/src/ui/mouse_down_handler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ui/mouse_move_handler.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ui/pre_draw_handler.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/include/pbd/constraint.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/pbd/deformable_mesh.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/pbd/edge_length_constraint.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/pbd/get_simple_cloth_model.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/pbd/green_strain_elastic_constraint.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/pbd/solve.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/pbd/tetrahedron_volume_constraint.h

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
    # range-v3::range-v3
)
