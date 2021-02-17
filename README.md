# Soft Body Virtual Cutting using Position Based Dynamics

## Overview

Academic prototyping project for soft body cutting using (X)PBD. 
Different constraint types and cutting methods will be implemented 
using [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) for matrix computations and [libigl](https://libigl.github.io/) for visualization and 
user interaction.

### Short term roadmap

#### Constraint types
- edge length
- tetrahedral volume
- fem based
- shape matching

#### Cutting methods
- element deletion
- swept surface

![Cloth simulation example](./doc/pbd-simple-cloth-example.gif)

## Dependencies

- C++17 compiler
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
- [libigl](https://libigl.github.io/)

[libigl](https://libigl.github.io/) is included in the projet using CMake's FetchContent and pulls in [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page), [glfw](https://www.glfw.org/) and [Dear ImGui](https://github.com/ocornut/imgui) with it.

## Building

```
# Download repository
$ git clone https://github.com/Q-Minh/position-based-dynamics
$ cd position-based-dynamics

# Configure and build projet
$ cmake -S . -B build
$ cmake --build build --target pbd --config Release

# Run the program
$ ./build/Release/pbd.exe
```
