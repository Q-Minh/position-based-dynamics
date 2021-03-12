# Soft Body Virtual Cutting using eXtended Position Based Dynamics

## Overview

Academic prototyping project for soft body cutting using XPBD.
Different constraint types and cutting methods will be implemented 
using [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) for matrix computations and [libigl](https://libigl.github.io/) for visualization and 
user interaction.

### Short term roadmap

#### Constraint types
- edge length
- tetrahedral volume
- fem based

#### Cutting methods
- element deletion
- swept surface

![Bar simulation example](./doc/bar-fem-pbd.gif)
![Cloth simulation example](./doc/pbd-simple-cloth-example.gif)
![Bunny geometric simulation example](./doc/bunny-pbd-edge-length.gif)

## Dependencies

- C++17 compiler
- [libigl](https://libigl.github.io/)

[libigl](https://libigl.github.io/) is included in the project using CMake's FetchContent and pulls in [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page), [glfw](https://www.glfw.org/), [Dear ImGui](https://github.com/ocornut/imgui) and [TetGen](http://wias-berlin.de/software/index.jsp?id=TetGen&lang=1) with it.

## Building

```
# Download repository
$ git clone https://github.com/Q-Minh/position-based-dynamics
$ cd position-based-dynamics

# Configure and build project
$ cmake -S . -B build
$ cmake --build build --target pbd --config Release

# Run the program
$ ./build/Release/pbd.exe
```
