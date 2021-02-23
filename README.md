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
![Soft body simulation example](./doc/bunny-pbd-edge-length.gif)

## Dependencies

- C++17 compiler
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
- [libigl](https://libigl.github.io/)
- [boost](https://www.boost.org/)

[libigl](https://libigl.github.io/) is included in the projet using CMake's FetchContent and pulls in [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page), [glfw](https://www.glfw.org/) and [Dear ImGui](https://github.com/ocornut/imgui) with it.

The [boost libraries](https://www.boost.org/doc/libs/) are required by [CGAL](https://www.cgal.org/) and are not provided automatically by [libigl](https://libigl.github.io/). There are multiple ways to download the boost libraries. Use the download links [found here](https://www.boost.org/users/download/) or use git to download [modular boost](https://github.com/boostorg/wiki/wiki/Getting-Started%3A-Overview). Use the [git-submodule interface](https://git-scm.com/docs/git-submodule) to download only the [boost libraries](https://www.boost.org/doc/libs/) required by [CGAL](https://www.cgal.org/).

Boost full installation example:
```
$ git clone --recursive https://github.com/boostorg/boost.git
$ cd boost
$ git checkout develop # or whatever branch you want to use
$ ./bootstrap.sh
$ ./b2 headers
```

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
