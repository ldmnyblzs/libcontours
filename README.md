# libcontours
This library implements the algorithm described in https://doi.org/10.3311/PPee.12313

## Dependencies
Boost, CGAL and Bliss are required, the latter of which can be found here: https://github.com/ldmnyblzs/bliss-cmake

## Compilation
```
mkdir build
cd build
bliss_DIR=../../bliss-cmake/build/lib/cmake/bliss/ cmake .. -DCMAKE_INSTALL_PREFIX=.
cmake --build .
cmake --build . --target install
```
