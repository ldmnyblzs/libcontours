cmake_minimum_required(VERSION 3.9)
project(SHAPE)
set (CMAKE_CXX_STANDARD 17)

find_package(CGAL REQUIRED)
find_package(Boost REQUIRED)

set(bliss_dir src/bliss-0.73)
add_library(bliss STATIC
  ${bliss_dir}/defs.cc
  ${bliss_dir}/graph.cc
  ${bliss_dir}/partition.cc
  ${bliss_dir}/orbit.cc
  ${bliss_dir}/uintseqhash.cc
  ${bliss_dir}/heap.cc
  ${bliss_dir}/timer.cc
  ${bliss_dir}/utils.cc
  ${bliss_dir}/bliss_C.cc)

add_library(shape INTERFACE)
target_include_directories(shape INTERFACE include)
target_include_directories(shape INTERFACE src)
target_link_libraries(shape INTERFACE CGAL::CGAL Boost::boost bliss)

add_subdirectory(example)
