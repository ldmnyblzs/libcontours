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
target_include_directories(bliss PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include/bliss>)

add_library(shape INTERFACE)
target_include_directories(shape INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/shape>)
target_link_libraries(shape INTERFACE CGAL::CGAL Boost::boost bliss)

install(DIRECTORY src/bliss-0.73 DESTINATION include)
install(TARGETS bliss EXPORT bliss_target DESTINATION lib)
install(EXPORT bliss_target DESTINATION lib/cmake/bliss-0.73)
  
install(DIRECTORY include/shape DESTINATION include)
install(TARGETS shape EXPORT shape_target DESTINATION lib)
install(EXPORT shape_target DESTINATION lib/cmake/shape)

add_subdirectory(example)
