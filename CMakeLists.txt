cmake_minimum_required(VERSION 3.9)
project(SHAPE)
set (CMAKE_CXX_STANDARD 17)

option(BUILD_EXAMPLES "Build examples" ON)

find_package(CGAL REQUIRED)
find_package(Boost REQUIRED)
find_package(bliss)

include(GNUInstallDirs)
add_library(shape INTERFACE)
target_include_directories(shape INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/shape>)
target_link_libraries(shape INTERFACE CGAL::CGAL Boost::boost bliss-0.73)
  
install(DIRECTORY include/shape DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS shape
	EXPORT shapeConfig
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT shapeConfig DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/shape)
export(EXPORT shapeConfig)

if(BUILD_EXAMPLES)
  add_subdirectory(example)
endif(BUILD_EXAMPLES)
