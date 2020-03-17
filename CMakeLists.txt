cmake_minimum_required(VERSION 3.9)
project(CONTOURS)
set (CMAKE_CXX_STANDARD 17)

find_package(CGAL REQUIRED)
find_package(Boost REQUIRED)
find_package(bliss REQUIRED)

include(GNUInstallDirs)
add_library(contours INTERFACE)
target_include_directories(contours INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/shape>)
target_link_libraries(contours INTERFACE CGAL::CGAL Boost::boost bliss-0.73)
  
install(DIRECTORY include/contours DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS contours
	EXPORT contoursConfig
	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT contoursConfig DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/contours)
export(EXPORT contoursConfig)

add_subdirectory(example)
