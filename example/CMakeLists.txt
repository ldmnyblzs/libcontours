find_package(Boost REQUIRED iostreams)

add_executable(mesh mesh.cpp)
target_link_libraries(mesh PUBLIC Boost::iostreams contours)
