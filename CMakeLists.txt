cmake_minimum_required(VERSION 3.9)
project(main)

set(CMAKE_CXX_STANDARD 11)

add_subdirectory(lib)
add_executable(main main.cpp)
target_link_libraries(main PUBLIC transformations)
