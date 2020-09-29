add_library(transformations STATIC transformations.cpp radian_degree.cpp)
target_include_directories(transformations PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
