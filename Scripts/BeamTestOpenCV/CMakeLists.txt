cmake_minimum_required(VERSION 2.8.12)

project(beam)

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

add_executable(beam beam.cpp)
target_link_libraries(beam ${OpenCV_LIBS})
