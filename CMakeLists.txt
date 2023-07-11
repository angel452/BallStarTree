cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

project(BallStarTree)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")  # Línea de optimización

include_directories($(CMAKE_CURRENT_SOURCE_DIR)/src)

find_package(Eigen3)

add_executable(BallStarTree src/main.cpp )

target_link_libraries(BallStarTree Eigen3::Eigen)

