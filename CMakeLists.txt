cmake_minimum_required(VERSION 3.20)
project(MAXWELL)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3 -fopenmp")

add_executable(MAXWELL main.cpp)
