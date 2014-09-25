project(stvk-model)
cmake_minimum_required(VERSION 2.8)


set(CMAKE_CXX_FLAGS "-fopenmp -fpermissive -std=c++0x ${CMAKE_CXX_FLAGS}")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/;${CMAKE_MODULE_PATH}")
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib/)

include_directories(${CMAKE_SOURCE_DIR})
link_directories(${LIBRARY_OUTPUT_PATH})

include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
include_geo_sim_sdk()
link_geo_sim_sdk()

#EIGEN3
find_package(Eigen3 REQUIRED)
if(EIGEN3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIR})
endif(EIGEN3_FOUND)

#BLAS and LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_subdirectory(src)
add_subdirectory(example)

