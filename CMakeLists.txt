cmake_minimum_required(VERSION 3.12)
project(filamentFields)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-O3 -DNDEBUG -Wno-deprecated-declarations")
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

# Set build type and flags
# set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

# Set output directory
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

# Set threading preference
set(THREADS_PREFER_PTHREAD_FLAG ON)

# Find necessary packages
# set(pybind11_DIR /home1/99999/yjung3/.local/lib/python3.9/site-packages/pybind11/share/cmake/pybind11)
set(pybind11_DIR /Users/yeonsu/anaconda3/envs/jax-env/lib/python3.8/site-packages/pybind11/share/cmake/pybind11)
find_package(pybind11 REQUIRED)

set(Eigen3_DIR /opt/homebrew/Cellar/eigen/3.4.0_1/share/eigen3/cmake)
find_package(Eigen3 3.4 REQUIRED)

find_package(TBB REQUIRED)


# find_package(Python3 REQUIRED COMPONENTS Development)
# include_directories(${Python3_INCLUDE_DIRS})

#/scratch/projects/compilers/intel24.0/oneapi/intelpython/python3.9/include/python3.9
#scratch/projects/compilers/intel24.0/oneapi/intelpython/python3.9/lib/python3.9


# include_directories(/home1/99999/yjung3/myenv/include/python3.9)
include_directories(${EIGEN_INCLUDE_DIRS})
include_directories(${pybind11_INCLUDE_DIRS})

add_executable(filamentFieldsTest filamentFieldsTest.cpp filamentFields.cpp)

# Create the pybind11 module
pybind11_add_module(filamentFields MODULE filamentFields.cpp wrapper.cpp)
target_link_libraries(filamentFields PRIVATE Eigen3::Eigen TBB::tbb)

target_link_libraries(filamentFieldsTest PRIVATE Eigen3::Eigen TBB::tbb)