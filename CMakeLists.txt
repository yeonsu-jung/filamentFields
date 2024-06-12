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
#set(pybind11_DIR /Library/Frameworks/Python.framework/Versions/3.11/lib/python3.11/site-packages/pybind11/include)
#find_package(pybind11 REQUIRED)

set(Eigen3_DIR /opt/homebrew/Cellar/eigen/3.4.0_1/share/eigen3/cmake)
find_package(Eigen3 3.4 REQUIRED)

set(TBB_DIR /opt/homebrew/Cellar/tbb/2021.12.0/share/tbb/cmake)
find_package(TBB REQUIRED)

include_directories(${EIGEN_INCLUDE_DIRS})
#include_directories(${pybind11_INCLUDE_DIRS})

#add_executable(filamentFieldsTest filamentFieldsTest.cpp filamentFields.cpp)
#add_executable(filamentFields filamentFields.cpp filamentFieldsTest.cpp)
#add_executable(filamentFields_entanglement_only filamentFields.cpp filamentFields_entanglement_only.cpp)
#add_executable(entanglement_ca filamentFields.cpp entanglement_ca.cpp)
add_executable(filamentFields_volume_fraction_only filamentFields.cpp filamentFields_volume_fraction_only.cpp)


# Create the pybind11 module
#pybind11_add_module(filamentFields MODULE filamentFields.cpp wrapper.cpp)
#target_link_libraries(filamentFields PRIVATE Eigen3::Eigen)

#target_link_libraries(filamentFieldsTest PRIVATE Eigen3::Eigen)
#target_link_libraries(filamentFields PRIVATE Eigen3::Eigen)
#target_link_libraries(filamentFields_entanglement_only PRIVATE Eigen3::Eigen TBB::tbb)
#target_link_libraries(entanglement_ca PRIVATE Eigen3::Eigen)
#target_link_libraries(entanglement_ca PRIVATE Eigen3::Eigen TBB::tbb)
target_link_libraries(filamentFields_volume_fraction_only PRIVATE Eigen3::Eigen TBB::tbb)





