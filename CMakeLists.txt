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
# Allow overriding pybind11 location via cache or env but do not hard-code path
if(NOT pybind11_DIR AND DEFINED ENV{PYBIND11_DIR})
	set(pybind11_DIR $ENV{PYBIND11_DIR})
endif()
if(NOT pybind11_DIR AND DEFINED ENV{Pybind11_DIR})
	set(pybind11_DIR $ENV{Pybind11_DIR})
endif()
find_package(pybind11 CONFIG REQUIRED)

# Respect user/env-provided Eigen3 location if available; otherwise rely on default FindEigen3
if(NOT Eigen3_DIR AND DEFINED ENV{Eigen3_DIR})
	set(Eigen3_DIR $ENV{Eigen3_DIR})
endif()
if(NOT Eigen3_DIR AND DEFINED ENV{EIGEN3_DIR})
	set(Eigen3_DIR $ENV{EIGEN3_DIR})
endif()
find_package(Eigen3 3.3 REQUIRED)

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
# Ensure the Python module outputs to the source dir for easy import
set_target_properties(filamentFields PROPERTIES
	LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	OUTPUT_NAME filamentFields
)

target_link_libraries(filamentFieldsTest PRIVATE Eigen3::Eigen TBB::tbb)

# Optional CUDA support and smoke test
option(FF_WITH_CUDA "Enable CUDA acceleration and build CUDA smoke test" OFF)
if(FF_WITH_CUDA)
	message(STATUS "FF_WITH_CUDA=ON: Enabling CUDA language")
	enable_language(CUDA)

	# Discover CUDA toolkit imported targets (CUDA::cudart, etc.)
	find_package(CUDAToolkit REQUIRED)

	# Target a reasonable set of architectures commonly present on FAS RC (V100/A100/H100)
	# Users can override by setting CMAKE_CUDA_ARCHITECTURES at configure time
	if(NOT CMAKE_CUDA_ARCHITECTURES)
		set(CMAKE_CUDA_ARCHITECTURES 70;80;90)
	endif()

		# Add CUDA sources (prototype matrix kernel and streaming kernel)
		set(CUDA_SOURCES
				compute_linking_matrix.cu
				streaming_gpu.cu)
	add_library(filamentFieldsCUDA STATIC ${CUDA_SOURCES})
			# Ensure PIC for linking into the shared pybind11 module
			set_target_properties(filamentFieldsCUDA PROPERTIES POSITION_INDEPENDENT_CODE ON)
			# Extra safety: pass -fPIC to host compiler via nvcc for CUDA sources
			target_compile_options(filamentFieldsCUDA PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-fPIC>)

	# Link CUDA lib to module if available (no-op currently, safe to link)
				target_link_libraries(filamentFields PRIVATE filamentFieldsCUDA CUDA::cudart)

		# Define a compile-time macro so C++ can enable the GPU pathway
		target_compile_definitions(filamentFields PRIVATE FILAMENTFIELDS_WITH_CUDA=1)

	# Simple CUDA smoke test executable
	add_executable(ff_cuda_smoke tests/cuda_smoke.cu)
	# No external deps needed for the smoke test
endif()