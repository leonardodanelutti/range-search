cmake_minimum_required(VERSION 3.10)

project(RangeSearchAlgorithms)

set(CMAKE_CXX_STANDARD 17)

# Default to Release build if no build type is specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

include_directories(include)

find_package(CGAL REQUIRED COMPONENTS Core OPTIONAL_COMPONENTS Qt5)
if(CGAL_FOUND)
	message(STATUS "CGAL found: ${CGAL_VERSION}")
	include(${CGAL_USE_FILE})
else()
	message(FATAL_ERROR "CGAL not found. Please install CGAL.")
endif()

if(CGAL_Qt5_FOUND)
	message(STATUS "CGAL Qt5 found")
	find_package(Qt5 REQUIRED COMPONENTS Core Widgets)
	set(CMAKE_AUTOMOC ON)
else()
	message(STATUS "CGAL Qt5 not found - GUI features will be disabled")
endif()

find_package(Eigen3 REQUIRED)
if(Eigen3_FOUND)
	message(STATUS "Eigen3 found: ${EIGEN3_VERSION}")
	include_directories(${EIGEN3_INCLUDE_DIR})
else()
	message(FATAL_ERROR "Eigen3 not found. Please install Eigen3.")
endif()

# Create visualization executable (requires Qt5)
if(CGAL_Qt5_FOUND)
	add_executable(range_search_viz src/visualization_main.cpp src/visualization.cpp)
	target_link_libraries(range_search_viz CGAL::CGAL CGAL::CGAL_Core CGAL::CGAL_Qt5 Qt5::Core Qt5::Widgets)
	target_include_directories(range_search_viz PRIVATE ${EIGEN3_INCLUDE_DIR})
	message(STATUS "Building visualization executable: range_search_viz")
else()
	message(STATUS "Skipping visualization executable - Qt5 not found")
endif()

# Create benchmark executable (console only, no Qt5 required)
add_executable(range_search_benchmark src/benchmark_main.cpp)
target_link_libraries(range_search_benchmark CGAL::CGAL CGAL::CGAL_Core)
target_include_directories(range_search_benchmark PRIVATE ${EIGEN3_INCLUDE_DIR})
message(STATUS "Building benchmark executable: range_search_benchmark")
