cmake_minimum_required(VERSION 3.10)

project(RangeSearchAlgorithms)

set(CMAKE_CXX_STANDARD 17)

include_directories(include)

find_package(CGAL REQUIRED COMPONENTS Core)
if(CGAL_FOUND)
	message(STATUS "CGAL found: ${CGAL_VERSION}")
	include(${CGAL_USE_FILE})
else()
	message(FATAL_ERROR "CGAL not found. Please install CGAL.")
endif()

find_package(Eigen3 REQUIRED)
if(Eigen3_FOUND)
	message(STATUS "Eigen3 found: ${EIGEN3_VERSION}")
	include_directories(${EIGEN3_INCLUDE_DIR})
else()
	message(FATAL_ERROR "Eigen3 not found. Please install Eigen3.")
endif()

file(GLOB SOURCES "src/*.cpp")
add_executable(range_search ${SOURCES})
target_link_libraries(range_search CGAL::CGAL CGAL::CGAL_Core)
target_include_directories(range_search PRIVATE ${EIGEN3_INCLUDE_DIR})

enable_testing()
add_subdirectory(tests)
