cmake_minimum_required (VERSION 3.5)
include(FindPackageHandleStandardArgs)

project (nfm LANGUAGES CXX VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${USER_CXX_FLAGS}")

# find packages
message(STATUS "Configured CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")
message(STATUS "Configured CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

# set header / library paths
include_directories(include/) # headers

enable_testing()

# continue with subdirectories
add_subdirectory(src) 
add_subdirectory(test)
add_subdirectory(examples)
