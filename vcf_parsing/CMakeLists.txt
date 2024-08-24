cmake_minimum_required(VERSION 3.20.0)

project(vcf_parsing)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose Release or Debug" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(BUILD_TYPE_FLAGS "-O3")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(BUILD_TYPE_FLAGS "-g")
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(BUILD_TYPE_FLAGS "-O3 -g -DNDEBUG")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BUILD_TYPE_FLAGS}")

link_directories(/apps/samtools/1.12/lib)
include_directories(/apps/samtools/1.12/include)
include_directories(/apps/samtools/1.12/include/htslib)

include_directories(${CMAKE_SOURCE_DIR}/thirdparty)
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_BINARY_DIR}/generated)

set(VCF_LIB_DEPS ${HTSlib_LIBRARIES} dl pthread ssl crypto hts)

add_executable(read_freq main.cpp)
target_link_libraries(read_freq ${VCF_LIB_DEPS})
