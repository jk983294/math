cmake_minimum_required(VERSION 3.5)

project(ornate)

set(CMAKE_INSTALL_PREFIX "/opt/version/latest/ornate/")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/)

find_package(OpenMP REQUIRED)
find_library(FST_LIB libfst.a HINTS /opt/version/latest/cppfst/lib REQUIRED)
find_library(LZ4_LIB libliblz4.a HINTS /opt/version/latest/cppfst/lib REQUIRED)
find_library(ZSTD_LIB liblibzstd.a HINTS /opt/version/latest/cppfst/lib REQUIRED)

SET(CMAKE_CXX_FLAGS "-std=c++17 -O3 -s -Wall -Werror")
#SET(CMAKE_CXX_FLAGS "-std=c++17 -g -Wall -Werror -fsanitize=address -fno-omit-frame-pointer")

message ("cxx Flags: " ${CMAKE_CXX_FLAGS})
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)

# Source code
include_directories(SYSTEM /opt/3rd/common/include)
include_directories(include)
include_directories(SYSTEM /opt/3rd/Eigen)
include_directories(SYSTEM /opt/3rd/gcem)
include_directories(SYSTEM /opt/3rd/stats)
include_directories(/opt/version/latest/zergtool/include)
include_directories(SYSTEM /opt/version/latest/cppfst/include)

link_directories(/opt/3rd/common/lib)

install(DIRECTORY include/ DESTINATION include)

add_subdirectory(example)
add_subdirectory(test)
add_subdirectory(src)
