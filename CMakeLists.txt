cmake_minimum_required(VERSION 3.5)

project(ornate)

set(CMAKE_INSTALL_PREFIX "/opt/version/latest/")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)


SET(CMAKE_CXX_FLAGS "-std=c++14 -O3 -Wall -Werror")
message ("cxx Flags: " ${CMAKE_CXX_FLAGS})
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)

# Source code
include_directories(SYSTEM /opt/3rd/common/include)
include_directories(base/include)
include_directories(tool/include)
include_directories(include)
include_directories(SYSTEM /opt/3rd/Eigen)

link_directories(/opt/3rd/common/lib)

install(DIRECTORY include/ DESTINATION ornate/include)

add_subdirectory(example)
#add_subdirectory(test)
