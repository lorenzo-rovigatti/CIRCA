cmake_minimum_required(VERSION 3.16)
project(CIRCA LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_library(circa_lib INTERFACE)
target_include_directories(circa_lib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(circa_1D src/main.cpp)
add_executable(circa_2D src/main.cpp)

target_link_libraries(circa_1D PRIVATE circa_lib)
target_compile_definitions(circa_1D PRIVATE DIM=1)

target_link_libraries(circa_2D PRIVATE circa_lib)
target_compile_definitions(circa_2D PRIVATE DIM=2)
