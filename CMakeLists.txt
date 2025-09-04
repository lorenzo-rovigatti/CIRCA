cmake_minimum_required(VERSION 3.10)
project(CIRCA LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include_directories(extern extern/spdlog/include)

# build spdlog
add_subdirectory(extern/spdlog)
set_target_properties(spdlog PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_library(circa_lib INTERFACE)

add_executable(circa_1D src/main.cpp)
add_executable(circa_2D src/main.cpp)

target_link_libraries(circa_1D PRIVATE spdlog::spdlog circa_lib)
target_compile_definitions(circa_1D PRIVATE DIM=1)

target_link_libraries(circa_2D PRIVATE spdlog::spdlog circa_lib)
target_compile_definitions(circa_2D PRIVATE DIM=2)
