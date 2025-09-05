cmake_minimum_required(VERSION 3.10)
project(CIRCA LANGUAGES CXX)

option(Debug "Set to ON to compile with debug symbols" OFF)
option(G "Set to ON to compile with optimisations and debug symbols" OFF)
option(NATIVE_COMPILATION "Set to OFF to compile without the -march=native flag. This may be required when compiling binaries to be used elsewhere" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_BUILD_TYPE Release)

if(Debug)
	set(CMAKE_BUILD_TYPE Debug)
	message(STATUS "Compiling in Debug mode")
endif(Debug)

if(G)
	set(CMAKE_BUILD_TYPE RelWithDebInfo)
	message(STATUS "Compiling in Release+Debug symbols mode")
endif(G)

if(NOT CMAKE_BUILD_TYPE MATCHES Debug)
	add_definitions(-ffast-math)
	
	if(NATIVE_COMPILATION)
		include(CheckCXXCompilerFlag)
        unset(COMPILER_SUPPORTS_MARCH_NATIVE CACHE)
        CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
        if(COMPILER_SUPPORTS_MARCH_NATIVE)
                ADD_DEFINITIONS(-march=native)
		else()
			message(STATUS "The '-march=native' flag is not supported, compiling without it")
        endif()
	else()
		message(STATUS "Compiling without the '-march=native' flag")
	endif()
endif()

# build spdlog
add_subdirectory(extern/spdlog)
set_target_properties(spdlog PROPERTIES POSITION_INDEPENDENT_CODE ON)

include_directories(extern extern/spdlog/include)

set(sources
    src/util/config.cpp
	src/util/strings.cpp
)

add_library(circa_lib ${sources})
target_link_libraries(circa_lib PRIVATE spdlog::spdlog)

add_executable(circa_1D src/main.cpp)
add_executable(circa_2D src/main.cpp)
add_executable(circa_3D src/main.cpp)

target_link_libraries(circa_1D PRIVATE circa_lib)
target_compile_definitions(circa_1D PRIVATE DIM=1)

target_link_libraries(circa_2D PRIVATE circa_lib)
target_compile_definitions(circa_2D PRIVATE DIM=2)

target_link_libraries(circa_3D PRIVATE circa_lib)
target_compile_definitions(circa_3D PRIVATE DIM=3)
