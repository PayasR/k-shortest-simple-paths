cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
project(k_shortest_paths CXX)

include(CMakeDetermineCXXCompiler)

# default to RelWithDebInfo
if("${CMAKE_BUILD_TYPE}" STREQUAL "")
        set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(EXTRA_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -pedantic" CACHE STRING "Extra flags used by the compiler during all build types.")
        set(EXTRA_EXE_LINKER_FLAGS "" CACHE STRING "Extra flags used by the linker.")
        set(EXTRA_EXE_LINKER_FLAGS_RELEASE "-flto" CACHE STRING "Extra flags used by the linker for the Release build type.")
        set(EXTRA_EXE_LINKER_FLAGS_RELWITHDEBINFO "-flto" CACHE STRING "Extra flags used by the linker for the RelWithDebInfo build type.")
else()
        set(EXTRA_CXX_FLAGS "" CACHE STRING "Extra flags used by the compiler during all build types.")
        set(EXTRA_EXE_LINKER_FLAGS "" CACHE STRING "Extra flags used by the linker.")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_CXX_FLAGS} -std=c++14")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${EXTRA_EXE_LINKER_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${EXTRA_EXE_LINKER_FLAGS_RELEASE}")
set(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO} ${EXTRA_EXE_LINKER_FLAGS_RELWITHDEBINFO}")

option(VERBOSE "Verbose logging" OFF)

if(NOT VERBOSE)
        add_definitions(-DNVERBOSE)
endif()

# add_library(common OBJECT
#         src/algorithms.cpp
#         src/basic_types.cpp
#         src/digraph.cpp
#         src/graph.cpp
#         src/hyperbolicity.cpp
#         src/io.cpp
#         src/naive_iterator.cpp
#         src/improved_iterator1.cpp
# )

add_executable(main
        src/main.cpp
#        $<TARGET_OBJECTS:common>
)

add_executable(experiments
        src/experiments.cpp
#         $<TARGET_OBJECTS:common>
)
