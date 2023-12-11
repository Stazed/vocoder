#CMake file for Vocoder LV2/CLAP

cmake_minimum_required(VERSION 3.5.1)

project (Vocoder)

set (PACKAGE_VERSION "1.0.0")

option (BuildForDebug "Include gdb debugging support" OFF)

find_package(Threads REQUIRED)
if(THREADS_HAVE_PTHREAD_ARG)
  set_property(TARGET Vocoder PROPERTY COMPILE_OPTIONS "-pthread")
  set_property(TARGET Vocoder PROPERTY INTERFACE_COMPILE_OPTIONS "-pthread")
endif()

set(ExternLibraries
    fftw3
)

set (BuildOptionsDebug
    "-std=c++11 -Wno-cast-function-type -O0 -g3 -ggdb -Wall -Wextra -Wpointer-arith"
    CACHE STRING "Debug build flags"
)

if (BuildForDebug)
    set (CMAKE_BUILD_TYPE "Debug")
    set (CMAKE_CXX_FLAGS_DEBUG ${BuildOptionsDebug})
endif(BuildForDebug)

add_subdirectory(plugins/Vocoder)

add_executable(lv2-ttl-generator ${CMAKE_SOURCE_DIR}/DPF/utils/lv2-ttl-generator/lv2_ttl_generator.c)

