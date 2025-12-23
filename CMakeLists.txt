#CMake file for Vocoder LV2/CLAP

cmake_minimum_required(VERSION 3.5.1...3.31.1)

project (Vocoder)

set (PACKAGE_VERSION "1.0.0")

option (BuildForDebug "Include gdb debugging support" OFF)
option (UseKISSfft "Use the KISSfft library instead of fftw" OFF)
option (UseCarrierBypass "Use the carrier as pass through on bypass (default voice)" OFF)

SET (PluginLibDir "lib" CACHE STRING
    "Install directory for plugin libraries PREFIX/PLUGIN_LIB_DIR/{lv2,clap}")

find_package(Threads REQUIRED)
if(THREADS_HAVE_PTHREAD_ARG)
  set_property(TARGET Vocoder PROPERTY COMPILE_OPTIONS "-pthread")
  set_property(TARGET Vocoder PROPERTY INTERFACE_COMPILE_OPTIONS "-pthread")
endif()

set(CMAKE_BUILD_TYPE "Release")

set (BuildOptionsBasic
    "-Wno-unused-parameter -O3 -fvisibility=hidden"
    CACHE STRING "basic X86 complier options"
)

set (BuildOptionsDebug
    "-std=c++11 -Wno-cast-function-type -O0 -g3 -ggdb -Wall -Wextra -Wpointer-arith"
    CACHE STRING "Debug build flags"
)

if (UseCarrierBypass)
    message(STATUS "Using carrier on plugin bypass")
    add_definitions(-D'USE_CARRIER_ON_BYPASS=1')
endif (UseCarrierBypass)

if (UseKISSfft)
    find_package(kissfft CONFIG REQUIRED COMPONENTS STATIC)
    if(TARGET kissfft::kissfft)   # Modern cmake
        message(STATUS "Have kissfft::kissfft TARGET")
    endif()
    message(STATUS "Using KISSfft library")
    add_definitions(-D'KISSFFT_SUPPORT=1')
    set(ExternLibraries
        libkissfft-float.a
    )
else (UseKISSfft)
    set(ExternLibraries
        fftw3
    )
endif (UseKISSfft)

if (BuildForDebug)
    set (CMAKE_BUILD_TYPE "Debug")
    set (CMAKE_CXX_FLAGS_DEBUG ${BuildOptionsDebug})
else(BuildForDebug)
    set(CMAKE_BUILD_TYPE "Release")
    set (CMAKE_CXX_FLAGS_RELEASE ${BuildOptionsBasic})
endif(BuildForDebug)

add_subdirectory(plugins/Vocoder)

add_executable(lv2-ttl-generator ${CMAKE_SOURCE_DIR}/DPF/utils/lv2-ttl-generator/lv2_ttl_generator.c)
target_link_libraries(lv2-ttl-generator dl)

# uninstall target
if(NOT TARGET uninstall)
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()

