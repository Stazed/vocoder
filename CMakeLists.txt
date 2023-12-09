#CMake file for voc LV2/CLAP

cmake_minimum_required(VERSION 3.5.1)

project (voc)

set (PACKAGE_VERSION "1.0.0")

set(ExternLibraries
    fftw3
)

add_subdirectory(plugins/voc)

add_executable(lv2-ttl-generator ${CMAKE_SOURCE_DIR}/DPF/utils/lv2-ttl-generator/lv2_ttl_generator.c)

