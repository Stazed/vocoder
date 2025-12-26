
include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR}/DPF/distrho .)

add_library(Vocoder_lv2 SHARED ${CMAKE_SOURCE_DIR}/DPF/distrho/DistrhoPluginMain.cpp DSP.cpp Vocoder.cpp)
add_library(Vocoder_ladspa SHARED ${CMAKE_SOURCE_DIR}/DPF/distrho/DistrhoPluginMain.cpp DSP.cpp Vocoder.cpp)
add_library(Vocoder_clap SHARED ${CMAKE_SOURCE_DIR}/DPF/distrho/DistrhoPluginMain.cpp DSP.cpp Vocoder.cpp)
add_library(Vocoder_vst2 SHARED ${CMAKE_SOURCE_DIR}/DPF/distrho/DistrhoPluginMain.cpp DSP.cpp Vocoder.cpp)
add_library(Vocoder_vst3 SHARED ${CMAKE_SOURCE_DIR}/DPF/distrho/DistrhoPluginMain.cpp DSP.cpp Vocoder.cpp)

set_target_properties(Vocoder_lv2 PROPERTIES COMPILE_DEFINITIONS "DISTRHO_PLUGIN_TARGET_LV2")
set_target_properties(Vocoder_lv2 PROPERTIES LIBRARY_OUTPUT_DIRECTORY "lv2")
set_target_properties(Vocoder_lv2 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "lv2")
set_target_properties(Vocoder_lv2 PROPERTIES OUTPUT_NAME "Vocoder")
set_target_properties(Vocoder_lv2 PROPERTIES PREFIX "")

set_target_properties(Vocoder_ladspa PROPERTIES COMPILE_DEFINITIONS "DISTRHO_PLUGIN_TARGET_LADSPA")
set_target_properties(Vocoder_ladspa PROPERTIES LIBRARY_OUTPUT_DIRECTORY "ladspa")
set_target_properties(Vocoder_ladspa PROPERTIES RUNTIME_OUTPUT_DIRECTORY "ladspa")
set_target_properties(Vocoder_ladspa PROPERTIES OUTPUT_NAME "Vocoder")
set_target_properties(Vocoder_ladspa PROPERTIES PREFIX "")

set_target_properties(Vocoder_clap PROPERTIES COMPILE_DEFINITIONS "DISTRHO_PLUGIN_TARGET_CLAP")
set_target_properties(Vocoder_clap PROPERTIES LIBRARY_OUTPUT_DIRECTORY "clap")
set_target_properties(Vocoder_clap PROPERTIES RUNTIME_OUTPUT_DIRECTORY "clap")
set_target_properties(Vocoder_clap PROPERTIES OUTPUT_NAME "Vocoder")
set_target_properties(Vocoder_clap PROPERTIES PREFIX "")
set_target_properties(Vocoder_clap PROPERTIES SUFFIX ".clap")

set_target_properties(Vocoder_vst2 PROPERTIES COMPILE_DEFINITIONS "DISTRHO_PLUGIN_TARGET_VST2")
set_target_properties(Vocoder_vst2 PROPERTIES LIBRARY_OUTPUT_DIRECTORY "vst")
set_target_properties(Vocoder_vst2 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "vst")
set_target_properties(Vocoder_vst2 PROPERTIES OUTPUT_NAME "Vocoder")
set_target_properties(Vocoder_vst2 PROPERTIES PREFIX "")
set_target_properties(Vocoder_vst2 PROPERTIES SUFFIX ".so")

# Determine architecture automatically
if(NOT DEFINED VST3_ARCH)
    execute_process(
        COMMAND uname -m
        OUTPUT_VARIABLE VST3_ARCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
endif()

# Linux VST3 bundle paths
set(VST3_BUNDLE_FOLDER "${CMAKE_CURRENT_BINARY_DIR}/vst3/Vocoder.vst3")
set(VST3_CONTENTS_ARCH_FOLDER "${VST3_BUNDLE_FOLDER}/Contents/${VST3_ARCH}-linux")

set_target_properties(Vocoder_vst3 PROPERTIES COMPILE_DEFINITIONS "DISTRHO_PLUGIN_TARGET_VST3")
set_target_properties(Vocoder_vst3 PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${VST3_CONTENTS_ARCH_FOLDER}")
set_target_properties(Vocoder_vst3 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "vst3")
set_target_properties(Vocoder_vst3 PROPERTIES OUTPUT_NAME "Vocoder")
set_target_properties(Vocoder_vst3 PROPERTIES PREFIX "")
set_target_properties(Vocoder_vst3 PROPERTIES SUFFIX ".so")

target_link_libraries (Vocoder_lv2 ${ExternLibraries})
target_link_libraries (Vocoder_ladspa ${ExternLibraries})
target_link_libraries (Vocoder_clap ${ExternLibraries})
target_link_libraries (Vocoder_vst2 ${ExternLibraries})
target_link_libraries (Vocoder_vst3 ${ExternLibraries})


install(TARGETS Vocoder_lv2 LIBRARY DESTINATION ${PluginLibDir}/lv2/Vocoder.lv2/)
install(TARGETS Vocoder_ladspa LIBRARY DESTINATION ${PluginLibDir}/ladspa/)
install(TARGETS Vocoder_clap LIBRARY DESTINATION ${PluginLibDir}/clap/)
install(TARGETS Vocoder_vst2 LIBRARY DESTINATION ${PluginLibDir}/vst/)

add_custom_command(TARGET Vocoder_lv2 POST_BUILD
    COMMAND ../../../lv2-ttl-generator $<TARGET_FILE:Vocoder_lv2>
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lv2)

add_dependencies(Vocoder_lv2 lv2-ttl-generator)

install(FILES
	${CMAKE_CURRENT_BINARY_DIR}/lv2/manifest.ttl
	${CMAKE_CURRENT_BINARY_DIR}/lv2/Vocoder.ttl
	DESTINATION ${PluginLibDir}/lv2/Vocoder.lv2/)


# Linux VST3 bundle structure
set(VST3_BUNDLE_FOLDER "${CMAKE_CURRENT_BINARY_DIR}/vst3/Vocoder.vst3")
set(VST3_CONTENTS_ARCH_FOLDER "${VST3_BUNDLE_FOLDER}/Contents/${VST3_ARCH}-linux")

# Post-build: create directories and copy the .so
add_custom_command(TARGET Vocoder_vst3 POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory "${VST3_CONTENTS_ARCH_FOLDER}"
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:Vocoder_vst3> "${VST3_CONTENTS_ARCH_FOLDER}/Vocoder.so"
)

install(DIRECTORY "${VST3_BUNDLE_FOLDER}/"
        DESTINATION "${PluginLibDir}/vst3/Vocoder.vst3"
)
