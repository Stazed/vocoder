
include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR}/DPF/distrho .)

add_library(Vocoder_lv2 SHARED ${CMAKE_SOURCE_DIR}/DPF/distrho/DistrhoPluginMain.cpp DSP.cpp Vocoder.cpp)
add_library(Vocoder_ladspa SHARED ${CMAKE_SOURCE_DIR}/DPF/distrho/DistrhoPluginMain.cpp DSP.cpp Vocoder.cpp)
add_library(Vocoder_clap SHARED ${CMAKE_SOURCE_DIR}/DPF/distrho/DistrhoPluginMain.cpp DSP.cpp Vocoder.cpp)

    
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

target_link_libraries (Vocoder_lv2 ${ExternLibraries})
target_link_libraries (Vocoder_ladspa ${ExternLibraries})
target_link_libraries (Vocoder_clap ${ExternLibraries})


install(TARGETS Vocoder_lv2 LIBRARY DESTINATION ${PluginLibDir}/lv2/Vocoder.lv2/)
install(TARGETS Vocoder_ladspa LIBRARY DESTINATION ${PluginLibDir}/ladspa/)
install(TARGETS Vocoder_clap LIBRARY DESTINATION ${PluginLibDir}/clap/)

add_custom_command(TARGET Vocoder_lv2 POST_BUILD
    COMMAND ../../../lv2-ttl-generator $<TARGET_FILE:Vocoder_lv2>
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lv2)

add_dependencies(Vocoder_lv2 lv2-ttl-generator)

install(FILES
	${CMAKE_CURRENT_BINARY_DIR}/lv2/manifest.ttl
	${CMAKE_CURRENT_BINARY_DIR}/lv2/Vocoder.ttl
	DESTINATION ${PluginLibDir}/lv2/Vocoder.lv2/)
