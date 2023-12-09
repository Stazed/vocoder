
include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR}/DPF/distrho .)

add_library(voc_lv2 SHARED ${CMAKE_SOURCE_DIR}/DPF/distrho/DistrhoPluginMain.cpp DSP.cpp voc.cpp)
add_library(voc_clap SHARED ${CMAKE_SOURCE_DIR}/DPF/distrho/DistrhoPluginMain.cpp DSP.cpp voc.cpp)

    
set_target_properties(voc_lv2 PROPERTIES COMPILE_DEFINITIONS "DISTRHO_PLUGIN_TARGET_LV2")
set_target_properties(voc_lv2 PROPERTIES LIBRARY_OUTPUT_DIRECTORY "lv2")
set_target_properties(voc_lv2 PROPERTIES RUNTIME_OUTPUT_DIRECTORY "lv2")
set_target_properties(voc_lv2 PROPERTIES OUTPUT_NAME "voc")
set_target_properties(voc_lv2 PROPERTIES PREFIX "")

set_target_properties(voc_clap PROPERTIES COMPILE_DEFINITIONS "DISTRHO_PLUGIN_TARGET_CLAP")
set_target_properties(voc_clap PROPERTIES LIBRARY_OUTPUT_DIRECTORY "clap")
set_target_properties(voc_clap PROPERTIES RUNTIME_OUTPUT_DIRECTORY "clap")
set_target_properties(voc_clap PROPERTIES OUTPUT_NAME "voc")
set_target_properties(voc_clap PROPERTIES PREFIX "")

target_link_libraries (voc_lv2 ${ExternLibraries})


install(TARGETS voc_lv2 LIBRARY DESTINATION ${PluginLibDir}/lv2/voc.lv2/)
install(TARGETS voc_clap LIBRARY DESTINATION ${PluginLibDir}/clap/)

add_custom_command(TARGET voc_lv2 POST_BUILD
    COMMAND ../../../lv2-ttl-generator $<TARGET_FILE:voc_lv2>
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lv2)

add_dependencies(voc_lv2 lv2-ttl-generator)

install(FILES
	${CMAKE_CURRENT_BINARY_DIR}/lv2/manifest.ttl
	${CMAKE_CURRENT_BINARY_DIR}/lv2/presets.ttl
	${CMAKE_CURRENT_BINARY_DIR}/lv2/voc.ttl
	DESTINATION ${PluginLibDir}/lv2/voc.lv2/)