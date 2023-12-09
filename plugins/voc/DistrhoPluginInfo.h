#ifndef DISTRHO_PLUGIN_INFO_H_INCLUDED
#define DISTRHO_PLUGIN_INFO_H_INCLUDED

#define DISTRHO_PLUGIN_NAME  "voc"
#define DISTRHO_PLUGIN_URI   "https://Stazed@github.com/Stazed/voc"
#define DISTRHO_PLUGIN_CLAP_ID "com.stazed.github.voc"

#define DISTRHO_PLUGIN_NUM_INPUTS   2
#define DISTRHO_PLUGIN_NUM_OUTPUTS  1
#define DISTRHO_PLUGIN_IS_RT_SAFE   1

enum Parameters {
    kBypass,
    kParameterCount
};

#endif
