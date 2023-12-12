/*
  Vocoder

  DistrhoPluginInfo.h - DPF information header

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#ifndef DISTRHO_PLUGIN_INFO_H_INCLUDED
#define DISTRHO_PLUGIN_INFO_H_INCLUDED

#define DISTRHO_PLUGIN_NAME  "Vocoder"
#define DISTRHO_PLUGIN_URI   "https://github.com/Stazed/vocoder.git"
#define DISTRHO_PLUGIN_CLAP_ID "com.stazed.github.vocoder"

#define DISTRHO_PLUGIN_NUM_INPUTS   2
#define DISTRHO_PLUGIN_NUM_OUTPUTS  1
#define DISTRHO_PLUGIN_IS_RT_SAFE   1
#define DISTRHO_PLUGIN_LV2_CATEGORY  "lv2:ModulatorPlugin"

enum Parameters {
    kBypass,
    kParameterCount
};

#endif
