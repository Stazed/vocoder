
Vocoder
=======

Vocoder is a simple LADSPA/LV2/CLAP/VST2/VST3 plugin for vocoding based on VocProc by Igor Brkic <igor@hyperglitch.com>.
The plugin is built using the DPF framework and you must download the DPF submodule with:

```bash
    git submodule update --init
```

Building
========

Dependencies:
  -  pthread
  -  fftw3

For cmake build:

```bash
    mkdir build
    cd build
    cmake ..
    make
    sudo make install
```
To uninstall:

```bash
    sudo make uninstall
```

By default only the LADSPA/LV2/CLAP plugins are built. You can enable/disable the plugin types using cmake as follows:

Example to build all types:
```bash
   cmake -DBuildLADSPA=ON -DBuildLV2=ON -DBuildCLAP=ON -DBuildVST2=ON -DBuildVST3=ON ..
```
To enable:

```bash
    -DBuildLADSPA=ON
    -DBuildLV2=ON
    -DBuildCLAP=ON
    -DBuildVST2=ON
    -DBuildVST3=ON
```

To disable:

```bash
    -DBuildLADSPA=OFF
    -DBuildLV2=OFF
    -DBuildCLAP=OFF
    -DBuildVST2=OFF
    -DBuildVST3=OFF
```

By default the Bypass will pass through the voice (modulator). You can build with -DUseCarrierBypass=ON if you prefer
the bypass use the carrier.

The default build uses the PFFFT library. You can build using the FFTW library using -DUsePFFFT=OFF.
