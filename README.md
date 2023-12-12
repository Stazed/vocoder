
Vocoder
=======

Vocoder is a simple LADSPA/LV2/CLAP plugin for vocoding based on VocProc by Igor Brkic <igor@hyperglitch.com>.
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
