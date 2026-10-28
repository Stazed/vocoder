/*
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* 
 * File:   DSP.hpp
 * Author: sspresto
 *
 * Created on December 9, 2023, 2:39 PM
 */

#ifndef DSP_HPP
#define DSP_HPP

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <pffft/pffft.h>

#define MAX_FRAME_LENGTH 4096

class VocProc
{
public:
    explicit VocProc(double rate);
    ~VocProc();

    void run(const float **inputs, float **outputs, uint32_t nframes);
    void set_bypass(float bypass);

private:
    // Parameters
    float fSamplingFreq;

    float sPitchFactor;
    float sEffect;
    float sOutputGain;
    bool  bypassed;
    bool  bypassedPrev;

    float cFormantVoco;
    float cEffect;
    float cAutoTune;

    float powerIn;

    // FFT parameters
    long fftFrameSize;
    long overlap;

    // Time-domain buffers
    float *gInFIFO;
    float *gIn2FIFO;
    float *gOutFIFO;
    float *gOutputAccum;
    float *window;

    float gLastPhase[MAX_FRAME_LENGTH/2 + 1];
    float gSumPhase [MAX_FRAME_LENGTH/2 + 1];
    float gAnaFreq  [MAX_FRAME_LENGTH];
    float gAnaMagn  [MAX_FRAME_LENGTH];
    float gSynFreq  [MAX_FRAME_LENGTH];
    float gSynMagn  [MAX_FRAME_LENGTH];

    long gRover;
    bool gInit;
    uint32_t rngState;

    typedef float fft_complex_t; // interleaved real/imag
    PFFFT_Setup* fftSetup;
    float *fftTmpR;
    float *fftTmpC;   // interleaved complex
    float *fftOldC;

    // DSP helpers
    void spectralEnvelope(float *env,
                          const fft_complex_t *fft,
                          uint32_t nframes);

    void phaseVocAnalysis(fft_complex_t *block,
                          float *lastPhase,
                          double freqPerBin,
                          double expct,
                          float *anaMagn,
                          float *anaFreq);

    void phaseVocSynthesis(fft_complex_t *block,
                           float *sumPhase,
                           const float *synMagn,
                           const float *synFreq,
                           double freqPerBin,
                           double expct);
    void resetDSPState();
    inline float noise();
};

#endif // DSP_HPP
