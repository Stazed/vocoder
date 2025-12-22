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

#include <fftw3.h>
#include <pthread.h>

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
    float sSwitch;

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

    // FFTW
    double       *fftTmpR;
    fftw_complex *fftTmpC;   // WORK BUFFER (modulator OR carrier)
    fftw_complex *fftOldC;   // saved modulator spectrum

    fftw_plan fftPlanFwd;
    fftw_plan fftPlanInv;

    // DSP helpers
    void spectralEnvelope(float *env,
                          const fftw_complex *fft,
                          uint32_t nframes);

    void phaseVocAnalysis(fftw_complex *block,
                          float *lastPhase,
                          double freqPerBin,
                          double expct,
                          float *anaMagn,
                          float *anaFreq);

    void phaseVocSynthesis(fftw_complex *block,
                           float *sumPhase,
                           const float *synMagn,
                           const float *synFreq,
                           double freqPerBin,
                           double expct);
};

#endif // DSP_HPP
