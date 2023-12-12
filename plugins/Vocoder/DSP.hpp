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

#include <stdint.h>
#include <fftw3.h>
#include <pthread.h>
#include <string.h>     // because of memset
#include <stdlib.h>
#include <math.h>

class VocProc
{
    private:
        float   fSamplingFreq;
        float 	sPitchFactor, sEffect, sOutputGain;
        float   cFormantVoco, cEffect;
        float   powerIn;
        float   sSwitch;
        float   pOffset[2];

        float   cAutoTune;

        float   *gInFIFO, *gIn2FIFO, *gOutFIFO, *gOutputAccum;
        float   *window;

        long    fftFrameSize, overlap;

        float   freqOld;

        double  *fftTmpR;
        fftw_complex *fftTmpC;
        fftw_complex *fftOldC;
        fftw_complex *fftCeps;
        fftw_plan fftPlan1;
        fftw_plan fftPlan2;
        fftw_plan fftPlanInv;

        void    spectralEnvelope(float *env, fftw_complex *fft, uint32_t nframes);

        float   pitchFrequency(fftw_complex *block);

        void    phaseVocAnalysis(fftw_complex *block, float *gLastPhase, double freqPerBin, double expct, float *gAnaMagn, float *gAnaFreq);
        void    phaseVocSynthesis(fftw_complex *block, float *gSumPhase, float *gSynMagn, float *gSynFreq, double freqPerBin, double expct);



    public:

        VocProc(double rate);
        ~VocProc();

        void run(const float **inputs, float **outputs, uint32_t nframes);
        void set_bypass(float bypass);
};


#endif /* DSP_HPP */

