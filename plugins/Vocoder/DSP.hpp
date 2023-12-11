/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
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

