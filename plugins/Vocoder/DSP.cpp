/* 
    LV2 plugin for pitch shifting, pitch correction, vocoding and harmonizing
    singing voice

    (c) Igor Brkic 2010

    Phase vocoder code adapted from:
    http://www.dspdimension.com/admin/pitch-shifting-using-the-ft/

    
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


#include "DSP.hpp"

#define M_2PI 6.2831853071795
#define MAX_FRAME_LENGTH 4096

#define ROUND(a) ((float)((int)(a+0.5)))

static pthread_mutex_t fftw_planner_lock = PTHREAD_MUTEX_INITIALIZER;

// initialization
VocProc::VocProc(double rate)
{
    fSamplingFreq = (float)rate;

    sPitchFactor = 1.0;

    sEffect=0.0;
    cEffect=0.0;

    sOutputGain=1.0;

    cFormantVoco=1;

    pOffset[0]=0; pOffset[1]=0;

    sSwitch=1.0;

    cAutoTune=0.0;

    powerIn=0;

    fftFrameSize=2048;  // pitch detection currently doesn't work for 1024
                        // and there is aliasing present for 1024 with formant correction when
                        // large pitch shifting factor is applied
    overlap=4;

    freqOld=0;

    window=(float*)malloc(fftFrameSize*sizeof(float));
    for(int k=0;k<fftFrameSize;k++)
        window[k] = -.5*cos(M_2PI*(float)k/(float)fftFrameSize)+.5;

    gInFIFO=(float*)calloc(fftFrameSize, sizeof(float));
    gIn2FIFO=(float*)calloc(fftFrameSize, sizeof(float));
    gOutFIFO=(float*)calloc(fftFrameSize, sizeof(float));
    gOutputAccum=(float*)calloc(2*fftFrameSize, sizeof(float));

    // FFTW stuff
    fftTmpR=(double*)fftw_malloc(sizeof(double)*fftFrameSize);
    fftTmpC=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*fftFrameSize);
    fftOldC=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*fftFrameSize);
    fftCeps=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*fftFrameSize);

    pthread_mutex_lock (&fftw_planner_lock);
    fftPlan1=fftw_plan_dft_r2c_1d(fftFrameSize, fftTmpR, fftTmpC, FFTW_ESTIMATE);
    fftPlan2=fftw_plan_dft_r2c_1d(fftFrameSize, fftTmpR, fftTmpC, FFTW_ESTIMATE);
    fftPlanInv=fftw_plan_dft_c2r_1d(fftFrameSize, fftTmpC, fftTmpR, FFTW_ESTIMATE);
    pthread_mutex_unlock (&fftw_planner_lock);
    
}
VocProc::~VocProc()
{
    pthread_mutex_lock (&fftw_planner_lock);
    fftw_destroy_plan(fftPlan1);
    fftw_destroy_plan(fftPlan2);
    fftw_destroy_plan(fftPlanInv);
    pthread_mutex_unlock (&fftw_planner_lock);

    free(gInFIFO);
    free(gIn2FIFO);
    free(gOutFIFO);
    free(gOutputAccum);

    fftw_free(fftTmpR);
    fftw_free(fftTmpC);
    fftw_free(fftOldC);
    fftw_free(fftCeps);
};


/*************************************************
    here be procezzing
*************************************************/

void VocProc::run(const float **inputs, float **outputs, uint32_t nframes)
{
    const float *const input0 = inputs[0];

    const float *const input1 = inputs[1];

    float *output0 = outputs[0];

    static float gLastPhase[MAX_FRAME_LENGTH/2+1];
    static float gSumPhase[MAX_FRAME_LENGTH/2+1];
    static float gAnaFreq[MAX_FRAME_LENGTH];
    static float gAnaMagn[MAX_FRAME_LENGTH];
    static float gSynFreq[MAX_FRAME_LENGTH];
    static float gSynMagn[MAX_FRAME_LENGTH];

    static long gRover = false, gInit = false;

    double freqPerBin, expct;
    long i,k, index, inFifoLatency, stepSize, fftFrameSize2;

    float *fPointer, *fPointer2, *fPointer3;
    double *dPointer;

    /* set up some handy variables */
    fftFrameSize2 = fftFrameSize/2;

    stepSize = fftFrameSize/overlap;
    freqPerBin = (double)fSamplingFreq/(double)fftFrameSize;
    expct = M_2PI*(double)stepSize/(double)fftFrameSize;

    inFifoLatency = fftFrameSize-stepSize;
    if (gRover == false) gRover = inFifoLatency;

    /* initialize our static arrays */
    if (gInit == false) {
        memset(gLastPhase, 0, (MAX_FRAME_LENGTH/2+1)*sizeof(float));
        memset(gSumPhase, 0, (MAX_FRAME_LENGTH/2+1)*sizeof(float));
        memset(gAnaFreq, 0, MAX_FRAME_LENGTH*sizeof(float));
        memset(gAnaMagn, 0, MAX_FRAME_LENGTH*sizeof(float));
        gInit = true;
    }

    /* main processing loop */
    for (i = 0; i < nframes; i++){

        // As long as we have not yet collected enough data just read in
        gInFIFO[gRover] = input0[i];
        gIn2FIFO[gRover] = input1[i];

        output0[i] = gOutFIFO[gRover-inFifoLatency];
        gRover++;

        // now we have enough data for processing
        if (gRover >= fftFrameSize) {
            gRover = inFifoLatency;

            float tmpPower=0.0;
            dPointer=fftTmpR;
            fPointer=gInFIFO;
            fPointer2=window;
            for (k = 0; k < fftFrameSize;k++) {
                *dPointer=*(fPointer++) * *(fPointer2++);
                tmpPower+= *dPointer * *dPointer;
                dPointer++;
            }

            powerIn=tmpPower/(float)fftFrameSize;

            // do transform
            fftw_execute(fftPlan1);
 
            memcpy(fftOldC, fftTmpC, fftFrameSize*sizeof(fftw_complex));

            // pitch shifting with phase vocoder
            phaseVocAnalysis(fftTmpC, gLastPhase, freqPerBin, expct, gAnaMagn, gAnaFreq);
            memset(gSynMagn, 0, fftFrameSize*sizeof(float));
            memset(gSynFreq, 0, fftFrameSize*sizeof(float));
            for (k = 0; k <= fftFrameSize2; k++) {
                index = k*sPitchFactor;
                if (index <= fftFrameSize2) {
                    gSynMagn[index] += gAnaMagn[k];
                    gSynFreq[index] = gAnaFreq[k] * sPitchFactor;
                    if(cEffect)
                        gSynFreq[index] = gSynFreq[index]*sEffect + sEffect*200*(float)rand()/RAND_MAX-100;
                }
            }
            phaseVocSynthesis(fftTmpC, gSumPhase, gSynMagn, gSynFreq, freqPerBin, expct);

            // formant correction + vocoder
            if(cFormantVoco)
            {
                float env1[fftFrameSize2], env2[fftFrameSize2];

                if(sSwitch)
                {
                    dPointer=fftTmpR; fPointer=gIn2FIFO; fPointer2=window;
                    for (k = 0; k < fftFrameSize;k++) {
                        *dPointer++=*(fPointer++) * *(fPointer2++);
                    }

                    fftw_execute(fftPlan2);
                }

                spectralEnvelope(env1, fftOldC, fftFrameSize2);
                spectralEnvelope(env2, fftTmpC, fftFrameSize2);

                // modify spectral envelope of spectrum in fftTmpC to look like spectral
                // envelope of spectrum in fftOldC
                float koef;
                fPointer2=env1; fPointer3=env2;
                for(k=0;k<fftFrameSize2;k++){
                    koef = *(fPointer2++) / (*(fPointer3++)+.02)*2;
                    fftTmpC[k][0] *= koef;
                    fftTmpC[k][1] *= koef;
                }
            }

            // do inverse transform
            fftw_execute(fftPlanInv);

            fPointer=gOutputAccum; dPointer=fftTmpR; fPointer2=window;
            for(k=0; k < fftFrameSize; k++) {
                *fPointer += 0.7 * *(dPointer++) / (fftFrameSize2*overlap) * sOutputGain * *(fPointer2++);
                fPointer++;
            }

            memcpy(gOutFIFO, gOutputAccum, stepSize*sizeof(float));

            // shift accumulator
            memmove(gOutputAccum, gOutputAccum+stepSize, fftFrameSize*sizeof(float));

            // move input FIFO
            for (k = 0; k < inFifoLatency; k++) gInFIFO[k] = gInFIFO[k+stepSize];

            for (k = 0; k < inFifoLatency; k++) gIn2FIFO[k] = gIn2FIFO[k+stepSize];
        }
    }
}


void VocProc::phaseVocAnalysis(fftw_complex *block, float *gLastPhase, double freqPerBin, double expct, float *gAnaMagn, float *gAnaFreq){

    double real, imag, magn, phase, tmp;
    long qpd, k;

    for (k = 0; k <= fftFrameSize/2; k++) {

        /* de-interlace FFT buffer */
        real = block[k][0];
        imag = block[k][1];

        /* compute magnitude and phase */
        magn = 2.*sqrt(real*real + imag*imag);
        phase = atan2(imag,real);

        /* compute phase difference */
        tmp = phase - gLastPhase[k];
        gLastPhase[k] = phase;

        /* subtract expected phase difference */
        tmp -= (double)k*expct;

        /* map delta phase into +/- Pi interval */
        qpd = tmp/M_PI;
        if (qpd >= 0) qpd += qpd&1;
        else qpd -= qpd&1;
        tmp -= M_PI*(double)qpd;

        /* get deviation from bin frequency from the +/- Pi interval */
        tmp = overlap*tmp/(M_2PI);

        /* compute the k-th partials' true frequency */
        tmp = (double)k*freqPerBin + tmp*freqPerBin;

        /* store magnitude and true frequency in analysis arrays */
        gAnaMagn[k] = magn;
        gAnaFreq[k] = tmp;
    }
}

void VocProc::phaseVocSynthesis(fftw_complex *block, float *gSumPhase, float *gSynMagn, float *gSynFreq, double freqPerBin, double expct){

    int k;
    double magn, tmp, phase;

    for (k = 0; k <= fftFrameSize/2; k++) {
        /* get magnitude and true frequency from synthesis arrays */
        magn = gSynMagn[k];
        tmp = gSynFreq[k];

        /* subtract bin mid frequency */
        tmp -= (double)k*freqPerBin;

        /* get bin deviation from freq deviation */
        tmp /= freqPerBin;

        /* take overlap into acnframes */
        tmp = M_2PI*tmp/overlap;

        /* add the overlap phase advance back in */
        tmp += (double)k*expct;

        /* accumulate delta phase to get bin phase */
        gSumPhase[k] += tmp;
        phase = gSumPhase[k];

        /* get real and imag part and re-interleave */
        block[k][0] = magn*cos(phase);
        block[k][1] = magn*sin(phase);
    }
}

void VocProc::spectralEnvelope(float *env, fftw_complex *fft, uint32_t nframes){

    unsigned int nTaps=20;
    unsigned int nTaps2=10;
    float tmp[nframes+nTaps];

    float h[]={ // h=(firls(20, [0 0.02 0.1 1], [1 1 0 0]));
        0.0180, 0.0243, 0.0310, 0.0378, 0.0445, 0.0508, 0.0564, 0.0611,
        0.0646, 0.0667, 0.0675, 0.0667, 0.0646, 0.0611, 0.0564, 0.0508,
        0.0445, 0.0378, 0.0310, 0.0243, 0.0180
    };

    // |H(w)|
    memset(tmp,  0, (nframes+nTaps)*sizeof(float));
    for(unsigned int k=0;k<nframes;k++)
        tmp[k]=sqrt(fft[k][0]*fft[k][0]+fft[k][1]*fft[k][1]);

    memset(env, 0, nframes*sizeof(float));

    // magnitude spectrum filtering
    unsigned int i, j;
    float *p_h, *p_z, accum;

    float z[2 * nTaps];
    memset(z, 0, 2*nTaps*sizeof(float));
    int state = 0;
    for (j = 0; j < nframes+nTaps2; j++) {
        z[state] = z[state + nTaps] = tmp[j];
        p_h = h;
        p_z = z + state;
        accum = 0;
        for (i = 0; i < nTaps; i++) accum += *p_h++ * *p_z++;
        if (--state < 0) state += nTaps;
        if(j>=nTaps2) env[j-nTaps2]=accum;
    }
}

void VocProc::set_bypass(float bypass)
{
    if ( bypass )
    {
        sSwitch = 0.0;
    }
    else
        sSwitch = 1.0;
}