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
#include <cassert>

#define M_2PI 6.283185307179586

#ifndef PFFFT_SUPPORT
static pthread_mutex_t fftw_planner_lock = PTHREAD_MUTEX_INITIALIZER;
#endif

// ---------------- FFT helpers ----------------
#ifdef PFFFT_SUPPORT
    #define FFT_REAL(buf,k) buf[2*(k)]
    #define FFT_IMAG(buf,k) buf[2*(k)+1]
#else
    #define FFT_REAL(buf,k) buf[k][0]
    #define FFT_IMAG(buf,k) buf[k][1]
#endif

// ============================================================
// Constructor / Destructor
// ============================================================

VocProc::VocProc(double rate)
{
    fSamplingFreq = (float)rate;

    sPitchFactor = 1.0f;
    sEffect      = 0.0f;
    sOutputGain  = 1.0f;
    bypassed     = false;
    bypassedPrev = false;

    cFormantVoco = 1.0f;
    cEffect      = 0.0f;
    cAutoTune    = 0.0f;

    powerIn = 0.0f;

    fftFrameSize = 2048;
    overlap      = 4;

    assert(fftFrameSize <= MAX_FRAME_LENGTH);

    // Hann window
    window = (float*)malloc(sizeof(float) * fftFrameSize);
    for (long k = 0; k < fftFrameSize; ++k)
        window[k] = 0.5f - 0.5f * cos(M_2PI * k / fftFrameSize);

    gInFIFO      = (float*)calloc(fftFrameSize, sizeof(float));
    gIn2FIFO     = (float*)calloc(fftFrameSize, sizeof(float));
    gOutFIFO     = (float*)calloc(fftFrameSize, sizeof(float));
    gOutputAccum = (float*)calloc(2 * fftFrameSize, sizeof(float));

    memset(gLastPhase, 0, sizeof(gLastPhase));
    memset(gSumPhase,  0, sizeof(gSumPhase));
    memset(gAnaFreq, 0, sizeof(gAnaFreq));
    memset(gAnaMagn, 0, sizeof(gAnaMagn));
    memset(gSynFreq, 0, sizeof(gSynFreq));
    memset(gSynMagn, 0, sizeof(gSynMagn));

    gRover = 0;
    gInit  = false;
    rngState = 0x12345678u ^ (uintptr_t)this;

#ifdef PFFFT_SUPPORT
    fftSetup = pffft_new_setup(fftFrameSize, PFFFT_REAL);
    fftTmpR  = (float*)pffft_aligned_malloc(sizeof(float) * fftFrameSize);
    fftTmpC  = (float*)pffft_aligned_malloc(sizeof(float) * fftFrameSize);
    fftOldC  = (float*)pffft_aligned_malloc(sizeof(float) * fftFrameSize);
#else
    fftTmpR = (double*)malloc(sizeof(double) * fftFrameSize);
    fftTmpC = (fft_complex_t*)malloc(sizeof(fft_complex_t) * (fftFrameSize/2+1));
    fftOldC = (fft_complex_t*)malloc(sizeof(fft_complex_t) * (fftFrameSize/2+1));

    pthread_mutex_lock(&fftw_planner_lock);
    fftPlanFwd = fftw_plan_dft_r2c_1d(
        fftFrameSize, fftTmpR, fftTmpC, FFTW_ESTIMATE);
    fftPlanInv = fftw_plan_dft_c2r_1d(
        fftFrameSize, fftTmpC, fftTmpR, FFTW_ESTIMATE);
    pthread_mutex_unlock(&fftw_planner_lock);
#endif
}

VocProc::~VocProc()
{
#ifdef PFFFT_SUPPORT
    pffft_destroy_setup(fftSetup);
    pffft_aligned_free(fftTmpR);
    pffft_aligned_free(fftTmpC);
    pffft_aligned_free(fftOldC);
#else
    pthread_mutex_lock(&fftw_planner_lock);
    fftw_destroy_plan(fftPlanFwd);
    fftw_destroy_plan(fftPlanInv);
    pthread_mutex_unlock(&fftw_planner_lock);
    free(fftTmpR);
    free(fftTmpC);
    free(fftOldC);
#endif

    free(window);
    free(gInFIFO);
    free(gIn2FIFO);
    free(gOutFIFO);
    free(gOutputAccum);
}

// ============================================================
// Main Processing
// ============================================================

void VocProc::run(const float **inputs, float **outputs, uint32_t nframes)
{
    if (bypassed && !bypassedPrev) {
        // Entering bypass: reset DSP so it doesn't run silently
        resetDSPState();
    }

    if (!bypassed && bypassedPrev) {
        // Leaving bypass: clean DSP restart
        resetDSPState();
    }

    bypassedPrev = bypassed;

    const float *inputMod = inputs[0];
    const float *inputCar = inputs[1];
    float *output = outputs[0];
 
    const long fftFrameSize2 = fftFrameSize / 2;
    const long stepSize = fftFrameSize / overlap;
    const long fifoLatency = fftFrameSize - stepSize;

    const double freqPerBin = fSamplingFreq / fftFrameSize;
    const double expct = M_2PI * stepSize / fftFrameSize;

    if (!gInit) {
        memset(gLastPhase, 0, sizeof(gLastPhase));
        memset(gSumPhase,  0, sizeof(gSumPhase));
        gRover = fifoLatency;
        gInit = true;
    }

    for (uint32_t i = 0; i < nframes; ++i) {
        if (bypassed) {
#ifdef USE_CARRIER_ON_BYPASS
            output[i] = inputCar[i];
#else
            output[i] = inputMod[i];
#endif
            continue;
        }

        gInFIFO [gRover] = inputMod[i];
        gIn2FIFO[gRover] = inputCar[i];

        output[i] = gOutFIFO[gRover - fifoLatency];
        gRover++;

        if (gRover < fftFrameSize)
            continue;

        gRover = fifoLatency;

        // ================= Modulator FFT =================
        double pwr = 0.0;
        for (long k = 0; k < fftFrameSize; ++k) {
            fftTmpR[k] = gInFIFO[k] * window[k];
            pwr += fftTmpR[k] * fftTmpR[k];
        }
        powerIn = (float)(pwr / fftFrameSize);

#ifdef PFFFT_SUPPORT
        pffft_transform_ordered(fftSetup, fftTmpR, fftTmpC, nullptr, PFFFT_FORWARD);
#else
        fftw_execute(fftPlanFwd);
#endif

#ifdef PFFFT_SUPPORT
        memcpy(fftOldC, fftTmpC, sizeof(float) * fftFrameSize);
#else
        memcpy(fftOldC, fftTmpC, sizeof(fft_complex_t) * (fftFrameSize2 + 1));
#endif

        // ================= Phase Vocoder =================
        phaseVocAnalysis(
            fftTmpC, gLastPhase,
            freqPerBin, expct,
            gAnaMagn, gAnaFreq
        );

        memset(gSynMagn, 0, sizeof(float) * fftFrameSize);
        memset(gSynFreq, 0, sizeof(float) * fftFrameSize);

        for (long k = 0; k <= fftFrameSize2; ++k) {
            long idx = (long)(k * sPitchFactor);
            if (idx <= fftFrameSize2) {
                gSynMagn[idx] += gAnaMagn[k];
                gSynFreq[idx]  = gAnaFreq[k] * sPitchFactor;
                if (cEffect)
                    gSynFreq[idx] += sEffect * 100.0f * noise();
            }
        }

        phaseVocSynthesis(
            fftTmpC, gSumPhase,
            gSynMagn, gSynFreq,
            freqPerBin, expct
        );

        // ================= Carrier FFT =================
        for (long k = 0; k < fftFrameSize; ++k)
            fftTmpR[k] = gIn2FIFO[k] * window[k];

#ifdef PFFFT_SUPPORT
        pffft_transform_ordered(fftSetup, fftTmpR, fftTmpC, nullptr, PFFFT_FORWARD);
#else
        fftw_execute(fftPlanFwd);
#endif

        // ================= Envelope Transfer =================
        if (cFormantVoco) {
            float envMod[MAX_FRAME_LENGTH/2];
            float envCar[MAX_FRAME_LENGTH/2];

            spectralEnvelope(envMod, fftOldC, fftFrameSize2);
            spectralEnvelope(envCar, fftTmpC, fftFrameSize2);

            for (long k = 0; k < fftFrameSize2; ++k) {
                float coef = envMod[k] / (envCar[k] + 0.02f) * 2.0f;
                FFT_REAL(fftTmpC,k) *= coef;
                FFT_IMAG(fftTmpC,k) *= coef;
            }
        }

        // ================= IFFT + OLA =================
#ifdef PFFFT_SUPPORT
        pffft_transform_ordered(fftSetup, fftTmpC, fftTmpR, nullptr, PFFFT_BACKWARD);
#else
        fftw_execute(fftPlanInv);
#endif

        for (long k = 0; k < fftFrameSize; ++k) {
            gOutputAccum[k] +=
                0.7f * fftTmpR[k] * window[k]
                / (fftFrameSize2 * overlap)
                * sOutputGain;
        }

        memcpy(gOutFIFO, gOutputAccum,
               sizeof(float) * stepSize);

        memmove(gOutputAccum,
                gOutputAccum + stepSize,
                sizeof(float) * fftFrameSize);

        memmove(gInFIFO,
                gInFIFO + stepSize,
                sizeof(float) * fifoLatency);
        memmove(gIn2FIFO,
                gIn2FIFO + stepSize,
                sizeof(float) * fifoLatency);
    }
}

// ============================================================
// Phase Vocoder Analysis
// ============================================================

void VocProc::phaseVocAnalysis(
    fft_complex_t *block,
    float *lastPhase,
    double freqPerBin,
    double expct,
    float *anaMagn,
    float *anaFreq)
{
    for (long k = 0; k <= fftFrameSize / 2; ++k) {
        double real = FFT_REAL(block,k);
        double imag = FFT_IMAG(block,k);

        double magn = 2.0 * sqrt(real*real + imag*imag);
        double phase = atan2(imag, real);

        double tmp = phase - lastPhase[k];
        lastPhase[k] = (float)phase;

        tmp -= k * expct;

        long qpd = (long)(tmp / M_PI);
        if (qpd >= 0) qpd += qpd & 1;
        else          qpd -= qpd & 1;
        tmp -= M_PI * qpd;

        tmp = overlap * tmp / M_2PI;
        tmp = (k + tmp) * freqPerBin;

        anaMagn[k] = (float)magn;
        anaFreq[k] = (float)tmp;
    }
}

// ============================================================
// Phase Vocoder Synthesis
// ============================================================

void VocProc::phaseVocSynthesis(
    fft_complex_t *block,
    float *sumPhase,
    const float *synMagn,
    const float *synFreq,
    double freqPerBin,
    double expct)
{
    for (long k = 0; k <= fftFrameSize / 2; ++k) {
        double magn = synMagn[k];
        double tmp  = synFreq[k] - k * freqPerBin;

        tmp /= freqPerBin;
        tmp = M_2PI * tmp / overlap;
        tmp += k * expct;

        sumPhase[k] += (float)tmp;
        FFT_REAL(block,k) = magn * cos(sumPhase[k]);
        FFT_IMAG(block,k) = magn * sin(sumPhase[k]);
    }
}

// ============================================================
// Spectral Envelope
// ============================================================

void VocProc::spectralEnvelope(
    float *env,
    const fft_complex_t *fft,
    uint32_t nframes)
{
    const unsigned int nTaps  = 20;
    const unsigned int nTaps2 = 10;

    static const float h[21] = {
        0.0180f, 0.0243f, 0.0310f, 0.0378f, 0.0445f,
        0.0508f, 0.0564f, 0.0611f, 0.0646f, 0.0667f,
        0.0675f, 0.0667f, 0.0646f, 0.0611f, 0.0564f,
        0.0508f, 0.0445f, 0.0378f, 0.0310f, 0.0243f,
        0.0180f
    };

    float tmp[MAX_FRAME_LENGTH + 21];
    memset(tmp, 0, sizeof(float) * (nframes + nTaps));
    memset(env, 0, sizeof(float) * nframes);

    for (uint32_t k = 0; k < nframes; ++k)
        tmp[k] = sqrt(FFT_REAL(fft,k)*FFT_REAL(fft,k) + FFT_IMAG(fft,k)*FFT_IMAG(fft,k));

    float z[40] = {};
    int state = 0;

    for (uint32_t j = 0; j < nframes + nTaps2; ++j) {
        z[state] = z[state + nTaps] = tmp[j];

        float acc = 0.0f;
        for (unsigned int i = 0; i < nTaps; ++i)
            acc += h[i] * z[state + i];

        if (--state < 0) state += nTaps;
        if (j >= nTaps2) env[j - nTaps2] = acc;
    }
}

// ============================================================
// Bypass
// ============================================================

void VocProc::set_bypass(float bypass)
{
    bypassed = bypass > 0.5f;
}

// ============================================================
// Reset
// ============================================================

void VocProc::resetDSPState()
{
    gRover = fftFrameSize / overlap;
    gInit  = false;
    memset(gLastPhase, 0, sizeof(gLastPhase));
    memset(gSumPhase, 0, sizeof(gSumPhase));
    memset(gInFIFO, 0, sizeof(float)*fftFrameSize);
    memset(gIn2FIFO, 0, sizeof(float)*fftFrameSize);
    memset(gOutFIFO, 0, sizeof(float)*fftFrameSize);
    memset(gOutputAccum, 0, sizeof(float)*2*fftFrameSize);
}

// ============================================================
// Noise Helper
// ============================================================

inline float VocProc::noise()
{
    rngState ^= rngState << 13;
    rngState ^= rngState >> 17;
    rngState ^= rngState << 5;
    return (rngState * (1.0f / UINT32_MAX)) * 2.0f - 1.0f;
}
