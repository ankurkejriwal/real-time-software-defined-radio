
#ifndef DY4_SUPPORT_H
#define DY4_SUPPORT_H

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

struct pll_state{
		float integrator, phaseEst, feedbackI, feedbackQ, trigOffset,ncoLast;
	};


void vectorSlice(const std::vector<float> &inp, std::vector<float> &out, int start, int end, int step);
void vectorSplit(const std::vector<float> &inp, std::vector<float> &Aout, std::vector<float> &Bout);
void fmDemodArctan(std::vector<float> &fm_demod, std::vector<float> &I, std::vector<float> &Q, std::vector<float> &prev_phase);
void normalized(std::vector<float> &data);
void fmPLL(std::vector<float> &ncoOut, std::vector<float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth, pll_state &state);
void pointwiseMult(std::vector<float> &a, std::vector<float> &b, std::vector<float> &out);

void estimatePSD(std::vector<float> &samples, int Fs, std::vector<float> &psd_est, std::vector<float> &freq);

#endif