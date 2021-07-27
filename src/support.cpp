#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include "support.h"
#include "dy4.h"
#include "fourier.h"
#include "iofunc.h"
#include "logfunc.h"


// for decimation
void vectorSlice(const std::vector<float> &inp, std::vector<float> &out, int start, int end, int step) 
{
    out.resize(inp.size()/step);
    int j =0;
    for(int i=start; i<end; i+=step) {
        out[j] = inp[i];
        j++;
    }
}

// for i and q samples from block
void vectorSplit(const std::vector<float> &inp, std::vector<float> &Aout, std::vector<float> &Bout) 
{
    Aout.resize(inp.size()/2);
    Bout.resize(inp.size()/2);
    int j =0;
    for(int i=0; i<inp.size(); i+=2) {
        Aout[j] = inp[i];
        Bout[j] = inp[i+1];
        j++;
    }
}

void fmDemodArctan(std::vector<float> &fm_demod, std::vector<float> &I, std::vector<float> &Q, std::vector<float> &prev_phase)
{
    float dQ;
    float dI;
    float scaling;
    float phase;
    // prev_phase = {0, 0};
    fm_demod.resize(I.size());
    for(int i =0; i<I.size(); i++){
        dQ = Q[i] - prev_phase[1];
        dI = I[i] - prev_phase[0];
        if(pow(I[i],2) + pow(Q[i],2)==0){
            phase = 0;
            // fm_demod.push_back(0);
            fm_demod[i] = 0;
        } else {
            scaling = 1/(pow(I[i], (float)2) + pow(Q[i], (float)2));
            phase = (I[i]*dQ - Q[i]*dI)*scaling;
            fm_demod[i]=phase;
        }
        
        prev_phase = {I[i], Q[i]};
    }
}


void normalized(std::vector<float> &data){
    
    for(int i=0; i < data.size(); i++){
        data[i] = (data[i]-(float)128)/(float)128;
    }
}

void fmPLL(std::vector<float> &ncoOut, std::vector<float> &pllIn, float freq, float Fs, float ncoScale, float phaseAdjust, float normBandwidth , pll_state &state){
    
    float Cp = 2.666;
	float Ci = 3.555;

	// gain for the proportional term
	float Kp = (normBandwidth)*Cp;
	// gain for the integrator term
	float Ki = (normBandwidth*normBandwidth)*Ci;

	// output array for the NCO
	ncoOut.resize(pllIn.size()+1);

    //Intiliazing variables used inside PLL to state saved values
    float integrator = state.integrator;
    float phaseEst = state.phaseEst;
    float feedbackI = state.feedbackI;
    float feedbackQ = state.feedbackQ;
    ncoOut[0] = state.ncoLast;
    float trigOffset = state.trigOffset;

    // float phaseAdjust = 0.0;

	for (int k=0; k<pllIn.size(); k++)
    {
		float errorI = pllIn[k] * (+feedbackI);  // complex conjugate of the
		float errorQ = pllIn[k] * (-feedbackQ); // feedback complex exponential
    
    	// four-quadrant arctangent discriminator for phase error detection
		float errorD = atan2(errorQ, errorI);

		// loop filter
		integrator += Ki*errorD;

		// update phase estimate
		phaseEst += Kp*errorD + integrator;

		// internal oscillator
		float trigArg = 2*PI*(freq/Fs)*(trigOffset+k+1)+phaseEst;
		feedbackI = cos(trigArg);
		feedbackQ = sin(trigArg);
		ncoOut[k+1] = cos(trigArg*ncoScale + phaseAdjust);
    }

    //Update State Variables so they are saved to the struct object in main    
    state.integrator = integrator;
    state.phaseEst = phaseEst;
    state.feedbackI = feedbackI;
    state.feedbackQ = feedbackQ;
    state.ncoLast= ncoOut[ncoOut.size()-1];
    state.trigOffset= trigOffset + pllIn.size();

    //Resize to return 1:end of array
    ncoOut = std::vector<float>(ncoOut.begin(), ncoOut.end()-1);

}

void pointwiseMult(std::vector<float> &a, std::vector<float> &b, std::vector<float> &out){
    out.resize(a.size());
    for (int i=0; i<a.size(); i++){
        out[i] = a[i]*b[i];
    }
}


void estimatePSD(std::vector<float> &samples, int Fs, std::vector<float> &psd_est, std::vector<float> &freq)
{
	int freq_bins = NFFT;
 
	float df = (float)Fs/(float)freq_bins;
 
	freq.resize(freq_bins/2);
 
	// populate vector to hold x values of plot
	for (int i = 0; i < freq.size(); i++)
	{
		freq[i] = i*df;
	}
 
	std::vector<float> hann(freq_bins);
        

 
	// print out vector for X values
	for (int j = 0; j < hann.size(); j++)
	{
		hann[j] = pow((float)sin(j*PI/(float)freq_bins), 2);
	}
 
	// empty list for PSD segments
	std::vector<float> psd_list;
 
	int no_segments = floor(samples.size()/(float)freq_bins);
	for (int k = 0; k < no_segments; k++) {
		// apply the hann window using pointwise multiplication
		std::vector<float> windowed_samples(freq_bins);
		std::vector<float> multiplier_vector = std::vector<float>(samples.begin()+(k*freq_bins), samples.begin()+((k+1)*freq_bins));
		// multiply each point in the range of values by hann
		for (int i = 0; i < hann.size(); i++){
			windowed_samples[i] = multiplier_vector[i]*hann[i];
		}

		std::vector<std::complex<float>> Xf;
		DFT(windowed_samples, Xf);

		// to have a better and more accurate PSD estimate when plotting
		Xf = std::vector<std::complex<float>>(Xf.begin(), Xf.begin()+((int)(freq_bins/2)));

		// compute signal power and add energy from the negative freq bins
		// then translate to the decibel scale
		std::vector<float> psd_seg(Xf.size());
		for (int i=0; i<Xf.size(); i++){
			psd_seg[i] = 10*log10(2* pow(abs(Xf[i]), 2)/((float)Fs*(float)freq_bins/2));
		}
		// append to the list where PSD for each segment is stored in sequential order
		psd_list.insert(psd_list.end(), psd_seg.begin(), psd_seg.end());
	}
 
	psd_est.resize((int)(freq_bins/2));
	std::cout<<" PSD size: " <<psd_est.size()<<"\n";
	for (int k = 0; k < (int)(freq_bins/2); k++) {
		for (int l=0; l<no_segments; l++) {
			psd_est[k] += psd_list[k+l*((int)(freq_bins/2))];
		}
		psd_est[k] = psd_est[k]/(float)no_segments;
	}
}