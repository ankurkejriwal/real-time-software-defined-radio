/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include <cmath>

void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	// allocate memory for the impulse response
	h.resize(num_taps, 0.0);

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	float norm_cutoff = Fc/(Fs/2);
	float Ndiv = (num_taps-1)/2;
	for(int i=0; i<num_taps; i++){
		if (i == Ndiv){
			h[i] = norm_cutoff;
		} else {
			h[i] = norm_cutoff*sin(PI*norm_cutoff*(i-Ndiv))/(PI*norm_cutoff*(i-Ndiv));
		}
		h[i] =h[i]*(float)pow(sin(i*PI/num_taps),2.0);
	}
}


void impulseResponseBPF(float Fb, float Fe, float Fs, int n_taps, std::vector<float> &h){
	
	h.resize(n_taps, 0.0);
	float norm_center = ((Fe+Fb)/2)/(Fs/2);
    float norm_pass = (Fe-Fb)/(Fs/2);
	float n_half = (n_taps-1)/2;
	for(int i=0; i<n_taps; i++){
		if(i == n_half){
			h[i] = norm_pass;
		} else {
			h[i] = norm_pass * sin(PI*(norm_pass/2)*(i-n_half)) / (PI*(norm_pass/2)*(i-n_half));
		}
		h[i] = h[i] * cos(i*PI*norm_center);
        h[i] = h[i] * pow(sin((i*PI)/n_taps), 2);
	}
}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"; this is based on
// your Python code from the take-home exercise from the previous lab
void convolution(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h,int decim)
{
	// allocate memory for the output (filtered) data
	y.resize(ceil((x.size()+h.size()-1)/decim), 0.0);

	int y_non_decim_size = x.size()+h.size()-1;
	int d = 0;
	for(int i = 0; i<y_non_decim_size; i+=decim){
		for(int k = 0; k<h.size(); k++){
			if(i-k>=0 && i-k < x.size()){
				y[d] += x[i-k]*h[k];
			}
		}
		d++;
	}
	y = std::vector<float>(y.begin(),y.begin()+(x.size()/decim));
}

void convolutionBlock(std::vector<float> &y,const std::vector<float> &x, const std::vector<float> &h,std::vector<float> &state,const int decim){
	y.resize(x.size()/decim);
	int y_non_decim_size = x.size();
	int d = 0;
	int s = state.size()-1;
	int count;
	 for(int n = 0 ; n < y_non_decim_size; n+=decim){
		 count=0;
		 for(int k = 0 ; k < h.size(); k++){
			if (n-k >= 0 && n-k <x.size()){
				y[d] += h[k]*x[n-k];
			}else{
				y[d] += h[k]*state[s-count];
				count++;
			}
		}
		d++;
	}

	auto first = (x.size()-(h.size()-1));
	auto last = x.end();

	state = std::vector<float>(x.begin()+first, last);
}

// filter for concurrent i and q filtering
void convolutionBlockIQ(std::vector<float> &I, std::vector<float> &Q, const std::vector<float> &I_data, const std::vector<float> &Q_data, const std::vector<float> &h,std::vector<float> &stateI, std::vector<float> &stateQ,int decim){
	I.resize(I_data.size()/decim);
	Q.resize(Q_data.size()/decim);
	int d = 0;
	int s = stateI.size()-1;
	int count;
	 for(int n = 0 ; n < I_data.size(); n+=decim){
		 count=0;
		 for(int k = 0 ; k < h.size(); k++){
			if (n-k >= 0 && n-k <I_data.size()){
				I[d] += h[k]*I_data[n-k];
				Q[d] += h[k]*Q_data[n-k];
			}else{
				I[d] += h[k]*stateI[s-count];
				Q[d] += h[k]*stateQ[s-count];
				count++;
			}
		}
		d++;	
	}

	auto first = (I_data.size()-(h.size()-1));

	stateI = std::vector<float>(I_data.begin()+first, I_data.end());
	stateQ = std::vector<float>(Q_data.begin()+first, Q_data.end());

}

void mode1Convolution(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state, int decim, int upsample)
{
	y.resize(x.size()*upsample/decim);

	int s = state.size()-1;
	int count;

	for(int i=0; i < y.size(); i++){

		count = 0;

		// iterate through h vector with a step size of the upsample factor (24)
		for(int j=0; j < h.size(); j += upsample) {

			if(0 <= decim*i-j && decim*i-j < x.size()*upsample){
				
				y[i] += x[(decim*i-j)/upsample]*h[j];

			}
			else
			{
			
				y[i] += state[s-count] * h[j];
				count++;
			
			}

		}
	}

	auto first = (x.size()-(h.size()-1));
	auto last = x.end();

	state = std::vector<float>(x.begin()+first, last);
}