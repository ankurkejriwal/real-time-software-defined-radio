/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>

// declaration of a function prototypes
void impulseResponseLPF(float, float, unsigned short int, std::vector<float> &);
void convolution(std::vector<float> &, const std::vector<float> &, const std::vector<float> &,int decim);
void convolutionBlock(std::vector<float> &y,const std::vector<float> &x, const std::vector<float> &h,std::vector<float> &state,int decim);
void convolutionBlockIQ(std::vector<float> &I, std::vector<float> &Q, const std::vector<float> &I_data, const std::vector<float> &Q_data, const std::vector<float> &h,std::vector<float> &stateI, std::vector<float> &stateQ,int decim);
void impulseResponseBPF(float Fb, float Fe, float Fs, int n_taps, std::vector<float> &h);
void mode1Convolution(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, std::vector<float> &state, int decim, int upsample);
#endif // DY4_FILTER_H
