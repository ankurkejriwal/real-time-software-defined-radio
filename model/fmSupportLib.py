#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import numpy as np
import math, cmath

#
# you should use the demodulator based on arctan given below as a reference
#
# in order to implement your OWN FM demodulator without the arctan function,
# a very good and to-the-point description is given by Richard Lyons at:
#
# https://www.embedded.com/dsp-tricks-frequency-demodulation-algorithms/
#
# the demodulator boils down to implementing equation (13-117) from above, where
# the derivatives are nothing else but differences between consecutive samples
#
# needless to say, you should not jump directly to equation (13-117)
# rather try first to understand the entire thought process based on calculus
# identities, like derivative of the arctan function or derivatives of ratios
#

#
# use the four quadrant arctan function for phase detect between a pair of
# IQ samples; then unwrap the phase and take its derivative to demodulate
#
def fmDemodArctan(I, Q, prev_phase = 0.0):
#
# the default prev_phase phase is assumed to be zero, however
# take note in block processing it must be explicitly controlled

	# empty vector to store the demodulated samples
	fm_demod = np.empty(len(I))
	fm_demod2 = np.empty(len(I))

	# iterate through each of the I and Q pairs
	for k in range(len(I)):

		# if(k>0):
		# 	dQ = Q[k] - Q[k-1]
		# 	dI = I[k] - I[k-1]
		# else:
		# 	dQ = Q[k]
		# 	dI = I[k]
		# scaling = 1/(I[k]**2 + Q[k]**2)
		# fm_demod2[k] = (I[k]*dQ - Q[k]*dI)*scaling

		# use the atan2 function (four quadrant version) to detect angle between
		# the imaginary part (quadrature Q) and the real part (in-phase I)
		current_phase = math.atan2(Q[k], I[k])

		# we need to unwrap the angle obtained in radians through arctan2
		# to deal with the case when the change between consecutive angles
		# is greater than Pi radians (unwrap brings it back between -Pi to Pi)
		[prev_phase, current_phase] = np.unwrap([prev_phase, current_phase])

		# take the derivative of the phase
		fm_demod[k] = current_phase - prev_phase
		# print("Initial:",fm_demod[k])
		# print("New Implementation:", fm_demod2[k])


		# save the state of the current phase
		# to compute the next derivative
		prev_phase = current_phase

	# return both the demodulated samples as well as the last phase
	# (the last phase is needed to enable continuity for block processing)
	return fm_demod, prev_phase

def fmDemodArctan2(I, Q, prev_phase = 0.0):
#
# the default prev_phase phase is assumed to be zero, however
# take note in block processing it must be explicitly controlled

	# empty vector to store the demodulated samples
	fm_demod = np.empty(len(I))

	# iterate through each of the I and Q pairs
	for k in range(len(I)):

		# use the atan2 function (four quadrant version) to detect angle between
		# the imaginary part (quadrature Q) and the real part (in-phase I)
		# current_phase = math.atan2(Q[k], I[k])
		if(k>0):
			dQ = Q[k] - Q[k-1]
			dI = I[k] - I[k-1]
		else:
			dQ = Q[k]
			dI = I[k]
		scaling = 1/(I[k]**2 + Q[k]**2)
		fm_demod[k] = (I[k]*dQ - Q[k]*dI)*scaling
		# print(fm_demod[k])
		# we need to unwrap the angle obtained in radians through arctan2
		# to deal with the case when the change between consecutive angles
		# is greater than Pi radians (unwrap brings it back between -Pi to Pi)
		# [prev_phase, current_phase] = np.unwrap([prev_phase, current_phase])

		# take the derivative of the phase
		# fm_demod[k] = current_phase - prev_phase

		# save the state of the current phase
		# to compute the next derivative
		# prev_phase = current_phase

	# return both the demodulated samples as well as the last phase
	# (the last phase is needed to enable continuity for block processing)
	return fm_demod, prev_phase

# custom function for DFT that can be used by the PSD estimate
def DFT(x):

	# number of samples
	N = len(x)

	# frequency bins
	Xf = np.zeros(N, dtype='complex')

	# iterate through all frequency bins/samples
	for m in range(N):
		for k in range(N):
			Xf[m] += x[k] * cmath.exp(1j * 2 * math.pi * ((-k) * m) / N)

	# return the vector that holds the frequency bins
	return Xf

# custom function to estimate PSD based on the Bartlett method
# this is less accurate than the Welch method from matplotlib
# however, as the visual inspections confirm, the estimate gives
# the user a "reasonably good" view of the power spectrum
def estimatePSD(samples, NFFT, Fs):

	# rename the NFFT argument (notation consistent with matplotlib.psd)
	# to freq_bins (i.e., frequency bins for which we compute the spectrum)
	freq_bins = NFFT
	# frequency increment (or resolution)
	df = Fs/freq_bins

	# create the frequency vector to be used on the X axis
	# for plotting the PSD on the Y axis (only positive freq)
	freq = np.arange(0, Fs/2, df)

	# design the Hann window used to smoothen the discrete data in order
	# to reduce the spectral leakage after the Fourier transform
	hann = np.empty(freq_bins)
	for i in range(len(hann)):
		hann[i] = pow(math.sin(i*math.pi/freq_bins),2)

	# create an empty list where the PSD for each segment is computed
	psd_list = []

	# samples should be a multiple of frequency bins, so
	# the number of segments used for estimation is an integer
	# note: for this to work you must provide an argument for the
	# number of frequency bins not greater than the number of samples!
	no_segments = int(math.floor(len(samples)/float(freq_bins)))

	# iterate through all the segments
	for k in range(no_segments):

		# apply the hann window (using pointwise multiplication)
		# before computing the Fourier transform on a segment
		windowed_samples = samples[k*freq_bins:(k+1)*freq_bins] * hann

		# compute the Fourier transform using the built-in FFT from numpy
		Xf = np.fft.fft(windowed_samples, freq_bins)

		# note, you can check how MUCH slower is DFT vs FFT by replacing the
		# above function call with the one that is commented below
		#
		# Xf = DFT(windowed_samples)
		#
		# note, the slow impelementation of the Fourier transform is not as
		# critical when computing a static power spectra when troubleshooting
		# note: time permitting a custom FFT can be implemented

		# since input is real, we keep only the positive half of the spectrum
		# however, we will also add the signal energy of negative frequencies
		# to have a better a more accurate PSD estimate when plotting
		Xf = Xf[0:int(freq_bins/2)] # keep only positive freq bins
		psd_seg = 1/(Fs*freq_bins/2) * abs(Xf)**2 # compute signal power
		psd_seg = 2*psd_seg # add the energy from the negative freq bins

		# translate to the decibel (dB) scale
		for i in range(len(psd_seg)):
			psd_seg[i] = 10*math.log10(psd_seg[i])

		# append to the list where PSD for each segment is stored
		# in sequential order (first segment, followed by the second one, ...)
		psd_list.extend(psd_seg)

	# compute the estimate to be returned by the function through averaging
	psd_est = np.zeros(int(freq_bins/2))

	# iterate through all the frequency bins (positive freq only)
	# from all segments and average them (one bin at a time ...)
	for k in range(int(freq_bins/2)):
		# iterate through all the segments
		for l in range(no_segments):
			psd_est[k] += psd_list[k + l*int(freq_bins/2)]
		# compute the estimate for each bin
		psd_est[k] = psd_est[k] / no_segments

	# the frequency vector and PSD estimate
	return freq, psd_est


def FIR(cutoffFreq,samplingFreq,NTaps):
	NormCutoff = cutoffFreq / (samplingFreq/2)
	h = np.zeros((NTaps,),dtype=np.complex128)

	for i in range(NTaps):
		if i == (NTaps -1)/2:
			h[i] = NormCutoff
		else:
			h[i] = NormCutoff * (np.sin(np.pi*NormCutoff*(i-(NTaps-1)/2)))/(np.pi*NormCutoff*(i-(NTaps-1)/2))
		h[i] = h[i]*(np.sin((i*np.pi)/NTaps)) **2 
	
	return h

# def convolution (h ,x):

# 	#From 0 to kernelsize-2 (doesn't do last index)
# 	y = np.zeros(x.shape[0])

# 	#print(x.shape[0],h.shape[0])
	
# 	#From kernelsize-1 to len(x)-1
# 	for i in range(h.shape[0]-1,x.shape[0]-1):
# 		y[i] = 0
# 		j = i
# 		for k in range(h.shape[0]):
# 			y[i] += x[j]*h[k]
# 			j-=1
			

# 	for i in range(h.shape[0]):
# 		y[i] = 0
# 		k = 0
# 		for j in range(i,-1,-1):
# 			y[i] += x[j]*h[k]
# 			k+=1
	
# 	return y

def convolution(b, x, s):
    N = len(b)
    Yn = [0]*len(x)
	
    for n in range(len(x)):
		count=0
		for i in range(N):
			if (n-i>0):
				Yn[n] += b[i]*x[n-i]
    		else:
				Yn[n] += s[len(s)-1-count]*b[i]
				count+=1

    return Yn, s


def optimizedConvolution(x, h, state, downsamp, upsamp):
    y = [0]*(len(x)*upsamp//downsamp)
    s = len(state)-1
    count = 0

    for i in range(len(y)):
		count = 0
		for j in range(0:len(h):upsamp):
			if (0 <= downsamp*i-j and downsamp*i-j < len(x)*upsamp):
				y[i] += x[(downsamp*i-j)/upsamp]*h[j]
			else:
				y[i] += state[s-count]*h[j]
			count += 1
	return y


if __name__ == "__main__":

	# do nothing when this module is launched on its own
	pass