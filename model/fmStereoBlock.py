#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#


import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math
from fmPll import fmPll

# use "custom" fmDemodArctan
from fmSupportLib import fmDemodArctan

def customlfilter(b, a, x, state):
    a = 1.0
    N = len(b)
    Yn = [0]*len(x)
    
    for n in range(len(x)):
    	for i in range(N):
    		count =0
    		if (n-i >=0 and n-i < len(x)):
    			Yn[n] += b[i]*x[n-i]
    		else:
    			Yn[n] += b[i]*state[len(state)-1 -count]
    			count+=1
    state = x[-(len(b)-1):]
    return Yn, state
def bandpassFilter(Fb, Fe, Fs, num_taps):
    norm_center = ((Fe+Fb)/2)/(Fs/2)
    norm_pass = (Fe-Fb)/(Fs/2)

    h = [0] * (num_taps-1) 
    for i in range(num_taps-1):
        if (i == (num_taps - 1)/2):
            h[i] = norm_pass
        else:
            h[i] = norm_pass * ((math.sin(math.pi*(norm_pass/2)*(i-(num_taps-1)/2))) / (math.pi*(norm_pass/2)*(i-(num_taps-1)/2)))
    
        h[i] = h[i] * math.cos(i*math.pi*norm_center)
        h[i] = h[i] * pow(math.sin((i*math.pi)/num_taps), 2)

    return h


# custom lowpass filter
def lowPassFilter(Fc, Fs, Ntaps):
	norm_Fc = Fc/(Fs/2)
	h = [0]*(Ntaps-1)

	for i in range(Ntaps-1):
		if i == (Ntaps-1)/2:
			h[i] = norm_Fc
		else:
			h[i] = norm_Fc*(math.sin(math.pi*norm_Fc*(i-(Ntaps-1)/2))/(math.pi*norm_Fc*(i-(Ntaps-1)/2)))

		h[i] = h[i]*(math.sin(i*math.pi/Ntaps))**2

	return h


def plotSpectrum(x, Fs, type = 'FFT'):

	n = len(x)             # length of the signal
	df = Fs/n              # frequency increment (width of freq bin)

	# compute Fourier transform, its magnitude and normalize it before plotting
	if type == 'FFT':
		Xfreq = np.fft.fft(x)
	#XMag = abs(Xfreq)/n
		
	if type == 'DFT':
		Xfreq = DFT(x)
    

	# Note: because x is real, we keep only the positive half of the spectrum
	# Note also: half of the energy is in the negative half (not plotted)
	XMag = abs(Xfreq)/n
	XMag = XMag[0:int(n/2)]

	# freq vector up to Nyquist freq (half of the sample rate)
	freq = np.arange(0, Fs/2, df)

	fig, ax = plt.subplots()
	ax.plot(freq, XMag)
	ax.set(xlabel='Frequency (Hz)', ylabel='Magnitude',
		title='Frequency domain plot')
	# fig.savefig("freq.png")
	plt.show()

def plotTime(x, time):

    fig, ax = plt.subplots()
    ax.plot(time, x)
    ax.set(xlabel='Time (sec)', ylabel='Amplitude',
            title='Time domain plot')
    # fig.savefig("time.png")
    plt.show()


rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151
rf_decim = 10

audio_Fs = 48e3
audio_Fc = 16e3
audio_taps = 151
audio_decim = 5

if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is assumed to be in 8-bits unsigned (and interleaved)
	# in_fname = "../data/test1.raw"
	in_fname = "radio.raw"
	raw_data = np.fromfile(in_fname, dtype='uint8')
	# IQ data is normalized between -1 and +1
	iq_data = (raw_data - 128.0)/128.0
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	# coefficients for the front-end low-pass filter
	rf_coeff = lowPassFilter(rf_Fc, rf_Fs, rf_taps)
    # coefficients for the filter to extract the mono audio
	audio_coeff = lowPassFilter(audio_Fc, audio_Fs, audio_taps)

	# set up the drawing
	fig, (ax0, ax1, ax2) = plt.subplots(nrows=3)
	fig.subplots_adjust(hspace = 1.0)

	# select a block_size that is a multiple of KB
	# and a multiple of decimation factors
	block_size = 1024 * rf_decim * audio_decim * 2
	block_count = 0

	# states needed for continuity in block processing
	state_i_lpf_100k = np.zeros(rf_taps-1)
	state_q_lpf_100k = np.zeros(rf_taps-1)
	state_lpf_16k = np.zeros(audio_taps-1)
	state_phase = 0
	state_pilot = np.zeros(rf_taps-1)
	state_stereo = np.zeros(rf_taps-1)
	state_stereo_lpf = np.zeros(rf_taps-1) 
	integrator = 0.0
	phaseEst = 0.0
	feedbackI = 1.0
	feedbackQ = 0.0
	ncoLast = 1.0
	trigOffset = 0
	state_pll = {
		'integrator' : 0.0,
		'phaseEst' : 0.0,
		'feedbackI' : 1.0,
		'feedbackQ' : 0.0,
		'ncoLast' : 1.0,
		'trigOffset' : 0
	}  
	# state_pll = (integrator, phaseEst, feedbackI, feedbackQ, trigOffset, ncoLast)

	# audio buffer that stores all the audio blocks
	audio_data = np.array([])
	recovery_pll_test = np.array([])
	# if the number of samples in the last block is less than the block size
	# it is fine to ignore the last few samples from the raw I/Q file
	recovery_state = [0.0, 0.0, 1.0, 0.0, 1.0, 0.0]
	while (block_count+1)*block_size < len(iq_data):
		# if you wish to have shorter runtimes while troubleshooting
		# you can control the above loop exit condition as you see fit

		print('Processing block ' + str(block_count))

        #Downsampled I and Q samples
		i_filt, state_i_lpf_100k = customlfilter(rf_coeff, 1.0, iq_data[(block_count)*block_size:(block_count+1)*block_size:2],state_i_lpf_100k)
		q_filt, state_q_lpf_100k = customlfilter(rf_coeff, 1.0, iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2], state_q_lpf_100k)

		# downsample the I/Q data from the FM channel
		i_ds = i_filt[::rf_decim]
		q_ds = q_filt[::rf_decim]

		# FM demodulator
		fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)

		# extract the mono audio data through filtering
		audio_filt, state_lpf_16k = customlfilter(audio_coeff, 1.0, fm_demod, state_lpf_16k)

		# downsample audio data
		audio_block = audio_filt[::audio_decim]

		# concatenate the most recently processed audio_block
		# to the previous blocks stored already in audio_data
		audio_data = np.concatenate((audio_data, audio_block))

        ##START STEREO PILOT RECOVERY
		pilot_coeff = bandpassFilter(18.5e3, 19.5e3, audio_Fs, audio_taps)
		pilot_recovery, state_pilot = customlfilter(pilot_coeff,1.0,fm_demod, state_pilot)
		recovery_pll, recovery_state = fmPll(pilot_recovery, 19e3, 240e3, recovery_state, 2)
        ##END OF STEREO PILOT RECOVERY

        
        ##START OF STEREO CHANNEL RECOVERY
		stereo_coeff = bandpassFilter(22e3, 54e3, audio_Fs, audio_taps)
		stereo_recovery, state_stereo = customlfilter(stereo_coeff,1.0,fm_demod, state_stereo)
        ##END OF STEREO CHANNEL RECOVERY

        ##START OF MIXING AND AUDIO FILTERING/DOWNSAMPLING
		mixed = np.multiply(recovery_pll[0:len(stereo_recovery):1],stereo_recovery)
		stereo_coeff = lowPassFilter(audio_Fc, audio_Fs, audio_taps)
		stereo_filt, state_stereo_lpf = customlfilter(stereo_coeff,1.0,mixed, state_stereo_lpf)
		stereo_data = stereo_filt[::5] 
        ## END OF MIXING AND DOWNSAMPLING

        ##BEGIN COMBINING 
		left, right = np.array(audio_block),np.array(audio_block)
		tmp_arr = np.zeros(10)
		if(block_count==0):
			tmp_arr = recovery_pll[-10:]
		if(block_count==1):
			tmp_arr = recovery_pll[0:10]
		recovery_pll_test = np.concatenate((recovery_pll_test, tmp_arr))
		for i in range(len(audio_block)):
			left[i] = (audio_block[i]+stereo_data[i])/2
			right[i] = (audio_block[i]-stereo_data[i])/2
        ##END COMBINING 

       

		# to save runtime select the range of blocks to log data
		# this includes both saving binary files as well plotting PSD
		# below we assume we want to plot for graphs for blocks 10 and 11
		if block_count >= 0 and block_count < 2:

            #  ##GRAPHS
			fig, (ax0,ax1,ax2,ax3,ax4,ax5,ax6) = plt.subplots(nrows=7)
			fig.subplots_adjust(hspace = 1.0)
			
			#PILOT FILTERING GRAPH
			ax0.psd(pilot_recovery, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
			ax0.set_ylabel('PSD (db/Hz)')
			ax0.set_title('Extracted Stereo Pilot Tone')
			
			#PLL MODULATED PILOT (38 Khz)
			ax1.psd(recovery_pll, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
			ax1.set_ylabel('PSD (db/Hz)')
			ax1.set_title('PLL Modulated Stereo Pilot Tone')

            #STEREO RECOVERY GRAP
			ax2.psd(stereo_recovery, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
			ax2.set_ylabel('PSD (db/Hz)')
			ax3.set_title('Stereo Recovery')
			
            
			

            #STEREO MIXED FILTERED AND DOWNSAMPLED
			ax3.psd(stereo_data, NFFT=512, Fs=(audio_Fs/audio_decim)/1e3)
			ax3.set_ylabel('PSD (db/Hz)')
			ax3.set_title('Downconversion and PLL mixed')

			plt.show()

		if block_count == 200:
			break

			# save figure to file
			fig.savefig("../data/fmAudio" + str(block_count) + ".png")

		block_count += 1
	step = 1/(len(recovery_pll_test)/2)
	t1 = np.arange(0, 2, step)
	plotTime(recovery_pll_test, t1)
	print('Finished processing the raw I/Q samples')

	# write audio data to a .wav file (assumes audio_data samples are -1 to +1)
	# wavfile.write("../data/fmAudio.wav", int(audio_Fs), np.int16((audio_data/2)*32767))
	left = np.int16((left)*32767)
	right = np.int16((right)*32767)
	stereo = np.array([left, right]).transpose()

	print("Writing Wav File")
	wavfile.write("../data/fmStereo.wav", int(48e3), stereo)

