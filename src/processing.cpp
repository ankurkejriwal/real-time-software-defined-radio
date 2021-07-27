#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include <cmath>
#include "support.h"


void stereo_block_processing(){
    
	int rf_Fs = 2400000;
	int rf_Fc = 100000;
	int rf_taps = 151;
	int rf_decim = 10;
	int audio_Fs = 240000;
	int audio_Fc = 16000;
	int audio_taps = 151; 
	int audio_decim = 5; 
	unsigned int block_size = 1024 * rf_decim * audio_decim * 2;
	int position = 0;
	unsigned int block_id = 0;
	std::vector<float> iq_data, i_data, q_data,iq_ds_coeff,i_state, q_state;
	std::vector<float> i_ds, q_ds, i_down, q_down;

    std::vector<float> state_recovery, state_extraction, stereo_pre_state;

    std::vector<float>pilot_coeff,bpf_recovery,recovery_coeff,stereo_recovered;


	// std::vector<float> i_ds, q_ds, i_down, q_down;

	std::vector<float> prev_phase,fm_demod; 
	std::vector<float> mono_coeff,audio_state,audio_block, audio_filter;
	std::vector<short int> audio_data;
	std::vector<std::complex<float>> Xf;
	std::vector<float> vector_index;
	std::vector<float> Xmag;
	std::vector<float> psd_est, freq;
	std::vector<float> psd_est1, freq1;

    std::vector<float> recovery_pll;

    pll_state statePLL;
    // data(statePLL);
	statePLL.integrator = 0.0;
	statePLL.phaseEst = 0.0;
	statePLL.feedbackI = 1.0;
	statePLL.feedbackQ = 0.0;
	statePLL.trigOffset = 0.0;
	statePLL.ncoLast = 1.0;

    std::vector<float> mixed;

    std::vector<float> stereo_coeff;

    std::vector<float> stereo_filt;


	//Sets some inital values
	prev_phase.resize(2,0.0);
	audio_state.resize(audio_taps-1,0.0);
	i_state.resize(rf_taps-1,0.0);
	q_state.resize(rf_taps-1,0.0);
	audio_state.resize(audio_taps-1,0.0);
	
	iq_data.resize(block_size);
	i_data.resize(block_size/2);
	q_data.resize(block_size/2);
	

    state_recovery.resize(rf_taps-1,0.0);
    state_extraction.resize(rf_taps-1,0.0);
    stereo_pre_state.resize(rf_taps-1,0.0);

    mixed.resize(5120,0.0);

    std::vector<float> left,right;

    std::vector<float> block99,block100;

    std::vector<float> audio_comb;


	impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, iq_ds_coeff);
	impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, mono_coeff);
	impulseResponseLPF(audio_Fs, audio_Fc, rf_taps, stereo_coeff);
	impulseResponseBPF(18.5e3,19.5e3,audio_Fs,151,pilot_coeff);
	impulseResponseBPF(22e3,54e3,audio_Fs,151,recovery_coeff);


	while(true)
	{

		//Now that data is in, need to seperate the I and Q data into seperate vectors
		//Maybe this index is wrong
		readStdinBlockData(block_size, block_id, iq_data);

		for(auto i = 0; i < (block_size)/2; i ++)
		{
			i_data[i] = iq_data[i*2];
			q_data[i] = iq_data[1+i*2];
		}

		//now generate the filter coefficents
		

		//Filtering 
		
        convolutionBlockIQ(i_ds, q_ds, i_data, q_data, iq_ds_coeff, i_state, q_state,rf_decim);


		//Demoadulate data
		fmDemodArctan(fm_demod,i_ds,q_ds,prev_phase);
		
		
		convolutionBlock(audio_block, fm_demod, mono_coeff, audio_state, audio_decim);
		//End of Audio Processing 


		convolutionBlock(bpf_recovery, fm_demod, pilot_coeff, state_recovery, 1);
        fmPLL(recovery_pll, bpf_recovery, 19e3, 240e3,2.0,0.0, 0.01,statePLL);

        
		convolutionBlock(stereo_recovered, fm_demod, recovery_coeff, state_extraction, 1);

        // pointwiseMult(stereo_recovered,recovery_pll,mixed);

        
        for (int i = 0; i <stereo_recovered.size();i++){
            mixed[i] = stereo_recovered[i]*recovery_pll[i];
        }
        

		convolutionBlock(stereo_filt, mixed, stereo_coeff, stereo_pre_state, 5);

        left.resize(audio_block.size(),0.0);
        right.resize(audio_block.size(),0.0);
          
        audio_comb.resize(2*audio_block.size());
        for(int i = 0 ; i<audio_block.size(); i++){
		    audio_comb[2*i] = (audio_block[i] + stereo_filt[i])/2; // left channel
		    audio_comb[2*i+1] = (audio_block[i] - stereo_filt[i])/2; // right channel
		
        }
	    //END OF COMBINER
		//Write blocks to stdout
		audio_data.resize(audio_comb.size());
		for(unsigned int l=0; l<audio_comb.size(); l++)
		{
			//If block is null or a weird value, set it to zero
			if(std::isnan(audio_comb[l])) 
			{
				audio_data[l] =0;
			}
			//Cast the value to short int for standard out and then scale it
			else 
			{
				audio_data[l] = static_cast<short int>(audio_comb[l] *16384);
			}
		}

		//Write to standard out
		fwrite(&audio_data[0],sizeof(short int),audio_data.size(), stdout);

		//Fill the elements with zero to ensure when filtering happens no weird numbers are added
		std::fill(i_ds.begin(), i_ds.end(), 0);
		std::fill(q_ds.begin(), q_ds.end(), 0);
		std::fill(audio_block.begin(), audio_block.end(), 0);
        std::fill(stereo_filt.begin(), stereo_filt.end(), 0);
        mixed.clear();

	
		//Will keep iterating the block till there is nothing coming in from standard in 
		if((std::cin.rdstate()) != 0)
		{
			break;
		}
	}
}


void mono_block_processing(){
	const int rf_Fs = 2400000;
	const int rf_Fc = 100e3;
	const int rf_taps = 151;
	const int rf_decim = 10;
	const int audio_Fs = 240e3;
	const int audio_Fc = 16e3;
	const int audio_taps = 151;
	const int audio_decim = 5;
	unsigned int block_size = 1024 * rf_decim * audio_decim * 2;

	unsigned int block = 0;

	std::vector<float> block_data,i_data, q_data,i_ds,q_ds;
	std::vector<float> i_state,q_state,state,rf_coeff,audio_coeff,audio_ds;
	
	std::vector<short int> audio_out;
	std::vector<float> fm_demod,phase;

	phase.resize(2,0.0);
	state.resize(audio_taps-1,0.0);
	i_state.resize(rf_taps-1,0.0);
	q_state.resize(rf_taps-1,0.0);
	
	block_data.resize(block_size);
	i_data.resize(block_size/2);
	q_data.resize(block_size/2);

	int position = 0;

	impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, rf_coeff); // filter coefficients for rf
	impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, audio_coeff); // filter coefficients for audio
	while (true){
		//  ----- Front End ----- //
		//Read and split stdin data 
		readStdinBlockData(block_size,block,block_data);
		vectorSplit(block_data, i_data, q_data);

		//Convolve both I and Q samples in the same loop
		convolutionBlockIQ(i_ds, q_ds, i_data, q_data, rf_coeff, i_state, q_state, rf_decim);
		
	
		//Demodulate Block Data
		fmDemodArctan(fm_demod,i_ds,q_ds,phase);
		//  ------------------- //


		//Audio Filtering and Decimation
		convolutionBlock(audio_ds, fm_demod, audio_coeff,state,audio_decim);
	
		audio_out.resize(audio_ds.size());

		for(unsigned int i=0; i<audio_ds.size(); i++){
			if(std::isnan(audio_ds[i])) {
				audio_out[i] =0;
			} else {
				audio_out[i] = static_cast<short int>(audio_ds[i] *16384);
			}
		}

		fwrite(&audio_out[0],sizeof(short int),audio_out.size(), stdout);

		std::fill(i_ds.begin(), i_ds.end(), 0);
		std::fill(q_ds.begin(), q_ds.end(), 0);
		std::fill(audio_ds.begin(), audio_ds.end(), 0);
		
		if((std::cin.rdstate()) != 0){
			break;
		}else{
			block++;
		}
	}
}