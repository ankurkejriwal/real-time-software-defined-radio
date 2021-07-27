#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include <cmath>
#include "support.h"
#include <functional>
#include <thread>
#include <condition_variable>
#include <queue>
#include <mutex>


void data(pll_state &samples){
	samples.integrator = 0.0;
	samples.phaseEst = 0.0;
	samples.feedbackI = 1.0;
	samples.feedbackQ = 0.0;
	samples.trigOffset = 0.0;
	samples.ncoLast = 1.0;
}

void front_end_producer(std::queue<std::vector<float>> &data_queue, std::mutex &queue_mutex, std::condition_variable &processing, bool &end_state,int mode){
	int rf_Fs;
	if (mode == 1){
		rf_Fs = 2500000;
	}
	else{
		rf_Fs = 2400000;
	}
	
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
	std::vector<float> iq_data, i_data, q_data,iq_ds_coeff;
	std::vector<float> i_ds, q_ds,i_state,q_state;


	std::vector<float> prev_phase, fm_demod; 
	std::vector<float> mono_coeff;


	//Sets some inital values
	prev_phase.resize(2,0.0);
	i_state.resize(rf_taps-1,0.0);
	q_state.resize(rf_taps-1,0.0);
	
	iq_data.resize(block_size);
	i_data.resize(block_size/2);
	q_data.resize(block_size/2);

	impulseResponseLPF(rf_Fs, rf_Fc, rf_taps, iq_ds_coeff);

	while(true)
	{

		//READ DATA 
		readStdinBlockData(block_size, block_id, iq_data);
		//SPLIT INTO I & Q SAMPLES
		vectorSplit(iq_data, i_data, q_data);

		//Filtering 
		convolutionBlockIQ(i_ds, q_ds, i_data, q_data, iq_ds_coeff, i_state, q_state,rf_decim);

		//Demoadulate data
		fmDemodArctan(fm_demod,i_ds,q_ds,prev_phase);

		//load queue with block data from front end
		std::unique_lock<std::mutex> queue_lock(queue_mutex);
		if(data_queue.size() == 5){
			processing.wait(queue_lock);
		}
		data_queue.push(fm_demod);
		// std::cerr<<"producer"<<std::endl;
		block_id++;
		queue_lock.unlock();
		processing.notify_one();

		std::fill(i_ds.begin(), i_ds.end(), 0);
		std::fill(q_ds.begin(), q_ds.end(), 0);

		if((std::cin.rdstate()) != 0)
		{
			end_state = true;
			break;
		}
	}
}

void stereo_consumer(std::queue<std::vector<float>> &data_queue, std::mutex &queue_mutex, std::condition_variable &processing, bool &end_state,int mode){
	int rf_taps = 151;
	int rf_decim = 10;
	int audio_up = 1;
	int audio_Fs = 240000;
	int audio_Fc = 16000;
	int audio_taps = 151; 
	const int audio_taps_1 = 151;
	int audio_decim = 5; 
	int pilot_b = 18.5e3;
	int pilot_e = 19.5e3;
	int stereo_b = 22e3;
	int stereo_e = 54e3;
	int stereo_taps = 151;
	int pilot_F = 19e3;

	if (mode) {
		audio_Fs = 6000000;
		audio_decim = 125; 
		audio_up = 24; 
		audio_taps = audio_taps*audio_up; 
	}

	unsigned int block_size = 1024 * rf_decim * audio_decim * 2;
	int position = 0;
	unsigned int block_id = 0;

	

	std::vector<float> iq_data, i_data, q_data,iq_ds_coeff,i_state, q_state;
	std::vector<float> i_ds, q_ds;

    std::vector<float> state_recovery, state_extraction, stereo_pre_state;
    std::vector<float>pilot_coeff,bpf_recovery,recovery_coeff,stereo_recovered;
    std::vector<float> recovery_pll;

	std::vector<float> prev_phase,fm_demod; 
	std::vector<float> mono_coeff,audio_state,audio_block, audio_filter;
	std::vector<short int> audio_data;



    pll_state statePLL;
    data(statePLL);

    std::vector<float> mixed;
	mixed.resize(5120,0.0);

    std::vector<float> stereo_coeff;

    std::vector<float> stereo_filt;

	//Sets some inital values
	prev_phase.resize(2,0.0);
	audio_state.resize(audio_taps-1,0.0);
	audio_state.resize(audio_taps-1,0.0);

    state_recovery.resize(rf_taps-1,0.0);
    state_extraction.resize(rf_taps-1,0.0);
    stereo_pre_state.resize(rf_taps-1,0.0);

    std::vector<float> left,right;
    std::vector<float> audio_comb;

	// LOW PASS COEFF FOR MONO AND STEREO
	impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, mono_coeff);
	impulseResponseLPF(audio_Fs, audio_Fc, audio_taps, stereo_coeff);
	// BAND PASS COEFF FOR PILOT AND STEREO RECOVERY
	impulseResponseBPF(pilot_b, pilot_e,audio_Fs/audio_up,stereo_taps,pilot_coeff);
	impulseResponseBPF(stereo_b, stereo_e,audio_Fs/audio_up,stereo_taps,recovery_coeff);


	while(true)
	{

		// UNLOCK WHEN DATA ON QUEUE
		std::unique_lock<std::mutex> queue_lock(queue_mutex);
		if(data_queue.empty()) {
			processing.wait(queue_lock);
		}
		fm_demod = data_queue.front();
		data_queue.pop();

		queue_lock.unlock();
		processing.notify_one();

		convolutionBlock(bpf_recovery, fm_demod, pilot_coeff, state_recovery, 1);
		fmPLL(recovery_pll, bpf_recovery, pilot_F, audio_Fs/audio_up,2.0,0.0, 0.01,statePLL);
		convolutionBlock(stereo_recovered, fm_demod, recovery_coeff, state_extraction, 1);
			

		if (mode) {
			mode1Convolution(audio_block,fm_demod,mono_coeff,audio_state,audio_decim,audio_up);
			for (int i = 0; i <stereo_recovered.size();i++){
            	mixed[i] = stereo_recovered[i]*recovery_pll[i];
        	}
			mode1Convolution(stereo_filt,mixed,stereo_coeff,stereo_pre_state,5,audio_up);
		} else {
			convolutionBlock(audio_block, fm_demod, mono_coeff, audio_state, audio_decim);
			for (int i = 0; i <stereo_recovered.size();i++){
            	mixed[i] = stereo_recovered[i]*recovery_pll[i]*2;
        	}
			convolutionBlock(stereo_filt, mixed, stereo_coeff, stereo_pre_state, 5);
		}

		// INTERLEAVING LEFT AND RIGHT CHANNELS
		audio_comb.resize(2*audio_block.size());
		for(int i = 0 ; i<audio_block.size(); i++){
			audio_comb[2*i] = (audio_block[i] + stereo_filt[i])/2; // left channel
			audio_comb[2*i+1] = (audio_block[i] - stereo_filt[i])/2; // right channel
		}
			
	    //END OF COMBINER
		audio_data.resize(audio_comb.size());
		for(unsigned int l=0; l<audio_comb.size(); l++)
		{
			if(std::isnan(audio_comb[l])) {
				audio_data[l] =0;
			} else {
				audio_data[l] = static_cast<short int>(audio_comb[l] *16384*audio_up);
			}
		}

		//WRITE DATA OUT
		fwrite(&audio_data[0],sizeof(short int),audio_data.size(), stdout);

		//FILL WITH 0s TO CLEAR OLD DATA
		std::fill(audio_block.begin(), audio_block.end(), 0);
        std::fill(stereo_filt.begin(), stereo_filt.end(), 0);
        mixed.clear();

		//Will keep iterating the block till there is nothing coming in from producer
		if(end_state && data_queue.empty())
		{
			break;
		}
	}
}