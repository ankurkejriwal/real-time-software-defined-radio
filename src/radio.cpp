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
#include "threadfunc.h"
#include "processing.h"
	
int main(int argc, char* argv[])
{
	int mode = 0;
	if(argc < 2)
	{
		std::cerr << "Operating in mode 0" << std::endl;
	}
	else if (argc == 2)
	{
		mode = atoi(argv[1]);
		if(mode!=1)
		{
			std::cerr << "Wrong mode " <<mode <<std::endl;
			exit(1);
		}
	}
	else
	{
		std::cerr << "Usage " <<argv[0] <<std::endl;
	}

	std::queue<std::vector<float>> data_queue;
	std::mutex queue_mutex;
	std::condition_variable processing;
	bool end_state = false; 

	// multi threading stereo
	std::thread producer = std::thread(front_end_producer, std::ref(data_queue), std::ref(queue_mutex), std::ref(processing), std::ref(end_state),mode);
	std::thread consumer= std::thread(stereo_consumer, std::ref(data_queue), std::ref(queue_mutex), std::ref(processing), std::ref(end_state),mode);

	producer.join();
	consumer.join();

	// FOR TESTING BLOCK PROCESSING OF STEREO AND MONO IMP
	// stereo_block_processing();

	// mono_block_processing();

	return 0;
}

