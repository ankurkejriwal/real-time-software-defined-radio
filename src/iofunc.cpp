/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "iofunc.h"

// some basic functions for printing information from vectors
// or to read from/write to binary files in 32-bit float format
void printRealVector(const std::vector<float> &x)
{
	std::cout << "Printing float vector of size " << x.size() << "\n";
	for (auto i = 0; i < x.size(); i++)
		std::cout << x[i] << " ";
	std::cout << "\n";
}

void printComplexVector(const std::vector<std::complex<float>> &X)
{
	std::cout << "Printing complex vector of size " << X.size() << "\n";
	for (auto i = 0; i < X.size(); i++)
		std::cout << X[i] << " ";
	std::cout << "\n";
}

// assumes data in the raw binary file is in 32-bit float format
void readBinData(const std::string in_fname, std::vector<float> &bin_data)
{
	std::ifstream fdin(in_fname, std::ios::binary);
	if(!fdin) {
		std::cout << "File " << in_fname << " not found ... exiting\n";
		exit(1);
	} else {
		std::cout << "Reading raw binary from \"" << in_fname << "\"\n";
	}
	fdin.seekg(0, std::ios::end);
	const unsigned int num_samples = fdin.tellg() / sizeof(uint8_t);

	bin_data.resize(num_samples);
	fdin.seekg(0, std::ios::beg);
	fdin.read(reinterpret_cast<char*>(&bin_data[0]), num_samples*sizeof(uint8_t));
	fdin.close();
}

void readBinDataUint8(const std::string in_fname, std::vector<uint8_t> &bin_data)
{
	std::ifstream fdin(in_fname, std::ios::binary);
	if(!fdin) {
		std::cout << "File " << in_fname << " not found ... exiting\n";
		exit(1);
	} else {
		std::cout << "Reading raw binary from \"" << in_fname << "\"\n";
	}
	fdin.seekg(0, std::ios::end);
	const unsigned int num_samples = fdin.tellg() / sizeof(uint8_t);

	bin_data.resize(num_samples);
	fdin.seekg(0, std::ios::beg);
	fdin.read(reinterpret_cast<char*>(&bin_data[0]), num_samples*sizeof(uint8_t));
	fdin.close();
}

//function to read std in block data
void readStdinBlockData(unsigned int num_samples, unsigned int block_id, std::vector<float> &block_data){
	std::vector<char> raw_data(num_samples);
	std::cin.read(reinterpret_cast<char*>(&raw_data[0]), num_samples*sizeof(char));
	for(unsigned int k=0; k < num_samples; k++){
		block_data[k] = float(((unsigned char)raw_data[k]-128)/128.0);
	}
}

void writeStdoutBlockData(unsigned int num_samples, std::vector<float> &block_data){
	std::vector<short int> audio_data(num_samples);
	for(unsigned int k=0; k<block_data.size(); k++){
		if(std::isnan(block_data[k])) audio_data[k]=0;
		else audio_data[k] = static_cast<short int>(block_data[k]*16384);
	}

	fwrite(&audio_data[0], sizeof(short int), audio_data.size(), stdout);
}

// assumes data in the raw binary file is 32-bit float format
void writeBinData(const std::string out_fname, const std::vector<float> &bin_data)
{
	std::cout << "Writing raw binary to \"" << out_fname << "\"\n";
	std::ofstream fdout(out_fname, std::ios::binary);
	for (auto i=0; i<bin_data.size(); i++) {
		fdout.write(reinterpret_cast<const char*>(&bin_data[i]),\
								sizeof(bin_data[i]));
	}
	fdout.close();
}
