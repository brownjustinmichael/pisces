#include "utils/cuda/utils_cuda.hpp"
#include <vector>
#include <iostream>

int main (int argc, char const *argv[])
{
	std::vector <double> data (10), data_out (10);
	
	for (int i = 0; i < 10; ++i) {
		data [i] = i;
	}
	
	cuda::vector <double> data_device (10, &data [0]);
	
	data_device.copy_to_host (10, &data_out [0]);
	
	for (int i = 0; i < 10; ++i) {
		std::cout << data_out [i] << " ";
	}
	std::cout << "\n";
	return 0;
}