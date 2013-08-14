#include "utils.hpp"
#include <stdio.h>

int main (int argc, char const *argv[])
{
	
	float a [100], b [100];
	
	for (int i = 0; i < 100; ++i) {
		a [i] = i;
	}
	
	utils::scale (100, 2.0, a);
	
	for (int i = 0; i < 100; ++i) {
		printf ("Final %d: %f\n", i, a [i]);
	}
	
	return 0;
}