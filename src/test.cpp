#include <stdio.h>

int main (int argc, char const *argv[])
{
	int i = 1;
	try {
		i = 100;
		throw 0;
	} catch (...) {
		printf ("i = %d\n", i);
	}
	return 0;
}