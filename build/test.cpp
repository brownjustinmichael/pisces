#include "mpi.h"
#include <vector>

int main(int argc, char *argv[])
{
	MPI::Init (argc, argv);

	std::vector <double> x (100), y (1000);
	for (int i = 0; i < 100000000; ++i)
	{
		MPI::COMM_WORLD.Gather (&x [0], 100, MPI::DOUBLE, &y [0], 100, MPI::DOUBLE, 0);
	}

	MPI::Finalize ();
	return 0;
}