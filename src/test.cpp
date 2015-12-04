// LAPACK test code

#include <iostream>
#include <vector>
#include "../src/linalg/linalg.hpp"
#include <omp.h>

using namespace std;

extern "C" void dgetrf_ (int *m, int *n, double *a, int *lda, int *ipiv, int *info);
extern "C" void dgetrs_ (char *trans, int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

int main()
{
    char trans = 'N';
    int dim = 1024;    
    int nrhs = 1024*1024;
    int LDA = dim;
    int LDB = dim;
    int info;

    vector<double> a (dim * dim), b (dim * nrhs);

    for (int i = 0; i < dim; ++i)
    {
    	for (int j = 0; j < dim; ++j)
    	{
    		a [i * dim + j] = rand () % 100 / 100.;
    	}
    	for (int j = 0; j < nrhs; ++j)
    	{
    		b [i * dim + j] = rand () % 100 / 100.;
    	}
    }

    int ipiv[dim];

    int threads = 8;
    int dnrhs = nrhs / threads;
    int ddim = dim / threads;
    dgetrf_(&dim, &dim, &*a.begin(), &LDA, ipiv, &info);
    for (int i = 0; i < 1; ++i)
    {
    	#pragma omp parallel for
    	for (int j = 0; j < threads; ++j)
    	{
    		dgetrs_(&trans, &dim, &dnrhs, & *a.begin(), &LDA, ipiv, & *b.begin() + dim * j * ddim, &LDB, &info);
    	}
    }

    std::cout << "solution is:";    
    std::cout << "[" << b[0] << ", " << b[1] << ", " << "]" << std::endl;
    std::cout << "Info = " << info << std::endl; 

    return(0);
}