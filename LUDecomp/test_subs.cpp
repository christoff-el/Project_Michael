#include "_blas.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "lufac.h"
#include "gemm.h"
using namespace std;

//+++++++++++++++++++++++++++
// Display functions
void
print_matrix(double *A, int m, int n, int incm, int incn)
{
	for (int i=0, i_A=0; i<m; ++i, i_A+=incm)
	{
		for (int j=0, j_A=0; j<n; ++j, j_A+=incn)
		{
			cout << setw(8) << A[i_A+j_A] << " ";
		}
		cout << endl;
	}
}

void
print_vector(double *x, int size, int incx)
{
	for (int i=0, i_x=0; i<size; ++i, i_x+=incx)
	{
		cout << x[i_x] << endl;
	}
}
//+++++++++++++++++++++++++++


//+++++++++++++++++++++++++++
// Main
int main(int argc, char** argv)
{
	if (argc != 2)
	{
		cerr << "Usage " << argv[0] << ": " << "size of quadratic matrix" << endl;
		exit(1); 
	}	
	double *A, *x, *y;
	int N, M;
	N = atoi(argv[1]);
	M = N;
	srand(0);//time(NULL));

	//+++++++++++++++++++++++++++
	// Setup matrices and vectors
	A = new double [M*N]; x = new double [N]; y = new double [N];
	for (int i=0; i<M; ++i) 
	{
		x[i] = (double)rand()/RAND_MAX;//1.;
		y[i] = (double)rand()/RAND_MAX;//1.;
		for (int j=0; j<N; ++j)
		{
			A[i*N+j] = (double)rand()/RAND_MAX;//(i==j) ? 1. : 0.;
		}
	}
	//+++++++++++++++++++++++++++

	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

	//+++++++++++++++++++++++++++
	// Display
	cout << "The setup is:" << endl;
	cout << "A= " << endl;
	print_matrix(A, M, N, N, 1);
	cout << "x= " << endl;
	print_vector(x, N, 1);
	cout << "y= " << endl;
	print_vector(y, N, 1);
	//+++++++++++++++++++++++++++

	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;
	
	//+++++++++++++++++++++++++++
	// Scaling
	double alpha = 0.5;
	scal(A, M*N, 1, alpha);
	scal(x, N, 1, alpha);
	cout << "Test scaling constant: " << alpha << endl;
	cout << "Test scaling matrix A: " << endl;
	print_matrix(A, M, N, N, 1);
	cout << "Test scaling vector x: " << endl;
	print_vector(x, N, 1);
	//+++++++++++++++++++++++++++
	
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

	//+++++++++++++++++++++++++++
	// Rank update
	double beta = 1.;
	ger(A, M, N, N, 1, x, 1, y, 1, beta);
	cout << "Test rank update: " << endl;
	print_matrix(A, M, N, N, 1);
	//+++++++++++++++++++++++++++

	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;

	//+++++++++++++++++++++++++++
	// LU factorization
	/*A[0] = 1; A[1] = 2; A[2] = 3; // N = 3!!
	A[3] = 2; A[4] = -4; A[5] = 6;
	A[6] = 3; A[7] = -9; A[8] = -3;*/
	double *L = new double [N*N]; double *U = new double [N*N];
	double *result = new double [N*N];
	cout << "Test LU: " << endl;
	cout << "A= " << endl;
	print_matrix(A, M, N, N, 1);
	cout << endl;
	lufac(A, N, N, 1);
	for (int i=0; i<N; ++i)
	{
		for (int j=0; j<N; ++j)
		{
			if (j>=i)
			{
				U[i*N+j] = A[i*N+j];
				if (i==j) L[i*N+j] = 1.;
				else L[i*N+j] = 0.;
			}
			else
			{
				L[i*N+j] = A[i*N+j];
				U[i*N+j] = 0.;
			}
		}
	}
	cout << "L= " << endl;
	print_matrix(L, M, N, N, 1);
	cout << "U= " << endl;
	print_matrix(U, M, N, N, 1);
# define blocksize 1
	hpc::dgemm(N,N,N,result,L,U,blocksize);
	cout << "L*U= " << endl;
	print_matrix(result, M, N, N, 1);
	//+++++++++++++++++++++++++++

	// Free memory
	delete[] A;
	delete[] x;
	delete[] y;
	delete[] L;
	delete[] U;
	delete[] result;
	return 0;
}