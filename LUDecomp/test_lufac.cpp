#include <iomanip>
#include <assert.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "lufac.h"
#include "timer.h"
using namespace std;

// Default blocksize
#ifndef BLOCKSIZE
#define BLOCKSIZE 4
#endif


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

// Norm
double
abs_diff(double *A, double *B, int m, int n)
{
	double res = 0;
	for (int i=0; i<m; ++i)
	{
		for (int j=0; j<n; ++j)
		{
			res += fabs(A[i*n+j]-B[i*n+j]); 
		}
	}
	return res;
}
//+++++++++++++++++++++++++++

// Main
int main (int argc, char** argv)
{
	// Read input
	if (argc!=4)
	{
		cerr << "Usage: #rows #columns #print matrices (1=yes)  " << endl;
		exit(1);
	}

	int N = atoi(argv[1]);
	int M = atoi(argv[2]);
	int print = atoi(argv[3]);
	double *A, *B, *L, *U, *C;

	// Check
	assert(N>1 && M>1 && BLOCKSIZE<=N && BLOCKSIZE<=M && BLOCKSIZE>0);

	// Cases: N>M, N<M, N=M
	int ld;
	A = new double [N*M]; B = new double [N*M];
	if (M>N)
	{
		L = new double [N*N]; U = new double [N*M];
		ld = N;
	}
	else
	{
		U = new double [M*M]; L = new double [N*M];
		ld = M;	
	}
	C = new double [N*M];

	/*A[0] = 1;	A[1] = 2;	A[2] = -1;	A[3] = 2;
	A[4] = 4;	A[5] = 3;	A[6] = 1;	A[7] = 3;
	A[8] = 2;	A[9] = 2;	A[10] = 3;	A[11] = 5; <--- simple test */

	// Fill with random values
	srand(0);
	for (int i=0; i<N; ++i)
	{
		for (int j=0; j<M; ++j)
		{
			A[i*M+j] = (double)rand()/RAND_MAX;
			B[i*M+j] = A[i*M+j];
			C[i*M+j] = A[i*M+j];
		}
	}


	// Compute decompositions and measure time

	if (print==1)
	{
		cout << "A=" << endl; print_matrix(A, N, M, M, 1);
		cout << "B=" << endl; print_matrix(B, N, M, M, 1);
	}

	// Unblocked
	Timer timer;
	double t1;
	timer.start();
	lufac(B, N, M, M);
	timer.stop();
	t1 = timer.elapsed();

	if (print==1)
	{
		cout << "(unbl) B=L|U=" << endl; print_matrix(B, N, M, M, 1);
	}

	for (int i=0; i<N; ++i)
	{
		for (int j=0; j<M; ++j)
		{
			if (j>i)
			{
				if (i<M) U[i*M+j] = B[i*M+j];
				if (j<N) L[i*ld+j] = 0.;
			}
			else if (j<i)
			{
				if (i<M) U[i*M+j] = 0.;
				if (j<N) L[i*ld+j] = B[i*M+j];
			}
			else if (j==i)
			{
				if (i<M) U[i*M+j] = B[i*M+j];
				if (j<N) L[i*ld+j] = 1.;	
			}
			B[i*M+j] = 0.;
		}
	}
	dgemm(N, ld, M, L, ld, U, M, B, M);
	double error_unbl = abs_diff(C, B, N, M);


	// Blocked
	double t2;
	timer.start();
	lufac_bl(A, N, M, M, BLOCKSIZE);
	timer.stop();
	t2 = timer.elapsed();
	
	if (print==1)
	{
		cout << "(bl) A=L|U=" << endl; print_matrix(A, N, M, M, 1);
	}
	
	for (int i=0; i<N; ++i)
	{
		for (int j=0; j<M; ++j)
		{
			if (j>i)
			{
				if (i<M) U[i*M+j] = A[i*M+j];
				if (j<N) L[i*ld+j] = 0.;
			}
			else if (j<i)
			{
				if (i<M) U[i*M+j] = 0.;
				if (j<N) L[i*ld+j] = A[i*M+j];
			}
			else if (j==i)
			{
				if (i<M) U[i*M+j] = A[i*M+j];
				if (j<N) L[i*ld+j] = 1.;	
			}
			A[i*M+j] = 0.;
		}
	}
	dgemm(N, ld, M, L, ld, U, M, A, M);
	double error_bl = abs_diff(C, A, N, M);

	// Output
	cout << endl;
	cout << setw(20) << "Error unblocked: " << setprecision(6) << error_unbl << endl;
	cout << setw(20) << "Error blocked: " << setprecision(6) << error_bl << endl;
	cout << setw(20) << "Time unblocked LU: " << setprecision(6) << t1 << " sec" << endl;
	cout << setw(20) << "Time blocked LU: " << setprecision(6) << t2 << " sec" << endl;
	cout << setw(20) << "Scaling factor: " << setprecision(6) << t1/t2 << endl;
	cout << endl;

	delete[] A;
	delete[] B;
	delete[] L;
	delete[] U;
	delete[] C;
	return 0;
}