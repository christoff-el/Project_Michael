#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "lufac.h"
#include "timer.h"
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
	if (argc!=3)
	{
		cerr << endl << "Usage: #order of quadratic matrix #print matrices (1=yes)  " << endl;
		exit(1);
	}
	int N = 3;//atoi(argv[1]);
	int print = atoi(argv[2]);
	double *A, double *B;
	A = new double [N*N]; B = new double [N*N];

	// Fill with random values
	srand(0);
	for (int i=0; i<N; ++i)
	{
		for (int j=0; j<N, ++j)
		{
			A[i*N+j] = (double)rand()/RAND_MAX;
			B[i*N+j] = A[i*N+j];
		}
	}

	// Compute decompositions and measure time
	#define blocksize 1
	if (print==1)
	{
		cout << "A=" << endl; print_matrix(A, N, N, N, 1);
		cout << "B=" << endl; print_matrix(B, N, N, N, 1);
	}
	Timer timer;
	double t1, t2;
	timer.start();
	lufac(B, N, N, 1);
	timer.stop();
	t1 = timer.elapsed();

	timer.start();
	lufac_bl(A, N, N, 1, blocksize);
	timer.stop();
	t2 = timer.elapsed();

	if (print==1)
	{
		cout << "A=LU=" << endl; print_matrix(A, N, N, N, 1);
		cout << "B=LU=" << endl; print_matrix(B, N, N, N, 1);
	}
	double error = abs_diff(A, B, N, N);
	cout << setw(20) << "Error: " << setprecision(6) << error << endl;
	cout << setw(20) << "Unblocked LU: " << setprecision(6) << t1 << " sec" << endl;
	cout << setw(20) << "Blocked LU: " << setprecision(6) << t2 << " sec" << endl;
	cout << setw(20) << "Scaling factor: " << setprecision(6) << t1/t2 << endl;

	delete[] A;
	delete[] B;
	return 0;
}