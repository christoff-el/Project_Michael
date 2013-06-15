#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "gemm.h"
#include "trsv_blk.h"
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

//+++++++++++++++++++++++++++
// Main
int main(int argc, char** argv)
{
	int N = 3;
	// Test TR Solver
	double *A, *B, *res;
	A = new double [N*N]; B = new double [N*N]; res = new double [N*N];
	// A=
	A[0] = 4;	A[1] = 2;	A[2] = 4;
	A[3] = 0; 	A[4] = 1;	A[5] = 3;
	A[6] = 0;	A[7] = 0;	A[8] = 3;
	// B=
	B[0] = 2; 	B[1] = 1;	B[2] = 3;
	B[3] = 2;	B[4] = 1;	B[5] = 3;
	B[6] = 0;	B[7] = 3;	B[8] = 3;
	cout << "A=" << endl; print_matrix(A, N, N, N, 1);
	cout << "B=" << endl; print_matrix(B, N, N, N, 1);
	// Solver
	trsv_blk('U', 'N', 'N', 'X', N, A, N, B, N, N);
	cout << "Solution= " << endl; print_matrix(B, N, N, N, 1);
#define blocksize 1
	hpc::dgemm(N,N,N,res,B,A,blocksize);
	cout << "X*A= " << endl; print_matrix(res, N, N, N, 1);
	delete[] A;
	delete[] B;
	delete[] res;
	return 0;
};