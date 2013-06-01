#ifndef HPC_DGEMM_H
#define HPC_DGEMM_H 1

#include <iostream>
#include <iomanip>
#include <math.h>

namespace hpc {
// Simple GEMM
void
dgemm(int m, int n, int k, double *C, int lC, double *A, int lA, double *B, int lB) {
	for (int i=0; i<m; ++i) {
		for (int j=0; j<k; ++j) {
			for (int l=0; l<n; ++l) {
				C[i*lC+j] += A[i*lA+l]*B[l*lB+j];
			}
		}
	}
}
// Blocked GEMM
void
dgemm(int m, int n, int k, double *C, double *A, double *B, int blocksize) {
	int row_bl = m/blocksize;
	int column_bl_left = n/blocksize;
	int column_bl_right = k/blocksize;
	for (int i=0; i<row_bl; ++i) {
		for (int j=0; j<column_bl_left; ++j) {
			for (int l=0; l<column_bl_right; ++l) {
				hpc::dgemm(blocksize, blocksize, blocksize, C+i*k*blocksize+j*blocksize, k, A+i*n*blocksize+l*blocksize, n, B+l*k*blocksize+j*blocksize, k);
			}
		}
	}
}
// Blocked GEMM alternative order
void
dgemm2(int m, int n, int k, double *C, double *A, double *B, int blocksize) {
	int row_bl = m/blocksize;
	int column_bl_left = n/blocksize;
	int column_bl_right = k/blocksize;
	for (int i=0; i<row_bl; ++i) {
		for (int l=0; l<column_bl_right; ++l) {
			for (int j=0; j<column_bl_left; ++j) {
				hpc::dgemm(blocksize, blocksize, blocksize, C+i*k*blocksize+j*blocksize, k, A+i*n*blocksize+l*blocksize, n, B+l*k*blocksize+j*blocksize, k);
			}
		}
	}
}
// Initialize matrix
void
init_matrix(int m, int n, double *A) {
	for (int i=0; i<m; ++i) {
		for (int j=0; j<n; ++j) {
			A[i*n+j] = 0.;
		}
	}
}
// Difference Norm
double
norm_diff_matrix(int m, int n, double *A, double *B) {
	double norm = 0;
	for (int i=0; i<m; ++i) {
		for (int j=0; j<n; ++j) {
			if (fabs(A[i*n+j]-B[i*n+j])>norm) norm = fabs(A[i*n+j]-B[i*n+j]);
		}
	}
	return norm;
}
// Print matrix
void
print_matrix(int m, int n, double *A) {
	for (int i=0; i<m; ++i) {
		for (int j=0; j<n; ++j) {
			std::cout << std::setw(8) << A[i*n+j] << " ";
		}
		std::cout << std::endl;
	}
}
// Print submatrix
void
print_sub_matrix(int m, int n, int sub_m, int sub_n, int index_i, int index_j, double *A) {
	for (int i=0; i<sub_m; ++i) {
		for (int j=0; j<sub_n; ++j) {
			std::cout << std::setw(8) << A[index_i*n+index_j+i*n+j] << " ";
		}
		std::cout << std::endl;
	}
}
};

#endif