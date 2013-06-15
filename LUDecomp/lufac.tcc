#ifndef LU_DECOMP_TCC
#define LU_DECOMP_TCC 1
#include "_blas.h"
#include "trsv_blk.h"
#include "../GEMM/dgemm.tcc"

#define MIN(a,b) a<b ? a : b

template <typename IndexType, typename T>
void
lufac(T *A, IndexType m, IndexType n, IndexType lda)
{
	T alpha;
	IndexType steps = MIN(m, n);
	const IndexType mone = (T) -1;
	IndexType i=0;
	for (; i<steps-1; ++i)
	{
		alpha = (T) 1/A[i*lda+i];
		scal(&A[(i+1)*lda+i], m-i-1, lda, alpha);
		ger(&A[(i+1)*lda+i+1], m-i-1, n-i-1, lda,
			&A[(i+1)*lda+i], lda, &A[i*lda+i+1], 1, mone);
	}
	// Over-determined system
	if (m>n)
	{
		alpha = (T) 1/A[i*lda+i];
		scal(&A[(i+1)*lda+i], m-i-1, lda, alpha);	
	}
};

template <typename IndexType, typename T>
void
lufac_bl(T *A, IndexType m, IndexType n, IndexType lda, IndexType blocksize)
{
	IndexType size = MIN(m,n);
	IndexType steps = size/blocksize;
	IndexType i = 0;
	for (;i<steps-1; ++i)
	{
		lufac(A+i*blocksize*lda+i*blocksize, blocksize, blocksize, lda); // A11
		trsv_blk('L', 'N', 'U', 'A', blocksize, A+i*blocksize*lda+i*blocksize, lda,
				A+i*blocksize*lda+(i+1)*blocksize, lda,
				n-(i+1)*blocksize);// A12 tr. solver
		trsv_blk('U', 'N', 'N', 'X', blocksize, A+i*blocksize*lda+i*blocksize, lda,
				A+(i+1)*blocksize*lda+i*blocksize, lda,
				m-(i+1)*blocksize);// A21 tr. solver
		dgemmb_minus_l1(m-(i+1)*blocksize, blocksize, n-(i+1)*blocksize,
						A+(i+1)*blocksize*lda+i*blocksize, lda,
						A+i*blocksize*lda+(i+1)*blocksize, lda,
						A+(i+1)*blocksize*lda+(i+1)*blocksize,
						blocksize);// A22 GEMM
	}
	// Remainder
	lufac(A+i*blocksize*lda+i*blocksize, m-i*blocksize, n-i*blocksize, lda);
};

#endif
