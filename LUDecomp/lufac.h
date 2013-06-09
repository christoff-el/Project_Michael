#ifndef LU_DECOMP_H
#define LU_DECOMP_H 1
#include "../LUsolver/trsv_blk.h"
// LU factorization
// Function:	A <- (L|U)
// Decription:
//				- A is overwritten with L and U
//				- L is lower triangular with one on the diagonal
// 				- U is upper triangular
//				- ones and zeros are not stored
// size:		row or column size, matrix assumed quadratic
// incR:		increment to next row index
// incC:		increment to next column index			
template <typename IndexType, typename T>
void
lufac(T *A, IndexType size, IndexType incR, IndexType incC);

/*template <typename IndexType, typename T>
void
lufac_bl(T *A, IndexType size, IndexType incR, IndexType incC, IndexType blocksize)
{
	IndexType steps = size/blocksize; // first assumed to be a multiple
	for (IndexType i=0; i<steps-1; ++i)
	{
		lufac(A+i*blocksize*incR+i*blocksize*incC, blocksize, incR, incC);
		trsv_blk('L', 'N', 'U', blocksize, A+i*blocksize*incR+i*blocksize*incC, incR,
				A+i*blocksize*incR+(i+1)*blockisze*incC, incC,
				size-(i+1)*blockisze);// A12 tr. solver
		trsv_blk();// A21 tr. solver
		gemm_minus(A+(i+1)*blocksize*incR+(i+1)*blocksize*incC, size-(i+1)*blocksize, size-(i+1)*blocksize,
					A+(i+1)*blocksize*incR+i*blocksize*incC, size-(i+1)*blocksize, size-i*blocksize,
					A+i*blocksize*incR+(i+1)*blocksize*incC, size-i*blocksize*incC, size-(i+1)*blocksize);// A22 GEMM
	}
};*/

#include "lufac.tcc"
#endif