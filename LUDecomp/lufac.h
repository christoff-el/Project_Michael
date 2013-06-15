#ifndef LU_DECOMP_H
#define LU_DECOMP_H 1
// LU factorization
// Function:	A <- (L|U)
// Decription:
//				- A is overwritten with L and U
//				- L is lower triangular with one on the diagonal
// 				- U is upper triangular
//				- ones and zeros are not stored
//
// Assumes row-major ordering
// m:			row size, m>1
// n:			column size, n>1
// lda:			increment to next row index			
template <typename IndexType, typename T>
void
lufac(T *A, IndexType m, IndexType n, IndexType lda);


// LU factorization, blocked
// Function:	A <- (L|U)
// Decription:
//				- A is overwritten with L and U
//				- L is lower triangular with one on the diagonal
// 				- U is upper triangular
//				- ones and zeros are not stored
//
// Assumes row-major ordering
// m:			row size, m>1
// n:			column size, n>1
// lda:			increment to next row index
// blocksize:	size of block, should be L1-cache friendly
template <typename IndexType, typename T>
void
lufac_bl(T *A, IndexType m, IndexType n, IndexType lda, IndexType blocksize);

#include "lufac.tcc"
#endif