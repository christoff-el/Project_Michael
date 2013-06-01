#ifndef LU_DECOMP_UNBL_H
#define LU_DECOMP_UNBL_H 1

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

#include "lufac.tcc"
#endif