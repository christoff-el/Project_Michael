#ifndef _BLAS_H
#define _BLAS_H 1

// Scaling
// Function:	A <- alpha*A
// size: 		total(!) size of array
// inca: 		increment to next element of A
template <typename IndexType, typename T, typename SCAL>
void
scal(T *A, IndexType size, IndexType incA, SCAL &alpha);


// Rank 1 update
// Function: 	A <- A + alpha* a*b^T
// m: 			row number
// n: 			column number
// incM:		increment to next row index
// incN:		increment to next column index
// ínca:		increment to next element of a
// incb:		============//================
template <typename IndexType, typename T, typename VEC, typename SCAL>
void
ger(T *A, IndexType m, IndexType n, IndexType incM, IndexType incN, 
	VEC *a, IndexType inca, VEC *b, IndexType incb, SCAL &alpha);


#include "_blas.tcc"
#endif