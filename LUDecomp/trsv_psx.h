#ifndef LU_DECOMP_TRSV_PSX_H
#define LU_DECOMP_TRSV_PSX_H 1

#include "trsv_psx.tcc"

/***** POSIX distributor for trsv_blk.h *****

Distributes workload of trsv_blk to multiple processors.

Processor count and cache size should be defined by:
#define NUMPROC
#define CACHESIZE

Defaults are: 
NUMPROC = 4
CACHESIZE = 32*1024

Assumes double precision.

***** Input parameters: *****

uplo:
	CHARACTER*1. 
	Specifies whether the matrix A is upper or lower triangular:
		if uplo = 'U' or 'u', then the matrix is upper triangular;
		if uplo = 'L' or 'l', then the matrix is low triangular.

trans:
	CHARACTER*1. 
	Specifies the systems of equations:
		if trans = 'N' or 'n', then A*x = b;
		if trans = 'T' or 't', then A'*x = b;
		if trans = 'C' or 'c', then oconjg(A')*x = b.		<--- NOT IMPLEMENTED

diag:
	CHARACTER*1. 
	Specifies whether the matrix A is unit triangular:
		if diag = 'U' or 'u' then the matrix is unit triangular;
		if diag = 'N' or 'n', then the matrix is not unit triangular.

n:
	INTEGER. 
	Specifies the order of the matrix A. The value of n must be at least zero.

a:
	REAL for strsv
	DOUBLE PRECISION for dtrsv
	COMPLEX for ctrsv
	DOUBLE COMPLEX for ztrsv
	Array, DIMENSION (lda,n). 
		Before entry with uplo = 'U' or 'u', the leading n-by-n upper triangular part of 
		the array a must contain the upper triangular matrix and the strictly lower 
		triangular part of a is not referenced. --> It contains junk?
		Before entry with uplo = 'L' or 'l', the leading n-by-n lower triangular part of 
		the array a must contain the lower triangular matrix and the strictly upper 
		triangular part of a is not referenced.
		When diag = 'U' or 'u', the diagonal elements of a are not referenced either, but are assumed to be unity.

lda:
	INTEGER. 
	Specifies the leading dimension of a as declared in the calling (sub)program. 
	The value of lda must be at least max(1, n).

x:
	REAL for strsv
	DOUBLE PRECISION for dtrsv
	COMPLEX for ctrsv
	DOUBLE COMPLEX for ztrsv
	Array, DIMENSION at least (1 + (n - 1)*abs(incx)). 
		Before entry, the incremented array x must contain the n-element right-hand 
		side vector b.

ldx:
	INTEGER
	Specified the leading dimension of x.
	The value of ldx must be at least max(1,m).
	
m:
	INTEGER
	Specifies the number of columns of x. Must be at least zero.
		
ASSUMES ROW-MAJOR STORAGE

****************************************/

template <typename IndexType, typename A, typename X>
void
trsv_psx(char uplo, char trans, char diag, IndexType n,
		const A *a, IndexType lda,
				X *x, IndexType ldx, 
					IndexType m);

#endif	//LU_DECOMP_TRSV_PSX_H