#ifndef LU_DECOMP_TRSV_BLK_H
#define LU_DECOMP_TRSV_BLK_H 1

#include "trsv_blk.tcc"

/***** Blocked Version of trsv_mat *****/

/***** Triangular Solver where b is a matrix *****
I.e. we solve Ax(:,i) = b(:,i) where i=0,...,#of cols of b -1

Implementation works bottom up, solving each system at the same time - should be most
cache friendly wrt accessing x.

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

order:
	CHARACTER*1.
	Specifies whether the system is of the form AX=B or XA=B
		if order = 'A' or 'a', then AX=B;
		if order = 'X' or 'x', then XA=B;
			(no change to transpose parameter is required).

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
trsv_blk(char uplo, char trans, char diag, char ordering, IndexType n,
		const A *a, IndexType lda,
				X *x, IndexType ldx, 
					IndexType m);

#endif	//LU_DECOMP_TRSV_BLK_H