#ifndef LU_DECOMP_UNBL_TCC
#define LU_DECOMP_UNBL_TCC 1
#include "_blas.h"

template <typename IndexType, typename T>
void
lufac(T *A, IndexType size, IndexType incR, IndexType incC)
{
	T alpha;
	const IndexType one = (IndexType) 1;
	const IndexType mone = (T) -1;
	for (IndexType i=0; i<size-1; ++i)
	{
		alpha = (T) 1/A[i*incR+i*incC];
		scal(&A[(i+one)*incR+i*incC], size-i-one, incR, alpha);
		ger(&A[(i+one)*incR+i*incC+one], size-i-one, size-i-one, incR, incC,
			&A[(i+one)*size+i], incR, &A[i*size+i+one], incC, mone);
	}
};

#endif