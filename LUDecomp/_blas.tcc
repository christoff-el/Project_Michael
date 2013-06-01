#ifndef _BLAS_TCC
#define _BLAS_TCC 1

template <typename IndexType, typename T, typename SCAL>
void
scal(T *A, IndexType size, IndexType incA, SCAL &alpha)
{
	for (IndexType pos=0, i=0; i<size; ++i, pos+=incA)
	{
		A[pos] *= alpha;
	}
};

template <typename IndexType, typename T, typename VEC, typename SCAL>
void
ger(T *A, IndexType m, IndexType n, IndexType incm, IndexType incn, 
	VEC *a, IndexType inca, VEC *b, IndexType incb, SCAL &alpha)
{
	for (IndexType i_A=0, i_a=0, i=0; i<m; i_A+=incm, i_a+=inca, ++i)
	{
		for (IndexType j_A=0, j_b=0, j=0; j<n; j_A+=incn, j_b+=incb, ++j)
		{
			A[i_A+j_A] += alpha*a[i_a]*b[j_b];
		}
	}	
};

#endif