#ifndef LU_DECOMP_TRSV_H
#define LU_DECOMP_TRSV_H 1

#include "trsv.tcc"

template <typename IndexType, typename A, typename X>
void
trsv(char uplo, char trans, char diag, IndexType n,
		const A *a, IndexType lda,
				X *x, IndexType incX);

#endif	//LU_DECOMP_TRSV_H