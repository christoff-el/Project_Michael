#ifndef LU_DECOMP_TRSV_BLK_H
#define LU_DECOMP_TRSV_BLK_H 1

#include "trsv_blk.tcc"

template <typename IndexType, typename A, typename X>
void
trsv_blk(char uplo, char trans, char diag, IndexType n,
		const A *a, IndexType lda,
				X *x, IndexType ldx, 
					IndexType m);

#endif	//LU_DECOMP_TRSV_BLK_H