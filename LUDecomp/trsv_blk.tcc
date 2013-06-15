#ifndef LU_DECOMP_TRSV_BLK_TCC
#define LU_DECOMP_TRSV_BLK_TCC 1

#include "trsv_mat.h"

template <typename IndexType, typename A, typename X>
void
trsv_blk(char uplo, char trans, char diag, char ordering, IndexType n,
		const A *a, IndexType lda,
				X *x, IndexType ldx, 
					IndexType m)
{

	/*#ifndef CACHESIZE
	#define CACHESIZE 32*1024
	#endif
	
	int blkSize = sqrt((CACHESIZE / sizeof(double)) / 2);
	
	int blkCount_i = n / blkSize;
	int blkCount_j = m / blkSize;
	
	int LO_i = n - blkCount_i * blkSize;
	int LO_j = m - blkCount_j * blkSize;
	
	for (int I=blkCount_i-1; I>=0; --I) {
	
		//Diagonal block:
		trsv_mat(uplo, trans, diag, ordering, blkSize, a[blkSize*I*lda + blkSize*I], lda, x[, IndexType ldx, 
					IndexType m);
		
		for (int J=1; J<blkCount_i; ++J) {
		
			
	
		}
	}*/


	//to keep everything working as before:
	trsv_mat(uplo, trans, diag, ordering, n, a, lda, x, ldx, m);



	//tbc..









}

#endif	//LU_DECOMP_TRSV_BLK_TCC