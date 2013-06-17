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

	#ifndef CACHESIZE
	#define CACHESIZE 32*1024
	#endif
	
	//Square blocksize such that A and B blocks fit in L1:
	int blkSize = sqrt((CACHESIZE / sizeof(double)) / 2);
	
	//Blocks in n and m dimensions:
	int blkCount_n = n / blkSize;
	int blkCount_m = m / blkSize;
	
	//Leftovers in n and m dimensions:
	int LO_n = n - blkCount_i * blkSize;
	int LO_m = m - blkCount_j * blkSize;
	
	
	for (int I=0; I<blkCount_n; ++I) {
	
		//Back-solve the triangular block on the diagonal:
		
	
		for (int J=0; J<I; J++) {
			
			
			
		}
	
	}
	
	
	
	
	
	
	
	
	
	
	for (int I=blkCount_i-1; I>=0; --I) {
	
		//Diagonal block:
		trsv_mat(uplo, trans, diag, ordering, blkSize, a[blkSize*I*lda + blkSize*I], lda, x[, IndexType ldx, 
					IndexType m);
		
		for (int J=1; J<blkCount_i; ++J) {
		
			
	
		}
	}


	//to keep everything working as before:
	trsv_mat(uplo, trans, diag, ordering, n, a, lda, x, ldx, m);



	//tbc..









}

#endif	//LU_DECOMP_TRSV_BLK_TCC