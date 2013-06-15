#ifndef LU_DECOMP_TRSV_BLK_TCC
#define LU_DECOMP_TRSV_BLK_TCC 1

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
	
	int blkSize = sqrt((CACHESIZE / sizeof(double)) / 2)
















}

#endif	//LU_DECOMP_TRSV_BLK_TCC