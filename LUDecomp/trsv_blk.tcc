#ifndef LU_DECOMP_TRSV_BLK_TCC
#define LU_DECOMP_TRSV_BLK_TCC 1

#include "trsv_mat.h"
#include "dgemm.h"

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
	
	//Square blocksize such that A, B and tmp blocks fit in L1:
	int blkSize = sqrt((CACHESIZE / sizeof(double)) / 3);
	
	//Blocks in n and m dimensions:
	int blkCount_n = n / blkSize;
	int blkCount_m = m / blkSize;
	
	//Leftovers in n and m dimensions:
	int LO_n = n - blkCount_n * blkSize;
	int LO_m = m - blkCount_m * blkSize;
	
	//For ease, cut out the leftovers from the top+left of A, and top of X:
	IndexType start_a = LO_n*lda + LO_n;
	IndexType start_x = LO_n*ldx;
	
	//For each block row of A:
	for (int I=blkCount_n; I>0; --I) {		
		
		//1. Main blocking bit:
		
		//For each block in the m dimension..
		for (int K=0; K<blkCount_m; ++K) {	
		
			//1. a. For each full block that needs updating with dgemm.. 
			for (int J=I; J<blkCount_n; ++J) {
			
				dgemm_minus(blkSize, blkSize, blkSize, 
								&a[start_a + lda*I*blkSize + J*blkSize], lda,
								&x[start_x + ldx*J*blkSize + K*blkSize], ldx,
								&x[start_x + ldx*I*blkSize + K*blkSize], ldx);
			
			}
		
			//1. b. Back-solve the triangular block on the diagonal.
			trsv_mat(uplo, trans, diag, ordering, blkSize,
						&a[start_a + lda*I*blkSize + I*blkSize],lda,
							&x[start_x + ldx*I*blkSize + K*blkSize],ldx,blkSize);
			
		}
	
	}
	
	//2. Leftovers in m dimension (currently done unblocked):

	trsv_mat(uplo, trans, diag, ordering,n,
				a,lda,
					&x[blkCount_m*blkSize],ldx,LO_m);
	
	//We are now left with the top LO_n rows of X unsolved (minus the last LO_m columns of these!).
	//To do these, we change the block dimensions to be short and wide:
	
	//or, at least, when all is working as it should first.. for now just botch it:
	for (int i=LO_n-1; i>=0; --i) {
		for (int j=0; j<m-LO_m; ++j) {
			for (int k=i+1; k<n; ++k) {
			
				x[i*ldx + j] -= a[i*lda + k] * x[k*ldx + j];
				
			}
			
			x[i*ldx + j] /= a[i*lda + i];
			
		}
	}
	
	/*int blkSize_i = LO_n;
	int blkSize_j = ((CACHESIZE/sizeof(double))/3)/blkSize_i;
	
	int LO_j = 
	
	*/

}

#endif	//LU_DECOMP_TRSV_BLK_TCC