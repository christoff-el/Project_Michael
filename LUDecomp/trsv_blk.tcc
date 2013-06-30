#ifndef LU_DECOMP_TRSV_BLK_TCC
#define LU_DECOMP_TRSV_BLK_TCC 1

#include "trsv_mat.h"
#include "dgemm.h"
#include "./GEMMcache/gemmMainMinus.h"

template <typename IndexType, typename A, typename X>
void
trsv_blk(char uplo, char trans, char diag, char ordering, IndexType n,
		const A *a, IndexType lda,
				X *x, IndexType ldx, 
					IndexType m)
{

	#ifndef CACHESIZE
	#define CACHESIZE 32*1000
	#endif
	
	//Square blocksize such that A, B and C blocks fit in L1:
	int blkSize = 16;//sqrt((CACHESIZE / sizeof(double)) / 3);
	int blkSize2 = 16;//sqrt((CACHESIZE / sizeof(A)) / 3);
	
	
	//Blocks in n and m dimensions:
	int blkCount_n = n / blkSize;
	int blkCount_m = m / blkSize;
	
	//Leftovers in n and m dimensions:
	int LO_n = n - blkCount_n * blkSize;
	int LO_m = m - blkCount_m * blkSize;
	
	if (ordering == 'A' || ordering == 'a') {
	
		if (uplo == 'U' || uplo == 'u') {
	
			//For ease, cut out the leftovers from the top+left of A, and top of X:
			IndexType start_a = LO_n*lda + LO_n;
			IndexType start_x = LO_n*ldx;
	
			//For each block row of A:
			for (int I=blkCount_n; I>0; --I) {		
			
				//1. Main blocking bit:
		
				//For each block in the m dimension..
				for (int K=0, inxX2=start_x+ldx*(I-1)*blkSize; K<blkCount_m; ++K, inxX2+=blkSize) {	
		
					//1. a. For each full block that needs updating with dgemm.. 
					for (int J=I, inxA=start_a+lda*(I-1)+I*blkSize,
									inxX1=start_x+(ldx*I+K)*blkSize; 
															J<blkCount_n; 
																			++J, inxA+=blkSize,
																					inxX1+=ldx*blkSize) {
						
						/*dgemmb_minus_l1(blkSize, blkSize, blkSize, 
										&a[inxA], lda,
										&x[inxX1], ldx,
										&x[inxX2], blkSize2);*/
										
						
						gemmMainMinus(blkSize, &a[inxA], lda, &x[inxX1], ldx, &x[inxX2], ldx, blkSize2);
						
						
					}
		
					//1. b. Back-solve the triangular block on the diagonal.
					trsv_mat(uplo, trans, diag, ordering, blkSize,
								&a[start_a + lda*(I-1)*blkSize + (I-1)*blkSize],lda,
									&x[start_x + ldx*(I-1)*blkSize + K*blkSize],ldx,blkSize);
			
				}
	
			}
			//std::cout << "Leftovers" << std::endl;
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
	
		}
		else if (uplo == 'L' || uplo == 'l') {
	
	
			//For ease, cut out the leftovers from the bottom+right of A, and bottom of X:
		
	
			//For each block row of A:
			for (int I=1; I<=blkCount_n; ++I) {		
		
				//1. Main blocking bit:
		
				//For each block in the m dimension..
				for (int K=0; K<blkCount_m; ++K) {	
		
					//1. a. For each full block that needs updating with dgemm.. 
					for (int J=0; J<I-1; ++J) {
			
						dgemm_minus(blkSize, blkSize, blkSize, 
										&a[lda*(I-1)*blkSize + J*blkSize], lda,
										&x[ldx*J*blkSize + K*blkSize], ldx,
										&x[ldx*(I-1)*blkSize + K*blkSize], ldx);
			
					}
		
					//1. b. Back-solve the triangular block on the diagonal.
					trsv_mat(uplo, trans, diag, ordering, blkSize,
								&a[lda*(I-1)*blkSize + (I-1)*blkSize],lda,
									&x[ldx*(I-1)*blkSize + K*blkSize],ldx,blkSize);
			
				}
	
			}
			//std::cout << "Leftovers" << std::endl;
			//2. Leftovers in m dimension (currently done unblocked):

			trsv_mat(uplo, trans, diag, ordering,n,
						a,lda,
							&x[blkCount_m*blkSize],ldx,LO_m);
	
			//We are now left with the bottom LO_n rows of X unsolved (minus the last LO_m columns of these!).
			//To do these, we change the block dimensions to be short and wide:
	
			//or, at least, when all is working as it should first.. for now just botch it:
			for (int i=n-LO_n; i<n; ++i) {
				for (int j=0; j<m-LO_m; ++j) {
					for (int k=0; k<i; ++k) {
			
						x[i*ldx + j] -= a[i*lda + k] * x[k*ldx + j];
				
					}
			
					x[i*ldx + j] /= a[i*lda + i];
			
				}
			}
	
		}
	
	}
	
	else if (ordering == 'X' || ordering == 'x') {

		if (uplo == 'U' || uplo == 'u') {

			//For ease, cut out the leftovers from the top+left of A, and left of X:
			IndexType start_a = LO_n*lda + LO_n;
			IndexType start_x = LO_n;
			
			//std::cout<<blkSize<<" "<<start_a<<" "<<start_x<<" "<<std::endl;
			
			for (int i=0; i<n; ++i) {
				for (int j=0; j<n; ++j) {
					//std::cout<<a[lda*i+j]<<" ";
				}
				//std::cout<<std::endl;
			}
			//std::cout<<std::endl;
			
			for (int i=0; i<m; ++i) {
				for (int j=0; j<n; ++j) {
					//std::cout<<x[ldx*i+j]<<" ";
				}
				//std::cout<<std::endl;
			}
			//std::cout<<std::endl;
			
	
			//For each block col of A:
			for (int I=1; I<=blkCount_n; ++I) {		
		
				//1. Main blocking bit:
		
				//For each block in the m dimension..
				for (int K=0; K<blkCount_m; ++K) {	
		
					//1. a. For each full block that needs updating with dgemm.. 
					for (int J=0; J<I-1; ++J) {
			
						//std::cout<< "I=" << I << " K=" << K << " J=" << J << std::endl;
						//std::cout<< std::endl;
						dgemm_minus(blkSize, blkSize, blkSize, 
										&a[start_a + lda*J*blkSize + (I-1)*blkSize], lda,
										&x[start_x + ldx*J*blkSize + K*blkSize], ldx,
										&x[start_x + ldx*(I-1)*blkSize + K*blkSize], ldx);
			
					}
		
					//1. b. Back-solve the triangular block on the diagonal.
					trsv_mat(uplo, trans, diag, ordering, blkSize,
								&a[start_a + lda*(I-1)*blkSize + (I-1)*blkSize],lda,
									&x[start_x + ldx*K*blkSize + (I-1)*blkSize],ldx,blkSize);
			
				}
	
			}
			//std::cout << "Leftovers" << std::endl;
			//2. Leftovers in m dimension (currently done unblocked):

			trsv_mat(uplo, trans, diag, ordering, n,
						a,lda,
							&x[blkCount_m*blkSize*ldx],ldx,LO_m);
	
			//We are now left with the left LO_n columns of X unsolved (minus the last LO_m rows of these!).
			//To do these, we change the block dimensions to be short and wide:
	
			//or, at least, when all is working as it should first.. for now just botch it:
			for (int i=0; i<n-LO_m; ++i) {
				for (int j=0; j<LO_n; ++j) {
					for (int k=0; k<j; ++k) {
					
						x[i*ldx + j] -= a[k*lda + j] * x[i*ldx + k];
				
					}
			
					x[i*ldx + j] /= a[j*lda + j];
			
				}
			}
			
			for (int i=0; i<m; ++i) {
				for (int j=0; j<n; ++j) {
					//std::cout<<x[ldx*i+j]<<" ";
				}
				//std::cout<<std::endl;
			}
			//std::cout<<std::endl;
	
		}
	
	}

}

#endif	//LU_DECOMP_TRSV_BLK_TCC