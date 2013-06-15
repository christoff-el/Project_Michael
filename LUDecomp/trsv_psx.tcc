#ifndef LU_DECOMP_TRSV_PSX_TCC
#define LU_DECOMP_TRSV_PSX_TCC 1

#include <thread>

#include "trsv_mat.h"

template <typename IndexType, typename A, typename X>
class SolveThread{
    
    char 		uplo, trans, diag, ordering;
    IndexType 	n, lda, ldx, m_sub;
    const A 	*a;
    X 			*x_begin;

public:
    SolveThread(char _uplo, char _trans, char _diag, char _ordering, IndexType _n, 
    				const A *_a, IndexType _lda, 
    					X *_x_begin, IndexType _ldx, 
    						IndexType _m_sub) 
    :	uplo(_uplo),
        trans(_trans),
        diag(_diag),
        ordering(_ordering),
        n(_n),
        lda(_lda),
        ldx(_ldx),
        m_sub(_m_sub),
        a(_a),
        x_begin(_x_begin)       
	{}
    
    void operator()(){        
        trsv_mat(uplo, trans, diag, ordering, n, a, lda, x_begin, ldx, m_sub);
    }
    
};

template <typename IndexType, typename A, typename X>
void
trsv_psx(char uplo, char trans, char diag, char ordering, IndexType n,
		const A *a, IndexType lda,
				X *x, IndexType ldx, 
					IndexType m)
{
	
	//Set defaults if need be:
	#ifndef NUMPROC
	#define NUMPROC 4
	#endif
	
	#ifndef CACHESIZE
	#define CACHESIZE 32*1024
	#endif
	
	//Initialise threads:
	std::thread threads[NUMPROC];
	
	/*	---Leaving out block optimisation for now, since m is unlikely to be so
			large that it will make any difference.
			
	//Optimum block size (such that a row of x and a row of A fit in L1):
	int blkSize = 0.5 * CACHESIZE / sizeof(double);

	//Number of blocks and leftovers:
	int blkCount = m/blkSize;
	int leftover_els = m - blkCount*blkSize;
	
	//Number of blocks given to each processor, and leftover blocks:
	int blkPerThread = blkCount / NUMPROC;
	int leftover_blks = blkCount - blkPerThread * NUMPROC;
	
	int elsPerCPU[NUMPROC], startEls[NUMPROC];
	
	for (int i=0; i<blkCount; ++i) {
	
		elsPerCPU[i - NUMPROC*((int)(i/NUMPROC))] += blkSize;
	
	}*/
	
	int elsPerThread = m / NUMPROC;
	int leftovers = m - elsPerThread * NUMPROC;
	
	int x_begin = 0;
	
	for (int i=0; i<NUMPROC; ++i) {
	
		int elCount = (i < leftovers) ? (elsPerThread+1) : (elsPerThread);
		threads[i] = std::thread(SolveThread<IndexType, A, X>
										(uplo, trans, diag, ordering, n, a, lda, &x[x_begin], ldx, elCount) );
		
		x_begin += elCount;
	
	}
	
	for (int i=0; i<NUMPROC; ++i) {
	
		threads[i].join();
		
	}

}

#endif	//LU_DECOMP_TRSV_PSX_TCC