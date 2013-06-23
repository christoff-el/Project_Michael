#ifndef DGEMM_H
#define DGEMM_H 1

#include "dgemm.tcc"

// Description:
// 					arow  - rows number of matrix A
// 					acol  - columns number of matrix A
// 					bcol  - columns number of matrix B
// 					*A    - pointer to matrix A
// 					*B    - pointer to matrix B
// 					*C    - pointer to matrix C
// 					astep - 
// 					bstep - 

// Matrix multiplication without block 
// C = A * B
void 
dgemm(int arow, int acol, int bcol, double *A, int astep, double *B, int bstep, double *C, int cstep);


// Blocked dgemm, taking advantage of L1
void
dgemmbl1(int arow, int acol, int bcol, double *A, int alength, double *B, int  blength, double *C, int bsize);


// Blocked dgemm, taking advantage of L1 & L2
void
dgemmbl2(int arow, int acol, int bcol, double *A, double *B, double *C, int bsize);

   
// C = C - A*B     without block
template <typename AA, typename BB, typename CC>
void 
dgemm_minus(int arow, int acol, int bcol, const AA *A, int astep, BB *B, int bstep, CC *C, int cstep);



//Blocked dgemm_minus, taking advantage of L1;
void
dgemmb_minus_l1(int arow, int acol, int bcol, double *A, int alength, double *B, int  blength, double *C, int bsize);


#endif	//DGEMM_H