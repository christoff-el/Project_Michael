#ifndef DGEMM_TCC
#define DGEMM_TCC 1


void 
dgemm(int arow, int acol, int bcol, double *A, int astep, double *B, int bstep, double *C, int cstep)
{ 
  for(int i=0; i<arow; ++i)
     for(int j=0; j<bcol; ++j)
        for(int l=0; l<acol; ++l)
	   C[i*cstep+j]+=A[i*astep+l]*B[l*bstep+j];
} 


void
dgemmbl1(int arow, int acol, int bcol, double *A, int alength, double *B, int  blength, double *C, int bsize)
{
  int arowbl=arow / bsize + std::min(arow % bsize, 1);
  int acolbl=acol / bsize + std::min(acol % bsize, 1);
  int bcolbl=bcol / bsize + std::min(bcol % bsize, 1);
  
  int acurrow, acurcol, bcurcol;
  
  for(int i=0; i<arowbl; ++i)
    for(int j=0; j<bcolbl; ++j)
      for(int l=0; l<acolbl; ++l){
	acurrow=(i==arowbl-1) ? (arow-(arowbl-1)*bsize) : bsize;
	acurcol=(l==acolbl-1) ? (acol-(acolbl-1)*bsize) : bsize;
	bcurcol=(j==bcolbl-1) ? (bcol-(bcolbl-1)*bsize) : bsize;
	dgemm(acurrow, acurcol, bcurcol, &A[i*bsize*alength+l*bsize], alength, &B[l*bsize*blength+j*bsize], blength, &C[i*bsize*blength+j*bsize],blength);
      }
}


void
dgemmbl2(int arow, int acol, int bcol, double *A, double *B, double *C, int bsize)
{
  int arowbl=arow / bsize + std::min(arow % bsize, 1);
  int acolbl=acol / bsize + std::min(acol % bsize, 1);
  int bcolbl=bcol / bsize + std::min(bcol % bsize, 1);
  
  int acurrow, acurcol, bcurcol;
  
  for(int i=0; i<arowbl; ++i)
    for(int j=0; j<acolbl; ++j)
      for(int l=0; l<bcolbl; ++l){
	acurrow=(i==arowbl-1) ? (arow-(arowbl-1)*bsize) : bsize; 
	acurcol=(j==acolbl-1) ? (acol-(acolbl-1)*bsize) : bsize; 
	bcurcol=(l==bcolbl-1) ? (bcol-(bcolbl-1)*bsize) : bsize; 
	dgemm(acurrow, acurcol, bcurcol, &A[i*bsize*acol+j*bsize], acol, &B[j*bsize*bcol+l*bsize], bcol, &C[i*bsize*bcol+l*bsize],bcol);
      }
}

template <typename AA, typename BB, typename CC>
void 
dgemm_minus(int arow, int acol, int bcol, const AA *A, int astep, BB *B, int bstep, CC *C, int cstep){ 
  for(int i=0; i<arow; ++i)
     for(int j=0; j<bcol; ++j)
        for(int l=0; l<acol; ++l)
	   C[i*cstep+j]=C[i*cstep+j]-A[i*astep+l]*B[l*bstep+j];
} 

template <typename AA, typename BB, typename CC>
void
dgemmb_minus_l1(int arow, int acol, int bcol, AA *A, int alength, BB *B, int  blength, CC *C, int bsize)
{
 int arowbl=arow / bsize + std::min(arow % bsize, 1);
  int acolbl=acol / bsize + std::min(acol % bsize, 1);
  int browbl=acol / bsize + std::min(acol % bsize, 1);
  int bcolbl=bcol / bsize + std::min(bcol % bsize, 1);
  
  int acurrow, acurcol, bcurcol;
  
  for(int i=0; i<arowbl; ++i)
    for(int j=0; j<bcolbl; ++j)
      for(int l=0; l<browbl; ++l){
	acurrow=(i==arowbl-1) ? (arow-(arowbl-1)*bsize) : bsize;
	acurcol=(l==acolbl-1) ? (acol-(acolbl-1)*bsize) : bsize;
	bcurcol=(j==bcolbl-1) ? (bcol-(bcolbl-1)*bsize) : bsize;
	dgemm_minus(acurrow, acurcol, bcurcol, &A[i*bsize*alength+l*bsize], alength, &B[l*bsize*blength+j*bsize], blength, &C[i*bsize*blength+j*bsize],blength);
      }
}


#endif	//DGEMM_TCC