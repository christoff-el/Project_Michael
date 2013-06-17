#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstdlib>
#include "timer.h"
#include "../dgemm.tcc"
#include "utility.tcc"

/*The progtamm needs 5 parameters when you start it:
   1: number of rows of the first matrix; 
   2: number of cols of the first matrix;
   3: number of rows of the second matrix; 
   4: number of cols of the second matrix;
   5: block size for L1 cashing
*/

int main(int argc, char *argv[]){
    
   // Parameter Evaluation
   if(argc != 6){
      std::cerr << "Usage: " << argv[0] << " Sizes of matrices are not defined" << std::endl;
      exit(1);
   }
   
   int arow = atoi(argv[1]);
   int acol = atoi(argv[2]);
   int brow = atoi(argv[3]);
   int bcol = atoi(argv[4]);
   int bsize = atoi(argv[5]);
   

   // Fill array with random values in [0,1)
   assert( (arow>0) & (acol>0) & (brow>0) & (bcol>0)  & (acol==brow) );
   srand(2);
    
   std::cout << "Fill Array ... " << std::endl;
   double *A = new double[arow*acol];
   double *B = new double[brow*bcol];    
   double *C = new double[arow*bcol];    
   double *D = new double[arow*bcol];
   double *F = new double[arow*bcol]; 
   double *E = new double[arow*bcol]; 
   for(int i = 0; i < arow*acol; ++i){
     A[i] = (rand() / (double)RAND_MAX)*10.0;
   }
    
   for(int i = 0; i < brow*bcol; ++i){
     B[i] = (rand() / (double)RAND_MAX)*10.0;
   }
   
   for(int i = 0; i < arow*bcol; ++i){
      C[i]=(rand() / (double)RAND_MAX)*10.0;
      D[i]=C[i];
      F[i]=0;
      E[i]=0;
   }  
   std::cout << "... done" << std::endl;
   
   double time1, time2;
   
   Timer timer;
   
   timer.start();
   dgemm_minus(arow, acol, bcol, A, acol, B, bcol, C, bcol);
   timer.stop();
   time1=timer.elapsed();
   std::cout << "Without blocking : " << std::setprecision(8) << timer.elapsed() << " seconds" << std::endl;   
   
   timer.start();
   dgemmb_minus_l1(arow, acol, bcol, A, acol, B, bcol, D, bsize);
   timer.stop();
   time2=timer.elapsed();
   std::cout << "With blocking : " << std::setprecision(8) << timer.elapsed() << " seconds" << std::endl;   
   
   
   //MatrixPrint(A,arow,acol); std::cout<<std::endl; std::cout<<std::endl;
   std::cout<<"Scale:   " <<time1/time2<<std::endl;
   std::cout<<" C==D  " << check(C,D,arow, bcol, 0.1) <<std::endl;
   
delete []A;
delete []B;
delete []C;
delete []D;
delete []E;
delete []F;

}