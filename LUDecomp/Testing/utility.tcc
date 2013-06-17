#ifndef UTILITY_TCC
#define UTILITY_TCC 1

#include <math.h>
#include "dgemm.tcc"

//Here gathered auxilary programms that can simplify test procedures


//Return "True" if matrices A and B are equal
bool check(double *A, double *B, int row, int col, double eps){
 bool l=true;
   
 for(int i=0; i<row; ++i)
    for(int j=0; j<col; ++j){
       l=(  (A[i*col+j]-B[i*col+j]) < eps);
       if(!l) {
	  std::cout<<"i=  "<<i <<"    j=" <<j << std::endl;
	  return l;
       }
    }
 return l;   
}


//Print the matrix A: dimension of the matrix A is arow x acol
template <typename ValueType>
void
MatrixPrint(ValueType *A, int arow, int acol){
   for(int i=0; i<arow; ++i){
      for(int j=0; j<acol; ++j){
         std::cout<<A[i*acol+j] <<"   ";
      }
   std::cout<<std::endl;
   }
} 


//Round the value a
int
Rounding(double a){
   if ( (a-round(a))>=0.5) {
      return round(a)+1;
   }else{
      return round(a);
   }
} 


// Matrix Operation: A=A-B
void
minus(int row, int col, double *A, double *B){
   for(int i=0; i<row; ++i){
      for(int j=0; j<col; ++j){
	 A[i*col+j]=A[i*col+j]- B[i*col+j];
      }
   }
}


// Matrix Operation: A=B;
void
equal(int row, int col, double *A, double *B){
   for(int i=0; i<row; ++i){
      for(int j=0; j<col; ++j){
	 A[i*col+j]= B[i*col+j];
      }
   }
}  

#endif