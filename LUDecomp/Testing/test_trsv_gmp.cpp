#define PREC 1e-5
#define NUMPROC 2

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <gmpxx.h>

#include "../trsv.h"
#include "../trsv_blk.h"
#include "../trsv_psx.h"
#include "timer.h"

using namespace std;

mpf_class errCalc(int n, mpf_class *exact, mpf_class *calc, int incX, int prec);
mpf_class errCalc2d(int n, int m, mpf_class *exact, mpf_class *calc, int prec);
void testTRSV_gmp(int n, int prec);
void testTRSV_blk_gmp(int n, int m, int prec);
void testTRSV_PSX_gmp(int n, int m, int prec);


int main(int argc, char **argv) {

	if (argc != 4) {
		cerr << "Usage: " << argv[0] << " n m precision" << endl;
		return 1;
	}
	
	//System size:
	int n = atoi(argv[1]);
	int m = atoi(argv[2]);
	int prec = atoi(argv[3]);
	

	
	testTRSV_gmp(n,prec);
	testTRSV_blk_gmp(n,m,prec);
	testTRSV_PSX_gmp(n,m,prec);
	
	
	return 0;
	
}


void testTRSV_gmp(int n, int prec) {

	mpf_set_default_prec(prec);

	//Build A and x, then compute b:
	mpf_class *A_U = new mpf_class[n*n];
	mpf_class *A_L = new mpf_class[n*n];
	mpf_class *A_U_D = new mpf_class[n*n];	//=A_U with 1s on diagonal
	mpf_class *A_L_D = new mpf_class[n*n];	//=A_L with 1s on diagonal
	mpf_class *x = new mpf_class[n];
	mpf_class *b_U = new mpf_class[n];
	mpf_class *b_L = new mpf_class[n];
	mpf_class *b_U_T = new mpf_class[n];		//=b_L
	mpf_class *b_L_T = new mpf_class[n];		//=b_U
	mpf_class *b_U_D = new mpf_class[n];
	mpf_class *b_L_D = new mpf_class[n];
	mpf_class *b_U_T_D = new mpf_class[n];
	mpf_class *b_L_T_D = new mpf_class[n];
	mpf_class *b_U_M = new mpf_class[n];			//_M same as above, but backwards.
	mpf_class *b_L_M = new mpf_class[n];
	mpf_class *b_U_T_M = new mpf_class[n];		//=b_L_M
	mpf_class *b_L_T_M = new mpf_class[n];		//=b_U_M
	mpf_class *b_U_D_M = new mpf_class[n];
	mpf_class *b_L_D_M = new mpf_class[n];
	mpf_class *b_U_T_D_M = new mpf_class[n];
	mpf_class *b_L_T_D_M = new mpf_class[n];
	
	//Fill A, x with random numbers:
	srand(time(NULL));
    for (int i=0; i<n; ++i) {
        for (int j=i; j<n; ++j) {
        
        	A_U[j*n +i] = 0;
        	A_L[i*n +j] = 0;
        	
        	mpf_class tmp=rand()/(mpf_class)RAND_MAX;
            A_U[i*n +j] = tmp;
            A_L[j*n +i] = tmp;
            
            A_U_D[i*n +j] = A_U[i*n +j];
            A_L_D[j*n +i] = A_L[j*n +i];
            
        }
        
        A_U_D[i*n +i] = 1;
        A_L_D[i*n +i] = 1;
        
        x[i] = rand()/(mpf_class)RAND_MAX;
        
    }
    
    //b = Ax
    for (int i=0; i<n; ++i) {
    
    	b_U[i] = 0;
    	b_L[i] = 0;
    	
    	b_U_D[i] = 0;
    	b_L_D[i] = 0;
    	
    	b_U_M[n-1-i] = 0;
    	b_L_M[n-1-i] = 0;
    	
    	b_U_D_M[n-1-i] = 0;
    	b_L_D_M[n-1-i] = 0;
    	
    	for (int j=0; j<n; ++j) {
    	
    		b_U[i] += A_U[i*n +j] * x[j];
    		b_L[i] += A_L[i*n +j] * x[j];
    		
    		b_U_D[i] += A_U_D[i*n +j] * x[j];
    		b_L_D[i] += A_L_D[i*n +j] * x[j];
    		
    		b_U_M[n-1-i] += A_U[i*n +j] * x[j];
    		b_L_M[n-1-i] += A_L[i*n +j] * x[j];
    		
    		b_U_D_M[n-1-i] += A_U_D[i*n +j] * x[j];
    		b_L_D_M[n-1-i] += A_L_D[i*n +j] * x[j];
    		
    	}
    	
    	b_U_T[i] = b_L[i];
    	b_L_T[i] = b_U[i];
    	
    	b_U_T_D[i] = b_L_D[i];
    	b_L_T_D[i] = b_U_D[i];
    	
    	b_U_T_M[n-1-i] = b_L_M[n-1-i];
    	b_L_T_M[n-1-i] = b_U_M[n-1-i];
    	
    	b_U_T_D_M[n-1-i] = b_L_D_M[n-1-i];
    	b_L_T_D_M[n-1-i] = b_U_D_M[n-1-i];
    	
    }

	
	//Solve using the solver:
	int lda = n;

	//incX=1:
	//Upper:
	trsv('U', 'N', 'N', n, A_U, lda, b_U, 1);
	trsv('U', 'T', 'N', n, A_U, lda, b_U_T, 1);
	trsv('U', 'N', 'U', n, A_U_D, lda, b_U_D, 1);
	trsv('U', 'T', 'U', n, A_U_D, lda, b_U_T_D, 1);

	//Lower:
	trsv('L', 'N', 'N', n, A_L, lda, b_L, 1);
	trsv('L', 'T', 'N', n, A_L, lda, b_L_T, 1);
	trsv('L', 'N', 'U', n, A_L_D, lda, b_L_D, 1);
	trsv('L', 'T', 'U', n, A_L_D, lda, b_L_T_D, 1);
	
	//incX=-1
	//Upper:
	trsv('U', 'N', 'N', n, A_U, lda, b_U_M, -1);
	trsv('U', 'T', 'N', n, A_U, lda, b_U_T_M, -1);
	trsv('U', 'N', 'U', n, A_U_D, lda, b_U_D_M, -1);
	trsv('U', 'T', 'U', n, A_U_D, lda, b_U_T_D_M, -1);

	//Lower:
	trsv('L', 'N', 'N', n, A_L, lda, b_L_M, -1);
	trsv('L', 'T', 'N', n, A_L, lda, b_L_T_M, -1);
	trsv('L', 'N', 'U', n, A_L_D, lda, b_L_D_M, -1);
	trsv('L', 'T', 'U', n, A_L_D, lda, b_L_T_D_M, -1);


	//Output results:
	cout << "\nResults for trsv:\n" << endl;
	
	mpf_class err_U = errCalc(n,x,b_U,1,prec);
	mpf_class err_U_T = errCalc(n,x,b_U_T,1,prec);
	mpf_class err_U_D = errCalc(n,x,b_U_D,1,prec);
	mpf_class err_U_T_D = errCalc(n,x,b_U_T_D,1,prec);
	
	mpf_class err_L = errCalc(n,x,b_L,1,prec);
	mpf_class err_L_T = errCalc(n,x,b_L_T,1,prec);
	mpf_class err_L_D = errCalc(n,x,b_L_D,1,prec);
	mpf_class err_L_T_D = errCalc(n,x,b_L_T_D,1,prec);
	
	mpf_class err_U_M = errCalc(n,x,b_U_M,-1,prec);
	mpf_class err_U_T_M = errCalc(n,x,b_U_T_M,-1,prec);
	mpf_class err_U_D_M = errCalc(n,x,b_U_D_M,-1,prec);
	mpf_class err_U_T_D_M = errCalc(n,x,b_U_T_D_M,-1,prec);
	
	mpf_class err_L_M = errCalc(n,x,b_L_M,-1,prec);
	mpf_class err_L_T_M = errCalc(n,x,b_L_T_M,-1,prec);
	mpf_class err_L_D_M = errCalc(n,x,b_L_D_M,-1,prec);
	mpf_class err_L_T_D_M = errCalc(n,x,b_L_T_D_M,-1,prec);
	
	cout << "incX = 1:" << endl;
	cout << fixed << setprecision(10) << "Upper: " << err_U << "\t\t\t ----> " << ((err_U<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(10) << "Upper transpose: " << err_U_T << "\t\t ----> " << ((err_U_T<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(10) << "Upper unit: " << err_U_D << "\t\t ----> " << ((err_U_D<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(10) << "Upper unit transpose: " << err_U_T_D << "\t ----> " << ((err_U_T_D<PREC)?"PASS":"FAIL") << endl;
		
	cout << fixed << setprecision(10) << "Lower: " << err_L << "\t\t\t ----> " << ((err_L<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(10) << "Lower transpose: " << err_L_T << "\t\t ----> " << ((err_L_T<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(10) << "Lower unit: " << err_L_D << "\t\t ----> " << ((err_L_D<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(10) << "Lower unit transpose: " << err_L_T_D << "\t ----> " << ((err_L_T_D<PREC)?"PASS":"FAIL") << endl;

	cout << "\nincX = -1:" << endl;
	cout << fixed << setprecision(10) << "Upper: " << err_U_M << "\t\t\t ----> " << ((err_U_M<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(10) << "Upper transpose: " << err_U_T_M << "\t\t ----> " << ((err_U_T_M<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(10) << "Upper unit: " << err_U_D_M << "\t\t ----> " << ((err_U_D_M<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(10) << "Upper unit transpose: " << err_U_T_D_M << "\t ----> " << ((err_U_T_D_M<PREC)?"PASS":"FAIL") << endl;
		
	cout << fixed << setprecision(10) << "Lower: " << err_L_M << "\t\t\t ----> " << ((err_L_M<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(10) << "Lower transpose: " << err_L_T_M << "\t\t ----> " << ((err_L_T_M<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(10) << "Lower unit: " << err_L_D_M << "\t\t ----> " << ((err_L_D_M<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(10) << "Lower unit transpose: " << err_L_T_D_M << "\t ----> " << ((err_L_T_D_M<PREC)?"PASS":"FAIL") << endl;


	//Clean up:
	delete[] A_U;
	delete[] A_L;
	delete[] A_U_D;
	delete[] A_L_D;
	delete[] x;
	delete[] b_U;
	delete[] b_L;
	delete[] b_U_T;
	delete[] b_L_T;
	delete[] b_U_D;
	delete[] b_L_D;
	delete[] b_U_T_D;
	delete[] b_L_T_D;
	delete[] b_U_M;
	delete[] b_L_M;
	delete[] b_U_T_M;
	delete[] b_L_T_M;
	delete[] b_U_D_M;
	delete[] b_L_D_M;
	delete[] b_U_T_D_M;
	delete[] b_L_T_D_M;
	
}

void testTRSV_blk_gmp(int n, int m, int prec) {
	
	mpf_set_default_prec(prec);
	
	//Build A and x, then compute b:
	mpf_class *A_U = new mpf_class[n*n];
	mpf_class *A_L = new mpf_class[n*n];
	mpf_class *A_U_D = new mpf_class[n*n];	//=A_U with 1s on diagonal
	mpf_class *A_L_D = new mpf_class[n*n];	//=A_L with 1s on diagonal
	mpf_class *x = new mpf_class[n*m];
	mpf_class *b_U = new mpf_class[n*m];
	mpf_class *b_L = new mpf_class[n*m];
	mpf_class *b_U_T = new mpf_class[n*m];		//=b_L
	mpf_class *b_L_T = new mpf_class[n*m];		//=b_U
	mpf_class *b_U_D = new mpf_class[n*m];
	mpf_class *b_L_D = new mpf_class[n*m];
	mpf_class *b_U_T_D = new mpf_class[n*m];
	mpf_class *b_L_T_D = new mpf_class[n*m];
	
	mpf_class *b_U_X = new mpf_class[n*n];		//For XA=B
	mpf_class *b_L_X = new mpf_class[n*n];
	mpf_class *b_U_T_X = new mpf_class[n*n];		//=b_L
	mpf_class *b_L_T_X = new mpf_class[n*n];		//=b_U
	mpf_class *b_U_D_X = new mpf_class[n*n];
	mpf_class *b_L_D_X = new mpf_class[n*n];
	mpf_class *b_U_T_D_X = new mpf_class[n*n];
	mpf_class *b_L_T_D_X = new mpf_class[n*n];
	mpf_class *x_X = new mpf_class[n*n];
	
	//Fill A, x with random numbers:
	srand(time(NULL));
    for (int i=0; i<n; ++i) {
        for (int j=i; j<n; ++j) {
        
        	A_U[j*n +i] = 0;
        	A_L[i*n +j] = 0;
        	A_U_D[j*n +i] = 0;
        	A_L_D[i*n +j] = 0;
        	
        	mpf_class tmp = rand()/(mpf_class)RAND_MAX;
            A_U[i*n +j] = tmp;
            A_L[j*n +i] = tmp;
            
            A_U_D[i*n +j] = A_U[i*n +j];
            A_L_D[j*n +i] = A_L[j*n +i];
            
        }
        
        A_U_D[i*n +i] = 1;
        A_L_D[i*n +i] = 1;
        
        mpf_class tmp = rand()/(mpf_class)RAND_MAX;
        for (int j=0; j<m; ++j) {
        
        	x[i*m +j] = tmp;						//Columns of x all the same.
        	
        }
        
        for (int j=0; j<n; ++j) {
        
        	x_X[i*n +j] = tmp;						//Columns of x all the same.
        	
        }
        
    }
    
    //b = Ax
    //cout << "b:" << endl;
    for (int k=0; k<m; ++k) {
    	for (int i=0; i<n; ++i) {
    
    		b_U[i*m +k] = 0;
	    	b_L[i*m +k] = 0;
    	
    		b_U_D[i*m +k] = 0;
    		b_L_D[i*m +k] = 0;
    	
	    	for (int j=0; j<n; ++j) {
    	
    			b_U[i*m +k] += A_U[i*n +j] * x[j*m +k];
    			b_L[i*m +k] += A_L[i*n +j] * x[j*m +k];
    		
    			b_U_D[i*m +k] += A_U_D[i*n +j] * x[j*m +k];
	    		b_L_D[i*m +k] += A_L_D[i*n +j] * x[j*m +k];
    		
    		}

    		b_U_T[i*m +k] = b_L[i*m +k];
	    	b_L_T[i*m +k] = b_U[i*m +k];
    	
    		b_U_T_D[i*m +k] = b_L_D[i*m +k];
    		b_L_T_D[i*m +k] = b_U_D[i*m +k];
    	
    	//cout << b_U_T[i] << endl;
    	}
    }
    
    for (int k=0; k<n; ++k) {
    	for (int i=0; i<n; ++i) {

    		b_U_X[k*n +i] = 0;
	    	b_L_X[k*n +i] = 0;
    	
    		b_U_D_X[k*n +i] = 0;
    		b_L_D_X[k*n +i] = 0;
    		
    		for (int j=0; j<n; ++j) {
    		
    			b_U_X[k*n +i] += A_U[j*n +i] * x_X[k*n +j];
    			b_L_X[k*n +i] += A_L[j*n +i] * x_X[k*n +j];
    		
    			b_U_D_X[k*n +i] += A_U_D[j*n +i] * x_X[k*n +j];
	    		b_L_D_X[k*n +i] += A_L_D[j*n +i] * x_X[k*n +j];
	    		
	    	}
    	
    		b_U_T_X[k*n +i] = b_L_X[k*n +i];
	    	b_L_T_X[k*n +i] = b_U_X[k*n +i];
    	
    		b_U_T_D_X[k*n +i] = b_L_D_X[k*n +i];
    		b_L_T_D_X[k*n +i] = b_U_D_X[k*n +i];
    		
    	}
    }

	//Solve using the solver:
	int lda = n;
	int ldx = m;

	//AX=B	
	//Upper:
	trsv_blk('U', 'N', 'N', 'A', n, A_U, lda, b_U, ldx, m);
	trsv_blk('U', 'T', 'N', 'A', n, A_U, lda, b_U_T, ldx, m);
	trsv_blk('U', 'N', 'U', 'A', n, A_U_D, lda, b_U_D, ldx, m);
	trsv_blk('U', 'T', 'U', 'A', n, A_U_D, lda, b_U_T_D, ldx, m);

	//Lower:
	trsv_blk('L', 'N', 'N', 'A', n, A_L, lda, b_L, ldx, m);
	trsv_blk('L', 'T', 'N', 'A', n, A_L, lda, b_L_T, ldx, m);
	trsv_blk('L', 'N', 'U', 'A', n, A_L_D, lda, b_L_D, ldx, m);
	trsv_blk('L', 'T', 'U', 'A', n, A_L_D, lda, b_L_T_D, ldx, m);
	
	//XA=B
	ldx=n;
	//Upper:
	trsv_blk('U', 'N', 'N', 'X', n, A_U, lda, b_U_X, ldx, n);
	trsv_blk('U', 'T', 'N', 'X', n, A_U, lda, b_U_T_X, ldx, n);
	trsv_blk('U', 'N', 'U', 'X', n, A_U_D, lda, b_U_D_X, ldx, n);
	trsv_blk('U', 'T', 'U', 'X', n, A_U_D, lda, b_U_T_D_X, ldx, n);

	//Lower:
	trsv_blk('L', 'N', 'N', 'X', n, A_L, lda, b_L_X, ldx, n);
	trsv_blk('L', 'T', 'N', 'X', n, A_L, lda, b_L_T_X, ldx, n);
	trsv_blk('L', 'N', 'U', 'X', n, A_L_D, lda, b_L_D_X, ldx, n);
	trsv_blk('L', 'T', 'U', 'X', n, A_L_D, lda, b_L_T_D_X, ldx, n);
	
    
	//Output results:
	cout << "\nResults for trsv_blk:\n" << endl;
	
	mpf_class err_U = errCalc2d(n,m,x,b_U,prec);
	mpf_class err_U_T = errCalc2d(n,m,x,b_U_T,prec);
	mpf_class err_U_D = errCalc2d(n,m,x,b_U_D,prec);
	mpf_class err_U_T_D = errCalc2d(n,m,x,b_U_T_D,prec);
	
	mpf_class err_L = errCalc2d(n,m,x,b_L,prec);
	mpf_class err_L_T = errCalc2d(n,m,x,b_L_T,prec);
	mpf_class err_L_D = errCalc2d(n,m,x,b_L_D,prec);
	mpf_class err_L_T_D = errCalc2d(n,m,x,b_L_T_D,prec);
	
	mpf_class err_U_X = errCalc2d(n,n,x_X,b_U_X,prec);
	mpf_class err_U_T_X = errCalc2d(n,n,x_X,b_U_T_X,prec);
	mpf_class err_U_D_X = errCalc2d(n,n,x_X,b_U_D_X,prec);
	mpf_class err_U_T_D_X = errCalc2d(n,n,x_X,b_U_T_D_X,prec);
	
	mpf_class err_L_X = errCalc2d(n,n,x_X,b_L_X,prec);
	mpf_class err_L_T_X = errCalc2d(n,n,x_X,b_L_T_X,prec);
	mpf_class err_L_D_X = errCalc2d(n,n,x_X,b_L_D_X,prec);
	mpf_class err_L_T_D_X = errCalc2d(n,n,x_X,b_L_T_D_X,prec);
	
	cout << "AX=B:" << endl;
	cout << fixed << setprecision(4) << "Upper: " << err_U << "\t\t\t ----> " << ((err_U<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Upper transpose: " << err_U_T << "\t\t ----> " << ((err_U_T<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Upper unit: " << err_U_D << "\t\t ----> " << ((err_U_D<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Upper unit transpose: " << err_U_T_D << "\t ----> " << ((err_U_T_D<PREC)?"PASS":"FAIL") << endl;
		
	cout << fixed << setprecision(4) << "Lower: " << err_L << "\t\t\t ----> " << ((err_L<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Lower transpose: " << err_L_T << "\t\t ----> " << ((err_L_T<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Lower unit: " << err_L_D << "\t\t ----> " << ((err_L_D<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Lower unit transpose: " << err_L_T_D << "\t ----> " << ((err_L_T_D<PREC)?"PASS":"FAIL") << endl;

	cout << "XA=B" << endl;
	cout << fixed << setprecision(4) << "Upper: " << err_U_X << "\t\t\t ----> " << ((err_U_X<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Upper transpose: " << err_U_T_X << "\t\t ----> " << ((err_U_T_X<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Upper unit: " << err_U_D_X << "\t\t ----> " << ((err_U_D_X<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Upper unit transpose: " << err_U_T_D_X << "\t ----> " << ((err_U_T_D_X<PREC)?"PASS":"FAIL") << endl;
		
	cout << fixed << setprecision(4) << "Lower: " << err_L_X << "\t\t\t ----> " << ((err_L_X<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Lower transpose: " << err_L_T_X << "\t\t ----> " << ((err_L_T_X<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Lower unit: " << err_L_D_X << "\t\t ----> " << ((err_L_D_X<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Lower unit transpose: " << err_L_T_D_X << "\t ----> " << ((err_L_T_D_X<PREC)?"PASS":"FAIL") << endl;


	//Clean up:
	delete[] A_U;
	delete[] A_L;
	delete[] A_U_D;
	delete[] A_L_D;
	delete[] x;
	delete[] b_U;
	delete[] b_L;
	delete[] b_U_T;
	delete[] b_L_T;
	delete[] b_U_D;
	delete[] b_L_D;
	delete[] b_U_T_D;
	delete[] b_L_T_D;
	delete[] b_U_X;
	delete[] b_L_X;
	delete[] b_U_T_X;
	delete[] b_L_T_X;
	delete[] b_U_D_X;
	delete[] b_L_D_X;
	delete[] b_U_T_D_X;
	delete[] b_L_T_D_X;
	
}

void testTRSV_PSX_gmp(int n, int m, int prec) {

	mpf_set_default_prec(prec);
	
	//Build A and x, then compute b:
	mpf_class *A_U = new mpf_class[n*n];
	mpf_class *A_L = new mpf_class[n*n];
	mpf_class *A_U_D = new mpf_class[n*n];	//=A_U with 1s on diagonal
	mpf_class *A_L_D = new mpf_class[n*n];	//=A_L with 1s on diagonal
	mpf_class *x = new mpf_class[n*m];
	mpf_class *b_U = new mpf_class[n*m];
	mpf_class *b_L = new mpf_class[n*m];
	mpf_class *b_U_T = new mpf_class[n*m];		//=b_L
	mpf_class *b_L_T = new mpf_class[n*m];		//=b_U
	mpf_class *b_U_D = new mpf_class[n*m];
	mpf_class *b_L_D = new mpf_class[n*m];
	mpf_class *b_U_T_D = new mpf_class[n*m];
	mpf_class *b_L_T_D = new mpf_class[n*m];
	
	mpf_class *b_U_X = new mpf_class[n*n];		//For XA=B
	mpf_class *b_L_X = new mpf_class[n*n];
	mpf_class *b_U_T_X = new mpf_class[n*n];		//=b_L
	mpf_class *b_L_T_X = new mpf_class[n*n];		//=b_U
	mpf_class *b_U_D_X = new mpf_class[n*n];
	mpf_class *b_L_D_X = new mpf_class[n*n];
	mpf_class *b_U_T_D_X = new mpf_class[n*n];
	mpf_class *b_L_T_D_X = new mpf_class[n*n];
	mpf_class *x_X = new mpf_class[n*n];
	
	//Fill A, x with random numbers:
	srand(time(NULL));
    for (int i=0; i<n; ++i) {
        for (int j=i; j<n; ++j) {
        
        	A_U[j*n +i] = 0;
        	A_L[i*n +j] = 0;
        	A_U_D[j*n +i] = 0;
        	A_L_D[i*n +j] = 0;
        	
        	mpf_class tmp = rand()/(mpf_class)RAND_MAX;
            A_U[i*n +j] = tmp;
            A_L[j*n +i] = tmp;
            
            A_U_D[i*n +j] = A_U[i*n +j];
            A_L_D[j*n +i] = A_L[j*n +i];
            
        }
        
        A_U_D[i*n +i] = 1;
        A_L_D[i*n +i] = 1;
        
        mpf_class tmp = rand()/(mpf_class)RAND_MAX;
        for (int j=0; j<m; ++j) {
        
        	x[i*m +j] = tmp;						//Columns of x all the same.
        	
        }
        
        for (int j=0; j<n; ++j) {
        
        	x_X[i*n +j] = tmp;						//Columns of x all the same.
        	
        }        
        
    }

    //b = Ax
    //cout << "b:" << endl;
    for (int k=0; k<m; ++k) {
    	for (int i=0; i<n; ++i) {
    
    		b_U[i*m +k] = 0;
	    	b_L[i*m +k] = 0;
    	
    		b_U_D[i*m +k] = 0;
    		b_L_D[i*m +k] = 0;
    	
	    	for (int j=0; j<n; ++j) {
    	
    			b_U[i*m +k] += A_U[i*n +j] * x[j*m +k];
    			b_L[i*m +k] += A_L[i*n +j] * x[j*m +k];
    		
    			b_U_D[i*m +k] += A_U_D[i*n +j] * x[j*m +k];
	    		b_L_D[i*m +k] += A_L_D[i*n +j] * x[j*m +k];
    		
    		}

    		b_U_T[i*m +k] = b_L[i*m +k];
	    	b_L_T[i*m +k] = b_U[i*m +k];
    	
    		b_U_T_D[i*m +k] = b_L_D[i*m +k];
    		b_L_T_D[i*m +k] = b_U_D[i*m +k];
    	
    	}
    }
    
    for (int k=0; k<n; ++k) {
    	for (int i=0; i<n; ++i) {
    	
    		b_U_X[k*n +i] = 0;
	    	b_L_X[k*n +i] = 0;
    	
    		b_U_D_X[k*n +i] = 0;
    		b_L_D_X[k*n +i] = 0;
    		
    		for (int j=0; j<n; ++j) {
    		
    			b_U_X[k*n +i] += A_U[j*n +i] * x_X[k*n +j];
    			b_L_X[k*n +i] += A_L[j*n +i] * x_X[k*n +j];
    		
    			b_U_D_X[k*n +i] += A_U_D[j*n +i] * x_X[k*n +j];
	    		b_L_D_X[k*n +i] += A_L_D[j*n +i] * x_X[k*n +j];
	    		
	    	}
    	
    		b_U_T_X[k*n +i] = b_L_X[k*n +i];
	    	b_L_T_X[k*n +i] = b_U_X[k*n +i];
    	
    		b_U_T_D_X[k*n +i] = b_L_D_X[k*n +i];
    		b_L_T_D_X[k*n +i] = b_U_D_X[k*n +i];
    		
    	}
    }

	//Solve using the solver:
	int lda = n;
	int ldx = m;

	//AX=B
	//Upper:
	trsv_psx('U', 'N', 'N', 'A', n, A_U, lda, b_U, ldx, m);
	trsv_psx('U', 'T', 'N', 'A', n, A_U, lda, b_U_T, ldx, m);
	trsv_psx('U', 'N', 'U', 'A', n, A_U_D, lda, b_U_D, ldx, m);
	trsv_psx('U', 'T', 'U', 'A', n, A_U_D, lda, b_U_T_D, ldx, m);

	//Lower:
	trsv_psx('L', 'N', 'N', 'A', n, A_L, lda, b_L, ldx, m);
	trsv_psx('L', 'T', 'N', 'A', n, A_L, lda, b_L_T, ldx, m);
	trsv_psx('L', 'N', 'U', 'A', n, A_L_D, lda, b_L_D, ldx, m);
	trsv_psx('L', 'T', 'U', 'A', n, A_L_D, lda, b_L_T_D, ldx, m);
	
	//XA=B
	ldx=n;
	//Upper:
	trsv_psx('U', 'N', 'N', 'X', n, A_U, lda, b_U_X, ldx, n);
	trsv_psx('U', 'T', 'N', 'X', n, A_U, lda, b_U_T_X, ldx, n);
	trsv_psx('U', 'N', 'U', 'X', n, A_U_D, lda, b_U_D_X, ldx, n);
	trsv_psx('U', 'T', 'U', 'X', n, A_U_D, lda, b_U_T_D_X, ldx, n);

	//Lower:
	trsv_psx('L', 'N', 'N', 'X', n, A_L, lda, b_L_X, ldx, n);
	trsv_psx('L', 'T', 'N', 'X', n, A_L, lda, b_L_T_X, ldx, n);
	trsv_psx('L', 'N', 'U', 'X', n, A_L_D, lda, b_L_D_X, ldx, n);
	trsv_psx('L', 'T', 'U', 'X', n, A_L_D, lda, b_L_T_D_X, ldx, n);

    //Output results:
	cout << "\nResults for trsv_psx:\n" << endl;
	
	mpf_class err_U = errCalc2d(n,m,x,b_U,prec);
	mpf_class err_U_T = errCalc2d(n,m,x,b_U_T,prec);
	mpf_class err_U_D = errCalc2d(n,m,x,b_U_D,prec);
	mpf_class err_U_T_D = errCalc2d(n,m,x,b_U_T_D,prec);
	
	mpf_class err_L = errCalc2d(n,m,x,b_L,prec);
	mpf_class err_L_T = errCalc2d(n,m,x,b_L_T,prec);
	mpf_class err_L_D = errCalc2d(n,m,x,b_L_D,prec);
	mpf_class err_L_T_D = errCalc2d(n,m,x,b_L_T_D,prec);
	
	mpf_class err_U_X = errCalc2d(n,n,x_X,b_U_X,prec);
	mpf_class err_U_T_X = errCalc2d(n,n,x_X,b_U_T_X,prec);
	mpf_class err_U_D_X = errCalc2d(n,n,x_X,b_U_D_X,prec);
	mpf_class err_U_T_D_X = errCalc2d(n,n,x_X,b_U_T_D_X,prec);
	
	mpf_class err_L_X = errCalc2d(n,n,x_X,b_L_X,prec);
	mpf_class err_L_T_X = errCalc2d(n,n,x_X,b_L_T_X,prec);
	mpf_class err_L_D_X = errCalc2d(n,n,x_X,b_L_D_X,prec);
	mpf_class err_L_T_D_X = errCalc2d(n,n,x_X,b_L_T_D_X,prec);
	
	cout << "AX=B:" << endl;
	cout << fixed << setprecision(4) << "Upper: " << err_U << "\t\t\t ----> " << ((err_U<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Upper transpose: " << err_U_T << "\t\t ----> " << ((err_U_T<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Upper unit: " << err_U_D << "\t\t ----> " << ((err_U_D<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Upper unit transpose: " << err_U_T_D << "\t ----> " << ((err_U_T_D<PREC)?"PASS":"FAIL") << endl;
		
	cout << fixed << setprecision(4) << "Lower: " << err_L << "\t\t\t ----> " << ((err_L<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Lower transpose: " << err_L_T << "\t\t ----> " << ((err_L_T<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Lower unit: " << err_L_D << "\t\t ----> " << ((err_L_D<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Lower unit transpose: " << err_L_T_D << "\t ----> " << ((err_L_T_D<PREC)?"PASS":"FAIL") << endl;

	cout << "XA=B" << endl;
	cout << fixed << setprecision(4) << "Upper: " << err_U_X << "\t\t\t ----> " << ((err_U_X<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Upper transpose: " << err_U_T_X << "\t\t ----> " << ((err_U_T_X<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Upper unit: " << err_U_D_X << "\t\t ----> " << ((err_U_D_X<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Upper unit transpose: " << err_U_T_D_X << "\t ----> " << ((err_U_T_D_X<PREC)?"PASS":"FAIL") << endl;
		
	cout << fixed << setprecision(4) << "Lower: " << err_L_X << "\t\t\t ----> " << ((err_L_X<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Lower transpose: " << err_L_T_X << "\t\t ----> " << ((err_L_T_X<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Lower unit: " << err_L_D_X << "\t\t ----> " << ((err_L_D_X<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Lower unit transpose: " << err_L_T_D_X << "\t ----> " << ((err_L_T_D_X<PREC)?"PASS":"FAIL") << endl;

	//Timing tests:
	Timer timer;
	ldx=m;
	timer.start();
	trsv_blk('U', 'N', 'N', 'A', n, A_U, lda, b_U, ldx, m);
	timer.stop();
	cout << "\nSingle threaded calculation took: " << timer.elapsed() << endl;
	
	timer.start();
	trsv_psx('U', 'N', 'N', 'A', n, A_U, lda, b_U, ldx, m);
	timer.stop();
	cout << "Using " << NUMPROC << " threads took: " << timer.elapsed() << endl;

	//Clean up:
	delete[] A_U;
	delete[] A_L;
	delete[] A_U_D;
	delete[] A_L_D;
	delete[] x;
	delete[] b_U;
	delete[] b_L;
	delete[] b_U_T;
	delete[] b_L_T;
	delete[] b_U_D;
	delete[] b_L_D;
	delete[] b_U_T_D;
	delete[] b_L_T_D;
	delete[] b_U_X;
	delete[] b_L_X;
	delete[] b_U_T_X;
	delete[] b_L_T_X;
	delete[] b_U_D_X;
	delete[] b_L_D_X;
	delete[] b_U_T_D_X;
	delete[] b_L_T_D_X;
	
}

mpf_class errCalc(int n, mpf_class *exact, mpf_class *calc, int incX, int prec) {

	mpf_set_default_prec(prec);
	
	mpf_class error = 0;
	
	int firstInxX = (copysign(1,incX) < 0) ? (n-1)*(-incX) : 0;
	
	for (int i=0, iX=firstInxX; i<n; ++i, iX+=incX) {
	
		error += abs(exact[i] - calc[iX]);
		
	}
	
	return error;
		
}

mpf_class errCalc2d(int n, int m, mpf_class *exact, mpf_class *calc, int prec) {

	mpf_set_default_prec(prec);
	
	mpf_class error = 0;
	
	for (int i=0; i<n; ++i) {
		for (int j=0; j<m; ++j) {
		
			error += abs(exact[i*m +j] - calc[i*m +j]);
			
		}
	}
	
	return error;
	
}