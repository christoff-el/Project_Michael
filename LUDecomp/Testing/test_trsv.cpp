#define PREC 1e-5
#define NUMPROC 2

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include "../trsv.h"
#include "../trsv_mat.h"
#include "../trsv_psx.h"
#include "timer.h"

using namespace std;

double errCalc(int n, double *exact, double *calc, int incX);
double errCalc2d(int n, int m, double *exact, double *calc);
void testTRSV(int n);
void testTRSV_mat(int n, int m);
void testTRSV_PSX(int n, int m);


int main(int argc, char **argv) {

	if (argc != 3) {
		cerr << "Usage: " << argv[0] << " n m" << endl;
		return 1;
	}
	
	//System size:
	int n=atoi(argv[1]);
	int m=atoi(argv[2]);
	
	if (n>40) {
		cout << "Warning: routine failures may be due to truncation errors!\n" << endl;
	}
	
	testTRSV(n);
	testTRSV_mat(n,m);
	testTRSV_PSX(n,m);
	
	
	return 0;
	
}


void testTRSV(int n) {

	//Build A and x, then compute b:
	double *A_U = new double[n*n];
	double *A_L = new double[n*n];
	double *A_U_D = new double[n*n];	//=A_U with 1s on diagonal
	double *A_L_D = new double[n*n];	//=A_L with 1s on diagonal
	double *x = new double[n];
	double *b_U = new double[n];
	double *b_L = new double[n];
	double *b_U_T = new double[n];		//=b_L
	double *b_L_T = new double[n];		//=b_U
	double *b_U_D = new double[n];
	double *b_L_D = new double[n];
	double *b_U_T_D = new double[n];
	double *b_L_T_D = new double[n];
	double *b_U_M = new double[n];			//_M same as above, but backwards.
	double *b_L_M = new double[n];
	double *b_U_T_M = new double[n];		//=b_L_M
	double *b_L_T_M = new double[n];		//=b_U_M
	double *b_U_D_M = new double[n];
	double *b_L_D_M = new double[n];
	double *b_U_T_D_M = new double[n];
	double *b_L_T_D_M = new double[n];
	
	//Fill A, x with random numbers:
	srand(time(NULL));
    for (int i=0; i<n; ++i) {
        for (int j=i; j<n; ++j) {
        
        	A_U[j*n +i] = 0;
        	A_L[i*n +j] = 0;
        	
        	double tmp=rand()/(double)RAND_MAX;
            A_U[i*n +j] = tmp;
            A_L[j*n +i] = tmp;
            
            A_U_D[i*n +j] = A_U[i*n +j];
            A_L_D[j*n +i] = A_L[j*n +i];
            
        }
        
        A_U_D[i*n +i] = 1;
        A_L_D[i*n +i] = 1;
        
        x[i] = rand()/(double)RAND_MAX;
        
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
	
	double err_U = errCalc(n,x,b_U,1);
	double err_U_T = errCalc(n,x,b_U_T,1);
	double err_U_D = errCalc(n,x,b_U_D,1);
	double err_U_T_D = errCalc(n,x,b_U_T_D,1);
	
	double err_L = errCalc(n,x,b_L,1);
	double err_L_T = errCalc(n,x,b_L_T,1);
	double err_L_D = errCalc(n,x,b_L_D,1);
	double err_L_T_D = errCalc(n,x,b_L_T_D,1);
	
	double err_U_M = errCalc(n,x,b_U_M,-1);
	double err_U_T_M = errCalc(n,x,b_U_T_M,-1);
	double err_U_D_M = errCalc(n,x,b_U_D_M,-1);
	double err_U_T_D_M = errCalc(n,x,b_U_T_D_M,-1);
	
	double err_L_M = errCalc(n,x,b_L_M,-1);
	double err_L_T_M = errCalc(n,x,b_L_T_M,-1);
	double err_L_D_M = errCalc(n,x,b_L_D_M,-1);
	double err_L_T_D_M = errCalc(n,x,b_L_T_D_M,-1);
	
	cout << "incX = 1:" << endl;
	cout << fixed << setprecision(4) << "Upper: " << err_U << "\t\t\t ----> " << ((err_U<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Upper transpose: " << err_U_T << "\t\t ----> " << ((err_U_T<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Upper unit: " << err_U_D << "\t\t ----> " << ((err_U_D<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Upper unit transpose: " << err_U_T_D << "\t ----> " << ((err_U_T_D<PREC)?"PASS":"FAIL") << endl;
		
	cout << fixed << setprecision(4) << "Lower: " << err_L << "\t\t\t ----> " << ((err_L<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Lower transpose: " << err_L_T << "\t\t ----> " << ((err_L_T<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Lower unit: " << err_L_D << "\t\t ----> " << ((err_L_D<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Lower unit transpose: " << err_L_T_D << "\t ----> " << ((err_L_T_D<PREC)?"PASS":"FAIL") << endl;

	cout << "\nincX = -1:" << endl;
	cout << fixed << setprecision(4) << "Upper: " << err_U_M << "\t\t\t ----> " << ((err_U_M<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Upper transpose: " << err_U_T_M << "\t\t ----> " << ((err_U_T_M<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Upper unit: " << err_U_D_M << "\t\t ----> " << ((err_U_D_M<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Upper unit transpose: " << err_U_T_D_M << "\t ----> " << ((err_U_T_D_M<PREC)?"PASS":"FAIL") << endl;
		
	cout << fixed << setprecision(4) << "Lower: " << err_L_M << "\t\t\t ----> " << ((err_L_M<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Lower transpose: " << err_L_T_M << "\t\t ----> " << ((err_L_T_M<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Lower unit: " << err_L_D_M << "\t\t ----> " << ((err_L_D_M<PREC)?"PASS":"FAIL") << endl;
	cout << fixed << setprecision(4) << "Lower unit transpose: " << err_L_T_D_M << "\t ----> " << ((err_L_T_D_M<PREC)?"PASS":"FAIL") << endl;


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

void testTRSV_mat(int n, int m) {

	//Build A and x, then compute b:
	double *A_U = new double[n*n];
	double *A_L = new double[n*n];
	double *A_U_D = new double[n*n];	//=A_U with 1s on diagonal
	double *A_L_D = new double[n*n];	//=A_L with 1s on diagonal
	double *x = new double[n*m];
	double *b_U = new double[n*m];
	double *b_L = new double[n*m];
	double *b_U_T = new double[n*m];		//=b_L
	double *b_L_T = new double[n*m];		//=b_U
	double *b_U_D = new double[n*m];
	double *b_L_D = new double[n*m];
	double *b_U_T_D = new double[n*m];
	double *b_L_T_D = new double[n*m];

	double *b_U_X = new double[n*n];		//For XA=B
	double *b_L_X = new double[n*n];
	double *b_U_T_X = new double[n*n];		//=b_L
	double *b_L_T_X = new double[n*n];		//=b_U
	double *b_U_D_X = new double[n*n];
	double *b_L_D_X = new double[n*n];
	double *b_U_T_D_X = new double[n*n];
	double *b_L_T_D_X = new double[n*n];
	
	//Fill A, x with random numbers:
	srand(time(NULL));
    for (int i=0; i<n; ++i) {
        for (int j=i; j<n; ++j) {
        
        	A_U[j*n +i] = 0;
        	A_L[i*n +j] = 0;
        	A_U_D[j*n +i] = 0;
        	A_L_D[i*n +j] = 0;
        	
        	double tmp = rand()/(double)RAND_MAX;
            A_U[i*n +j] = tmp;
            A_L[j*n +i] = tmp;
            
            A_U_D[i*n +j] = A_U[i*n +j];
            A_L_D[j*n +i] = A_L[j*n +i];
            
        }
        
        A_U_D[i*n +i] = 1;
        A_L_D[i*n +i] = 1;
        
        double tmp = rand()/(double)RAND_MAX;
        for (int j=0; j<m; ++j) {
        
        	x[i*m +j] = tmp;						//Columns of x all the same.
        	
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
    		
    			b_U_X[k*n +i] += A_U[j*n +i] * x[k*n +j];
    			b_L_X[k*n +i] += A_L[j*n +i] * x[k*n +j];
    		
    			b_U_D_X[k*n +i] += A_U_D[j*n +i] * x[k*n +j];
	    		b_L_D_X[k*n +i] += A_L_D[j*n +i] * x[k*n +j];
	    		
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
	
	for (int i=0; i<n; ++i) {
		for (int j=0; j<n; ++j) {
		
			cout << A_U[i*lda+j] << " ";
			
		}
		cout<<endl;
	}

	for (int i=0; i<n; ++i) {
		for (int j=0; j<n; ++j) {
		
			cout << b_U_X[i*ldx+j] << " ";
			
		}
		cout<<endl;
	}
	
	//AX=B	
	//Upper:
	trsv_mat('U', 'N', 'N', 'A', n, A_U, lda, b_U, ldx, m);
	trsv_mat('U', 'T', 'N', 'A', n, A_U, lda, b_U_T, ldx, m);
	trsv_mat('U', 'N', 'U', 'A', n, A_U_D, lda, b_U_D, ldx, m);
	trsv_mat('U', 'T', 'U', 'A', n, A_U_D, lda, b_U_T_D, ldx, m);

	//Lower:
	trsv_mat('L', 'N', 'N', 'A', n, A_L, lda, b_L, ldx, m);
	trsv_mat('L', 'T', 'N', 'A', n, A_L, lda, b_L_T, ldx, m);
	trsv_mat('L', 'N', 'U', 'A', n, A_L_D, lda, b_L_D, ldx, m);
	trsv_mat('L', 'T', 'U', 'A', n, A_L_D, lda, b_L_T_D, ldx, m);
	
	//XA=B
	//Upper:
	trsv_mat('U', 'N', 'N', 'X', n, A_U, lda, b_U_X, ldx, m);
	trsv_mat('U', 'T', 'N', 'X', n, A_U, lda, b_U_T_X, ldx, m);
	trsv_mat('U', 'N', 'U', 'X', n, A_U_D, lda, b_U_D_X, ldx, m);
	trsv_mat('U', 'T', 'U', 'X', n, A_U_D, lda, b_U_T_D_X, ldx, m);

	//Lower:
	trsv_mat('L', 'N', 'N', 'X', n, A_L, lda, b_L_X, ldx, m);
	trsv_mat('L', 'T', 'N', 'X', n, A_L, lda, b_L_T_X, ldx, m);
	trsv_mat('L', 'N', 'U', 'X', n, A_L_D, lda, b_L_D_X, ldx, m);
	trsv_mat('L', 'T', 'U', 'X', n, A_L_D, lda, b_L_T_D_X, ldx, m);

	for (int i=0; i<n; ++i) {
		for (int j=0; j<n; ++j) {
		
			cout << x[i*ldx+j] << " ";
			
		}
		cout<<endl;
	}
    
	//Output results:
	cout << "\nResults for trsv_mat:\n" << endl;
	
	double err_U = errCalc2d(n,m,x,b_U);
	double err_U_T = errCalc2d(n,m,x,b_U_T);
	double err_U_D = errCalc2d(n,m,x,b_U_D);
	double err_U_T_D = errCalc2d(n,m,x,b_U_T_D);
	
	double err_L = errCalc2d(n,m,x,b_L);
	double err_L_T = errCalc2d(n,m,x,b_L_T);
	double err_L_D = errCalc2d(n,m,x,b_L_D);
	double err_L_T_D = errCalc2d(n,m,x,b_L_T_D);
	
	double err_U_X = errCalc2d(n,m,x,b_U_X);
	double err_U_T_X = errCalc2d(n,m,x,b_U_T_X);
	double err_U_D_X = errCalc2d(n,m,x,b_U_D_X);
	double err_U_T_D_X = errCalc2d(n,m,x,b_U_T_D_X);
	
	double err_L_X = errCalc2d(n,m,x,b_L_X);
	double err_L_T_X = errCalc2d(n,m,x,b_L_T_X);
	double err_L_D_X = errCalc2d(n,m,x,b_L_D_X);
	double err_L_T_D_X = errCalc2d(n,m,x,b_L_T_D_X);
	
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

void testTRSV_PSX(int n, int m) {

	//Build A and x, then compute b:
	double *A_U = new double[n*n];
	double *A_L = new double[n*n];
	double *A_U_D = new double[n*n];	//=A_U with 1s on diagonal
	double *A_L_D = new double[n*n];	//=A_L with 1s on diagonal
	double *x = new double[n*m];
	double *b_U = new double[n*m];
	double *b_L = new double[n*m];
	double *b_U_T = new double[n*m];		//=b_L
	double *b_L_T = new double[n*m];		//=b_U
	double *b_U_D = new double[n*m];
	double *b_L_D = new double[n*m];
	double *b_U_T_D = new double[n*m];
	double *b_L_T_D = new double[n*m];
	
	double *b_U_X = new double[n*n];		//For XA=B
	double *b_L_X = new double[n*n];
	double *b_U_T_X = new double[n*n];		//=b_L
	double *b_L_T_X = new double[n*n];		//=b_U
	double *b_U_D_X = new double[n*n];
	double *b_L_D_X = new double[n*n];
	double *b_U_T_D_X = new double[n*n];
	double *b_L_T_D_X = new double[n*n];
	
	//Fill A, x with random numbers:
	srand(time(NULL));
    for (int i=0; i<n; ++i) {
        for (int j=i; j<n; ++j) {
        
        	A_U[j*n +i] = 0;
        	A_L[i*n +j] = 0;
        	A_U_D[j*n +i] = 0;
        	A_L_D[i*n +j] = 0;
        	
        	double tmp = rand()/(double)RAND_MAX;
            A_U[i*n +j] = tmp;
            A_L[j*n +i] = tmp;
            
            A_U_D[i*n +j] = A_U[i*n +j];
            A_L_D[j*n +i] = A_L[j*n +i];
            
        }
        
        A_U_D[i*n +i] = 1;
        A_L_D[i*n +i] = 1;
        
        double tmp = rand()/(double)RAND_MAX;
        for (int j=0; j<m; ++j) {
        
        	x[i*m +j] = tmp;						//Columns of x all the same.
        	
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
    		
    			b_U_X[k*n +i] += A_U[j*n +i] * x[k*n +j];
    			b_L_X[k*n +i] += A_L[j*n +i] * x[k*n +j];
    		
    			b_U_D_X[k*n +i] += A_U_D[j*n +i] * x[k*n +j];
	    		b_L_D_X[k*n +i] += A_L_D[j*n +i] * x[k*n +j];
	    		
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
	//Upper:
	trsv_mat('U', 'N', 'N', 'X', n, A_U, lda, b_U_X, ldx, m);
	trsv_mat('U', 'T', 'N', 'X', n, A_U, lda, b_U_T_X, ldx, m);
	trsv_mat('U', 'N', 'U', 'X', n, A_U_D, lda, b_U_D_X, ldx, m);
	trsv_mat('U', 'T', 'U', 'X', n, A_U_D, lda, b_U_T_D_X, ldx, m);

	//Lower:
	trsv_mat('L', 'N', 'N', 'X', n, A_L, lda, b_L_X, ldx, m);
	trsv_mat('L', 'T', 'N', 'X', n, A_L, lda, b_L_T_X, ldx, m);
	trsv_mat('L', 'N', 'U', 'X', n, A_L_D, lda, b_L_D_X, ldx, m);
	trsv_mat('L', 'T', 'U', 'X', n, A_L_D, lda, b_L_T_D_X, ldx, m);

    
	//Output results:
	cout << "\nResults for trsv_psx:\n" << endl;
	
	double err_U = errCalc2d(n,m,x,b_U);
	double err_U_T = errCalc2d(n,m,x,b_U_T);
	double err_U_D = errCalc2d(n,m,x,b_U_D);
	double err_U_T_D = errCalc2d(n,m,x,b_U_T_D);
	
	double err_L = errCalc2d(n,m,x,b_L);
	double err_L_T = errCalc2d(n,m,x,b_L_T);
	double err_L_D = errCalc2d(n,m,x,b_L_D);
	double err_L_T_D = errCalc2d(n,m,x,b_L_T_D);
	
	double err_U_X = errCalc2d(n,m,x,b_U_X);
	double err_U_T_X = errCalc2d(n,m,x,b_U_T_X);
	double err_U_D_X = errCalc2d(n,m,x,b_U_D_X);
	double err_U_T_D_X = errCalc2d(n,m,x,b_U_T_D_X);
	
	double err_L_X = errCalc2d(n,m,x,b_L_X);
	double err_L_T_X = errCalc2d(n,m,x,b_L_T_X);
	double err_L_D_X = errCalc2d(n,m,x,b_L_D_X);
	double err_L_T_D_X = errCalc2d(n,m,x,b_L_T_D_X);
	
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
	
	timer.start();
	trsv_mat('U', 'N', 'N', 'A', n, A_U, lda, b_U, ldx, m);
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

double errCalc(int n, double *exact, double *calc, int incX) {

	double error = 0;
	
	int firstInxX = (copysign(1,incX) < 0) ? (n-1)*(-incX) : 0;
	
	for (int i=0, iX=firstInxX; i<n; ++i, iX+=incX) {
	
		error += fabs(exact[i] - calc[iX]);
		
	}
	
	return error;
		
}

double errCalc2d(int n, int m, double *exact, double *calc) {

	double error = 0;
	
	for (int i=0; i<n; ++i) {
		for (int j=0; j<m; ++j) {
		
			error += fabs(exact[i*m +j] - calc[i*m +j]);
			
		}
	}
	
	return error;
	
}