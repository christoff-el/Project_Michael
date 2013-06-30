#define PREC 1e-5
#define NUMPROC 2
#define CACHESIZE 32000

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <gmpxx.h>

#include "../trsv.h"
#include "../trsv_mat.h"
#include "../trsv_blk.h"
#include "../trsv_psx.h"
#include "timer.h"

#include "test_trsv_tests.h"

using namespace std;

/*
Input:
	n - size of the system (i.e. A is n x n)
	m - number of columns of X
	
	test_mode:
		1 - test trsv
		2 - test trsv_mat
		3 - test trsv_blk
		4 - test trsv_psx
		5 - test all
		6 - test performance
		7 - EVERYTHING!
	
	precision_mode:
		0 - standard double precision
		x>0 - gmp precision of x

	g++ test_trsv.cpp -o b.trsv -Wall -O3 -std=gnu++11 -lgmp -lgmpxx
*/


int main(int argc, char **argv) {

	if (argc != 5) {
		cerr << "Usage: " << argv[0] << " n m test_mode precision_mode" << endl;
		return 1;
	}
	
	//System size:
	int n = atoi(argv[1]);
	int m = atoi(argv[2]);
	
	//Test mode:
	int tm = atoi(argv[3]);
	
	//Precision mode:
	int prec = atoi(argv[4]);
	
	if (prec <= 0) {
	
		if (n>40) {
			cout << "Warning: routine failures may be due to truncation errors!\n" << endl;
		}
	
		double dummy;
		
		if (tm == 1 || tm == 5 || tm == 7) {testTRSV(dummy,n,1);}
		if (tm == 2 || tm == 5 || tm == 7) {testTRSV_MAT(dummy,n,m,1);}
		if (tm == 3 || tm == 5 || tm == 7) {testTRSV_BLK(dummy,n,m,1);}
		if (tm == 4 || tm == 5 || tm == 7) {testTRSV_PSX(dummy,n,m,1);}
		if (tm == 6 || tm == 7) {performance(dummy,n,m,1);}	
	
	}
	else {
	
		mpf_set_default_prec(prec);
		mpf_class dummy;
		
		if (tm == 1 || tm == 5 || tm == 7) {testTRSV(dummy,n,prec);}
		if (tm == 2 || tm == 5 || tm == 7) {testTRSV_MAT(dummy,n,m,prec);}
		if (tm == 3 || tm == 5 || tm == 7) {testTRSV_BLK(dummy,n,m,prec);}
		if (tm == 4 || tm == 5 || tm == 7) {testTRSV_PSX(dummy,n,m,prec);}
		if (tm == 6 || tm == 7) {performance(dummy,n,m,prec);}	
	
	}
	
	return 0;
		
}


