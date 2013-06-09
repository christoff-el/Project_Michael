#ifndef LU_DECOMP_TRSV_TCC
#define LU_DECOMP_TRSV_TCC 1

#include <cmath>

template <typename IndexType, typename A, typename X>
void
trsv(char uplo, char trans, char diag, IndexType n,
		const A *a, IndexType lda,
				X *x, IndexType incX)
{

IndexType firstInxX = (copysign(1,incX) < 0) ? (n-1)*(-incX) : 0;

/*****UPPER TRIANGULAR*****/
if (uplo == 'U' || uplo == 'u') {

	/*****UPPER TRIANGULAR / UNIT TRIANGULAR*****/
	if (diag == 'U' || diag == 'u') {

		/*****UPPER TRIANGULAR / UNIT TRIANGULAR / NO TRANSPOSE*****/
		if (trans == 'N' || trans == 'n') {

			for (int i=n-1, iX=firstInxX+((n-1)*incX); i>=0; --i, iX-=incX) {
				for (int j=i+1, jX=iX+incX; j<n; ++j, jX+=incX) {

					x[iX] -= a[i*lda +j] * x[jX];

				}
			}

		}

		/*****UPPER TRIANGULAR / UNIT TRIANGULAR / TRANSPOSE*****/
		else if (trans == 'T' || trans == 't' || trans == 'C' || trans == 'c') {

			for (int i=0, iX=firstInxX; i<n; ++i, iX+=incX) {
				for (int j=0, jX=firstInxX; j<i; ++j, jX+=incX) {

					x[iX] -= a[j*lda +i] * x[jX];

				}
			}

		}

	}	//UNIT TRIANGULAR

	/*****UPPER TRIANGULAR / NOT UNIT TRIANGULAR*****/
	else if (diag == 'N' || diag == 'n') {

		/*****UPPER TRIANGULAR / NOT UNIT TRIANGULAR / NO TRANSPOSE*****/
		if (trans == 'N' || trans == 'n') {

			for (int i=n-1, iX=firstInxX+((n-1)*incX); i>=0; --i, iX-=incX) {

				for (int j=i+1, jX=iX+incX; j<n; ++j, jX+=incX) {

					x[iX] -= a[i*lda +j] * x[jX];

				}

				x[iX] /= a[i*lda +i];

			}

		}

		/*****UPPER TRIANGULAR / NOT UNIT TRIANGULAR / TRANSPOSE*****/
		else if (trans == 'T' || trans == 't' || trans == 'C' || trans == 'c') {

			for (int i=0, iX=firstInxX; i<n; ++i, iX+=incX) {

				for (int j=0, jX=firstInxX; j<i; ++j, jX+=incX) {

					x[iX] -= a[j*lda +i] * x[jX];

				}

				x[iX] /= a[i*lda +i];

			}

		}

	}	//NOT UNIT TRIANGULAR

}	//UPPER TRIANGULAR

/*****LOWER TRIANGULAR*****/
else if (uplo == 'L' || uplo == 'l') {

	/*****LOWER TRIANGULAR / UNIT TRIANGULAR*****/
	if (diag == 'U' || diag == 'u') {
		/*****LOWER TRIANGULAR / UNIT TRIANGULAR / NO TRANSPOSE*****/
		if (trans == 'N' || trans == 'n') {

			for (int i=0, iX=firstInxX; i<n; ++i, iX+=incX) {
				for (int j=0, jX=firstInxX; j<i; ++j, jX+=incX) {
			
					x[iX] -= a[i*lda +j] * x[jX];

				}
			}

		}
	
		/*****LOWER TRIANGULAR / UNIT TRIANGULAR / TRANSPOSE*****/
		else if (trans == 'T' || trans == 't' || trans == 'C' || trans == 'c') {

			for (int i=n-1, iX=firstInxX+((n-1)*incX); i>=0; --i, iX-=incX) {
				for (int j=i+1, jX=iX+incX; j<n; ++j, jX+=incX) {

					x[iX] -= a[j*lda +i] * x[jX];

				}
			}

		}
		
	}	//UNIT TRIANGULAR
	
	/*****LOWER TRIANGULAR / NOT UNIT TRIANGULAR*****/
	else if (diag == 'N' || diag == 'n') {
	
		/*****LOWER TRIANGULAR / NOT UNIT TRIANGULAR / NO TRANSPOSE*****/
		if (trans == 'N' || trans == 'n') {

			for (int i=0, iX=firstInxX; i<n; ++i, iX+=incX) {

				for (int j=0, jX=firstInxX; j<i; ++j, jX+=incX) {
			
					x[iX] -= a[i*lda +j] * x[jX];
				
				}

				x[iX] /= a[i*lda +i];

			}

		}
	
		/*****LOWER TRIANGULAR / NOT UNIT TRIANGULAR / TRANSPOSE*****/
		else if (trans == 'T' || trans == 't' || trans == 'C' || trans == 'c') {
	
			for (int i=n-1, iX=firstInxX+((n-1)*incX); i>=0; --i, iX-=incX) {

				for (int j=i+1, jX=iX+incX; j<n; ++j, jX+=incX) {

					x[iX] -= a[j*lda +i] * x[jX];

				}

				x[iX] /= a[i*lda +i];

			}

		}

	}	//NOT UNIT TRIANGULAR

}	//LOWER TRIANGULAR

}

#endif	// LU_DECOMP_TRSV_TCC
