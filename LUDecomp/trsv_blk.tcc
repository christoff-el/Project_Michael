#ifndef LU_DECOMP_TRSV_BLK_TCC
#define LU_DECOMP_TRSV_BLK_TCC 1

template <typename IndexType, typename A, typename X>
void
trsv_blk(char uplo, char trans, char diag, char ordering, IndexType n,
		const A *a, IndexType lda,
				X *x, IndexType ldx, 
					IndexType m)
{

/***** AX = B ORDERING *****/
if (ordering == 'A' || ordering == 'a') {

	/*****UPPER TRIANGULAR*****/
	if (uplo == 'U' || uplo == 'u') {

		/*****UPPER TRIANGULAR / UNIT TRIANGULAR*****/
		if (diag == 'U' || diag == 'u') {

			/*****UPPER TRIANGULAR / UNIT TRIANGULAR / NO TRANSPOSE*****/
			if (trans == 'N' || trans == 'n') {

				for (int i=n-1; i>=0; --i) {
					for (int j=i+1; j<n; ++j) {
						for (int k=0; k<m; ++k) {			//For each system

							x[i*ldx +k] -= a[i*lda +j] * x[j*ldx +k];
						
						}
					}
				}

			}

			/*****UPPER TRIANGULAR / UNIT TRIANGULAR / TRANSPOSE*****/
			else if (trans == 'T' || trans == 't') {

				for (int i=0; i<n; ++i) {
					for (int j=0; j<i; ++j) {
						for (int k=0; k<m; ++k) {			//For each system

							x[i*ldx +k] -= a[j*lda +i] * x[j*ldx +k];
						
						}

					}
				}

			}

		}	//UNIT TRIANGULAR

		/*****UPPER TRIANGULAR / NOT UNIT TRIANGULAR*****/
		else if (diag == 'N' || diag == 'n') {

			/*****UPPER TRIANGULAR / NOT UNIT TRIANGULAR / NO TRANSPOSE*****/
			if (trans == 'N' || trans == 'n') {

				for (int i=n-1; i>=0; --i) {

					for (int j=i+1; j<n; ++j) {
						for (int k=0; k<m; ++k) {			//For each system

							x[i*ldx +k] -= a[i*lda +j] * x[j*ldx +k];
						
						}
					}

					for (int k=0; k<m; ++k) {			//For each system
					
						x[i*ldx +k] /= a[i*lda +i];
					
					}

				}

			}

			/*****UPPER TRIANGULAR / NOT UNIT TRIANGULAR / TRANSPOSE*****/
			else if (trans == 'T' || trans == 't') {

				for (int i=0; i<n; ++i) {

					for (int j=0; j<i; ++j) {
						for (int k=0; k<m; ++k) {			//For each system

							x[i*ldx +k] -= a[j*lda +i] * x[j*ldx +k];
						
						}
					}

					for (int k=0; k<m; ++k) {			//For each system
				
						x[i*ldx +k] /= a[i*lda +i];
					
					}

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

				for (int i=0; i<n; ++i) {
					for (int j=0; j<i; ++j) {
						for (int k=0; k<m; ++k) {			//For each system

							x[i*ldx +k] -= a[i*lda +j] * x[j*ldx +k];

						}
					}
				}

			}
	
			/*****LOWER TRIANGULAR / UNIT TRIANGULAR / TRANSPOSE*****/
			else if (trans == 'T' || trans == 't') {

				for (int i=n-1; i>=0; --i) {
					for (int j=i+1; j<n; ++j) {
						for (int k=0; k<m; ++k) {			//For each system

							x[i*ldx +k] -= a[j*lda +i] * x[j*ldx +k];

						}
					}
				}

			}
		
		}	//UNIT TRIANGULAR
	
		/*****LOWER TRIANGULAR / NOT UNIT TRIANGULAR*****/
		else if (diag == 'N' || diag == 'n') {
	
			/*****LOWER TRIANGULAR / NOT UNIT TRIANGULAR / NO TRANSPOSE*****/
			if (trans == 'N' || trans == 'n') {

				for (int i=0; i<n; ++i) {

					for (int j=0; j<i; ++j) {
						for (int k=0; k<m; ++k) {			//For each system
			
							x[i*ldx +k] -= a[i*lda +j] * x[j*ldx +k];
						
						}
					}

					for (int k=0; k<m; ++k) {			//For each system
				
						x[i*ldx +k] /= a[i*lda +i];
					
					}

				}

			}
	
			/*****LOWER TRIANGULAR / NOT UNIT TRIANGULAR / TRANSPOSE*****/
			else if (trans == 'T' || trans == 't') {
	
				for (int i=n-1; i>=0; --i) {

					for (int j=i+1; j<n; ++j) {
						for (int k=0; k<m; ++k) {			//For each system

							x[i*ldx +k] -= a[j*lda +i] * x[j*ldx +k];

						}
					}

					for (int k=0; k<m; ++k) {			//For each system
				
						x[i*ldx +k] /= a[i*lda +i];
					
					}

				}

			}

		}	//NOT UNIT TRIANGULAR

	}	//LOWER TRIANGULAR

}	//AXB

/****************************************/
/****************************************/
/****************************************/

/***** XA = B ORDERING *****/
else if (ordering == 'X' || ordering == 'x') {

	/*****UPPER TRIANGULAR*****/
	if (uplo == 'U' || uplo == 'u') {

		/*****UPPER TRIANGULAR / UNIT TRIANGULAR*****/
		if (diag == 'U' || diag == 'u') {

			/*****UPPER TRIANGULAR / UNIT TRIANGULAR / NO TRANSPOSE*****/
			if (trans == 'T' || trans == 't') {

				for (int i=n-1; i>=0; --i) {
					for (int j=i+1; j<n; ++j) {
						for (int k=0; k<m; ++k) {			//For each system

							x[k*ldx +i] -= a[i*lda +j] * x[k*ldx +j];
						
						}
					}
				}

			}

			/*****UPPER TRIANGULAR / UNIT TRIANGULAR / TRANSPOSE*****/
			else if (trans == 'N' || trans == 'n') {

				for (int i=0; i<n; ++i) {
					for (int j=0; j<i; ++j) {
						for (int k=0; k<m; ++k) {			//For each system

							x[k*ldx +i] -= a[j*lda +i] * x[k*ldx +j];
						
						}

					}
				}

			}

		}	//UNIT TRIANGULAR

		/*****UPPER TRIANGULAR / NOT UNIT TRIANGULAR*****/
		else if (diag == 'N' || diag == 'n') {

			/*****UPPER TRIANGULAR / NOT UNIT TRIANGULAR / NO TRANSPOSE*****/
			if (trans == 'T' || trans == 't') {
	
				for (int i=n-1; i>=0; --i) {

					for (int j=i+1; j<n; ++j) {
						for (int k=0; k<m; ++k) {			//For each system

							x[k*ldx +i] -= a[i*lda +j] * x[k*ldx +j];
						
						}
					}

					for (int k=0; k<m; ++k) {			//For each system
					
						x[k*ldx +i] /= a[i*lda +i];
					
					}

				}

			}

			/*****UPPER TRIANGULAR / NOT UNIT TRIANGULAR / TRANSPOSE*****/
			else if (trans == 'N' || trans == 'n') {

				for (int i=0; i<n; ++i) {

					for (int j=0; j<i; ++j) {
						for (int k=0; k<m; ++k) {			//For each system

							x[k*ldx +i] -= a[j*lda +i] * x[k*ldx +j];
						
						}
					}

					for (int k=0; k<m; ++k) {			//For each system
				
						x[k*ldx +i] /= a[i*lda +i];
					
					}

				}

			}

		}	//NOT UNIT TRIANGULAR

	}	//UPPER TRIANGULAR

	/*****LOWER TRIANGULAR*****/
	else if (uplo == 'L' || uplo == 'l') {

		/*****LOWER TRIANGULAR / UNIT TRIANGULAR*****/
		if (diag == 'U' || diag == 'u') {
			/*****LOWER TRIANGULAR / UNIT TRIANGULAR / NO TRANSPOSE*****/
			if (trans == 'T' || trans == 't') {

				for (int i=0; i<n; ++i) {
					for (int j=0; j<i; ++j) {
						for (int k=0; k<m; ++k) {			//For each system

							x[k*ldx +i] -= a[i*lda +j] * x[k*ldx +j];

						}
					}
				}

			}
	
			/*****LOWER TRIANGULAR / UNIT TRIANGULAR / TRANSPOSE*****/
			else if (trans == 'N' || trans == 'n') {

				for (int i=n-1; i>=0; --i) {
					for (int j=i+1; j<n; ++j) {
						for (int k=0; k<m; ++k) {			//For each system

							x[k*ldx +i] -= a[j*lda +i] * x[k*ldx +j];

						}
					}
				}

			}
		
		}	//UNIT TRIANGULAR
	
		/*****LOWER TRIANGULAR / NOT UNIT TRIANGULAR*****/
		else if (diag == 'N' || diag == 'n') {
	
			/*****LOWER TRIANGULAR / NOT UNIT TRIANGULAR / NO TRANSPOSE*****/
			if (trans == 'T' || trans == 't') {

				for (int i=0; i<n; ++i) {

					for (int j=0; j<i; ++j) {
						for (int k=0; k<m; ++k) {			//For each system
			
							x[k*ldx +i] -= a[i*lda +j] * x[k*ldx +j];
						
						}
					}

					for (int k=0; k<m; ++k) {			//For each system
				
						x[k*ldx +i] /= a[i*lda +i];
					
					}

				}

			}
	
			/*****LOWER TRIANGULAR / NOT UNIT TRIANGULAR / TRANSPOSE*****/
			else if (trans == 'N' || trans == 'n') {
	
				for (int i=n-1; i>=0; --i) {

					for (int j=i+1; j<n; ++j) {
						for (int k=0; k<m; ++k) {			//For each system

							x[k*ldx +i] -= a[j*lda +i] * x[k*ldx +j];

						}
					}

					for (int k=0; k<m; ++k) {			//For each system
				
						x[k*ldx +i] /= a[i*lda +i];
					
					}

				}

			}

		}	//NOT UNIT TRIANGULAR

	}	//LOWER TRIANGULAR

}	//XAB

}


#endif	// LU_DECOMP_TRSV_BLK_TCC
