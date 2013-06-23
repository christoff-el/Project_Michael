#ifndef ERR_CALC_H
#define ERR_CALC_H 1

template <typename T>
T
errCalc(int n, T *exact, T *calc, int incX, int prec) {

	mpf_set_default_prec(prec);
	
	T error = 0;
	
	int firstInxX = (copysign(1,incX) < 0) ? (n-1)*(-incX) : 0;
	
	for (int i=0, iX=firstInxX; i<n; ++i, iX+=incX) {
	
		error += abs(exact[i] - calc[iX]);
		
	}
	
	return error;
		
}

template <typename T>
T
errCalc2d(int n, int m, T *exact, T *calc, int prec) {

	mpf_set_default_prec(prec);
	
	T error = 0;
	
	for (int i=0; i<n; ++i) {
		for (int j=0; j<m; ++j) {
		
			error += abs(exact[i*m +j] - calc[i*m +j]);
			
		}
	}
	
	return error;
	
}

#endif	//ERR_CALC_H