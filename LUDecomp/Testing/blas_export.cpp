//#include <iostream>

//#include "../level2/gbmv.tcc"			<-- NA
//#include "../level2/gemv.tcc"
//#include "../level2/ger.tcc"
//#include "../level2/sbmv.tcc"			<-- NA
//#include "../level2/spmv.tcc"
//#include "../level2/spr2.tcc"
//#include "../level2/spr.tcc"
//#include "../level2/symv.tcc"			<-- NA
//#include "../level2/syr2.tcc"			<-- NA
//#include "../level2/syr.tcc"			<-- NA
//#include "../level2/tbmv.tcc"			<-- NA
//#include "../level2/tbsv.tcc"			<-- NA
//#include "../level2/tpmv.tcc"			<-- NA
//#include "../level2/tpsv.tcc"			<-- NA
//#include "../level2/trmv.tcc"
#include "../trsv.tcc"


extern "C" {

//-----------------------------------------------------------------//
//----------------------------Level 2------------------------------//
//-----------------------------------------------------------------//

//-------Single Precision-------//

//sgbmv									<-- Not available
void
(sgbmv_)(char			*TRANS,
		 int			*M,
		 int			*N,
		 int			*KL,
		 int			*KU,
		 const float	*ALPHA,
		 const float	*A,
		 int			*LDA,
		 const float	*X,
		 int			*INCX,
		 const float	*BETA,
		 float			*Y,
		 int			*INCY)
{
	//gbmv(*TRANS,*M,*N,*KL,*KU,*ALPHA,A,*LDA,X,*INCX,*BETA,Y,*INCY);
}

//sgemv
void
(sgemv_)(char			*TRANS,
		 int			*M,
		 int			*N,
		 const float	*ALPHA,
		 const float	*A,
		 int			*LDA,
		 const float	*X,
		 int			*INCX,
		 const float	*BETA,
		 float			*Y,
		 int			*INCY)
{
	//gemv(*TRANS,*M,*N,*ALPHA,A,*LDA,X,*INCX,*BETA,Y,*INCY);
}

//ger
void
(sger_)(int				*M,
		int				*N,
		const float		*ALPHA,
		const float		*X,
		int				*INCX,
		const float		*Y,
		int				*INCY,
		float			*A,
		int				*LDA)
{
	//ger(*M,*N,*ALPHA,X,*INCX,Y,*INCY,A,*LDA);
}

//sbmv
void
(ssbmv_)(char			*UPLO,
		 int			*N,
		 int			*K,
		 const float	*ALPHA,
		 const float	*A,
		 int			*LDA,
		 const float	*X,
		 int			*INCX,
		 const float	*BETA,
		 float			*Y,
		 int			*INCY)
{
	//sbmv(*UPLO,*N,*K,ALPHA,*LDA,X,*INCX,*BETA,Y,*INCY);
}

//spmv
void
(sspmv_)(char			*UPLO,
		 int			*N,
		 const float	*ALPHA,
		 const float	*AP,
		 const float	*X,
		 int			*INCX,
		 const float	*BETA,
		 float			*Y,
		 int			*INCY)
{
	//spmv(*UPLO,*N,*ALPHA,AP,X,*INCX,*BETA,Y,*INCY);
}

//spr2
void
(sspr2_)(char			*UPLO,
		 int			*N,
		 const float	*ALPHA,
		 const float	*X,
		 int			*INCX,
		 const float	*Y,
		 int			*INCY,
		 float			*AP)
{
	//spr2(*UPLO,*N,*ALPHA,X,*INCX,Y,*INCY,AP);
}

//spr
void
(sspr_)(char			*UPLO,
		int				*N,
		const float		*ALPHA,
		const float		*X,
		int				*INCX,
		float			*AP)
{
	//spr(*UPLO,*N,*ALPHA,X,*INCX,AP);
}

//symv
void
(ssymv_)(char			*UPLO,
		 int			*N,
		 const float	*ALPHA,
		 const float	*A,
		 int			*LDA,
		 const float	*X,
		 int			*INCX,
		 const float	*BETA,
		 float			*Y,
		 int			*INCY)
{
	//symv(*UPLO,*N,*ALPHA,A,*LDA,X,*INCX,*BETA,Y,*INCY);
}

//syr2
void
(ssyr2_)(char			*UPLO,
		 int			*N,
		 const float	*ALPHA,
		 const float	*X,
		 int			*INCX,
		 const float	*Y,
		 int			*INCY,
		 float			*A,
		 int			*LDA)
{
	//syr2(*UPLO,*N,*ALPHA,X,*INCX,Y,*INCY,A,*LDA);
}
		 
//syr
void
(ssyr_)(char			*UPLO,
		int				*N,
		const float		*ALPHA,
		const float		*X,
		int				*INCX,
		float			*A,
		int				*LDA)
{
	//syr(*UPLO,*N,*ALPHA,X,*INCX,A,*LDA);
}

//tbmv
void
(stbmv_)(char			*UPLO,
		 char			*TRANS,
		 char			*DIAG,
		 int			*N,
		 int			*K,
		 const float	*A,
		 int			*LDA,
		 float			*X,
		 int			*INCX)
{
	//tbmv(*UPLO,*TRANS,*DIAG,*N,*K,A,*LDA,X,*INCX);
}

//tbsv
void
(stbsv_)(char			*UPLO,
		 char			*TRANS,
		 char			*DIAG,
		 int			*N,
		 int			*K,
		 const float	*A,
		 int			*LDA,
		 float			*X,
		 int			*INCX)
{
	//tbsv(*UPLO,*TRANS,*DIAG,*N,*K,A,*LDA,X,*INCX);
}

//tpmv
void
(stpmv_)(char			*UPLO,
		 char			*TRANS,
		 char			*DIAG,
		 int			*N,
		 const float	*AP,
		 float			*X,
		 int			*INCX)
{
	//tpmv(*UPLO,*TRANS,*DIAG,*N,AP,X,*INCX);
}

//tpsv
void
(stpsv_)(char			*UPLO,
		 char			*TRANS,
		 char			*DIAG,
		 int			*N,
		 const float	*AP,
		 float			*X,
		 int			*INCX)
{
	//tpsv(*UPLO,*TRANS,*DIAG,*N,AP,X,*INCX);
}

//trmv
void
(strmv_)(char			*UPLO,
		 char			*TRANS,
		 char			*DIAG,
		 int			*N,
		 const float	*A,
		 int			*LDA,
		 float			*X,
		 int			*INCX)
{
	//trmv(*UPLO,*TRANS,*DIAG,*N,A,*LDA,X,*INCX);
}

//trsv
void
(strsv_)(char			*UPLO,
		 char			*TRANS,
		 char			*DIAG,
		 int			*N,
		 const float	*A,
		 int			*LDA,
		 float			*X,
		 int			*INCX)
{
	//trsv(*UPLO,*TRANS,*DIAG,*N,A,*LDA,X,*INCX);
}


//-------Double Precision-------//

//dgbmv									<-- Not available
void
(dgbmv_)(char			*TRANS,
		 int			*M,
		 int			*N,
		 int			*KL,
		 int			*KU,
		 const double	*ALPHA,
		 const double	*A,
		 int			*LDA,
		 const double	*X,
		 int			*INCX,
		 const double	*BETA,
		 double			*Y,
		 int			*INCY)
{
	//gbmv(*TRANS,*M,*N,*KL,*KU,*ALPHA,A,*LDA,X,*INCX,*BETA,Y,*INCY);
}

//dgemv
void
(dgemv_)(char			*TRANS,
		 int			*M,
		 int			*N,
		 const double	*ALPHA,
		 const double	*A,
		 int			*LDA,
		 const double	*X,
		 int			*INCX,
		 const double	*BETA,
		 double			*Y,
		 int			*INCY)
{
	//gemv(*TRANS,*M,*N,*ALPHA,A,*LDA,X,*INCX,*BETA,Y,*INCY);
}

//ger
void
(dger_)(int				*M,
		int				*N,
		const double	*ALPHA,
		const double	*X,
		int				*INCX,
		const double	*Y,
		int				*INCY,
		double			*A,
		int				*LDA)
{
	//ger(*M,*N,*ALPHA,X,*INCX,Y,*INCY,A,*LDA);
}

//sbmv
void
(dsbmv_)(char			*UPLO,
		 int			*N,
		 int			*K,
		 const double	*ALPHA,
		 const double	*A,
		 int			*LDA,
		 const double	*X,
		 int			*INCX,
		 const double	*BETA,
		 double			*Y,
		 int			*INCY)
{
	//sbmv(*UPLO,*N,*K,ALPHA,*LDA,X,*INCX,*BETA,Y,*INCY);
}

//spmv
void
(dspmv_)(char			*UPLO,
		 int			*N,
		 const double	*ALPHA,
		 const double	*AP,
		 const double	*X,
		 int			*INCX,
		 const double	*BETA,
		 double			*Y,
		 int			*INCY)
{
	//spmv(*UPLO,*N,*ALPHA,AP,X,*INCX,*BETA,Y,*INCY);
}

//spr2
void
(dspr2_)(char			*UPLO,
		 int			*N,
		 const double	*ALPHA,
		 const double	*X,
		 int			*INCX,
		 const double	*Y,
		 int			*INCY,
		 double			*AP)
{
	//spr2(*UPLO,*N,*ALPHA,X,*INCX,Y,*INCY,AP);
}

//spr
void
(dspr_)(char			*UPLO,
		int				*N,
		const double	*ALPHA,
		const double	*X,
		int				*INCX,
		double			*AP)
{
	//spr(*UPLO,*N,*ALPHA,X,*INCX,AP);
}

//symv
void
(dsymv_)(char			*UPLO,
		 int			*N,
		 const double	*ALPHA,
		 const double	*A,
		 int			*LDA,
		 const double	*X,
		 int			*INCX,
		 const double	*BETA,
		 double			*Y,
		 int			*INCY)
{
	//symv(*UPLO,*N,*ALPHA,A,*LDA,X,*INCX,*BETA,Y,*INCY);
}

//syr2
void
(dsyr2_)(char			*UPLO,
		 int			*N,
		 const double	*ALPHA,
		 const double	*X,
		 int			*INCX,
		 const double	*Y,
		 int			*INCY,
		 double			*A,
		 int			*LDA)
{
	//syr2(*UPLO,*N,*ALPHA,X,*INCX,Y,*INCY,A,*LDA);
}
		 
//syr
void
(dsyr_)(char			*UPLO,
		int				*N,
		const double	*ALPHA,
		const double	*X,
		int				*INCX,
		double			*A,
		int				*LDA)
{
	//syr(*UPLO,*N,*ALPHA,X,*INCX,A,*LDA);
}

//tbmv
void
(dtbmv_)(char			*UPLO,
		 char			*TRANS,
		 char			*DIAG,
		 int			*N,
		 int			*K,
		 const double	*A,
		 int			*LDA,
		 double			*X,
		 int			*INCX)
{
	//tbmv(*UPLO,*TRANS,*DIAG,*N,*K,A,*LDA,X,*INCX);
}

//tbsv
void
(dtbsv_)(char			*UPLO,
		 char			*TRANS,
		 char			*DIAG,
		 int			*N,
		 int			*K,
		 const double	*A,
		 int			*LDA,
		 double			*X,
		 int			*INCX)
{
	//tbsv(*UPLO,*TRANS,*DIAG,*N,*K,A,*LDA,X,*INCX);
}

//tpmv
void
(dtpmv_)(char			*UPLO,
		 char			*TRANS,
		 char			*DIAG,
		 int			*N,
		 const double	*AP,
		 double			*X,
		 int			*INCX)
{
	//tpmv(*UPLO,*TRANS,*DIAG,*N,AP,X,*INCX);
}

//tpsv
void
(dtpsv_)(char			*UPLO,
		 char			*TRANS,
		 char			*DIAG,
		 int			*N,
		 const double	*AP,
		 double			*X,
		 int			*INCX)
{
	//tpsv(*UPLO,*TRANS,*DIAG,*N,AP,X,*INCX);
}

//trmv
void
(dtrmv_)(char			*UPLO,
		 char			*TRANS,
		 char			*DIAG,
		 int			*N,
		 const double	*A,
		 int			*LDA,
		 double			*X,
		 int			*INCX)
{
	//trmv(*UPLO,*TRANS,*DIAG,*N,A,*LDA,X,*INCX);
}

//trsv
void
(dtrsv_)(char			*UPLO,
		 char			*TRANS,
		 char			*DIAG,
		 int			*N,
		 double			*A,
		 int			*LDA,
		 double			*X,
		 int			*INCX)
{

	// ! swapping FORTRAN col-major to c++ row-major:
	double Acopy[(*LDA)*(*N)];		//Carbon copy of A to stop it whinging
	double tmp[(*LDA)*(*N)];		//Preserves LDA
	
	//Make copy of A:
	for (int i=0; i<(*N)*(*LDA); ++i) {
	
		Acopy[i] = A[i];
		
	}
	
	//To row-major:
	for (int i=0; i<(*N); ++i) {
		for (int j=0; j<(*N); ++j) {
	
			tmp[i*(*LDA) +j] = A[j*(*LDA) +i];
			
		}
	}
	
	//Copy to A:
	for (int i=0; i<(*N); ++i) {
		for (int j=0; j<(*N); ++j) {
	
			A[i*(*LDA) +j] = tmp[i*(*LDA) +j];
		
		}
	}
	
	//Do some testings:
	trsv(*UPLO,*TRANS,*DIAG,*N,A,*LDA,X,*INCX);
	
	//Restore A:
	for (int i=0; i<(*N)*(*LDA); ++i) {
	
		A[i] = Acopy[i];
		
	}
}







}	//extern "C"