/////////////////////////////////////////////////////////////////////
//		BLAS
/////////////////////////////////////////////////////////////////////
extern "C" void dgemv_(char *TRANS,int *M, int *N, double *ALPHA, double *A, int *LDA, double *X, int *INCX, double *BETA,double *Y, int *INCY);

extern "C" double dnrm2_(int *N,double *X, int *INCX);
/////////////////////////////////////////////////////////////////////
//		LAPACK
/////////////////////////////////////////////////////////////////////
extern "C" void dsyevd_(char *JOBZ,char *UPLO,int *N,double *A,int *LDA,double *W
		,double *WORK,int *LWORK
		,int *IWORK,int *LIWORK,int *INFO );

// DGELSY computes the minimum-norm solution to a real linear least
// squares problem: minimize || A * X - B ||
extern "C" void dgelsy_(int *M, int *N, int *NRHS, double  *A, int *LDA,
		double *B, int *LDB,int *JPVT, double *RCOND, int *RANK,
		double *WORK, int *LWORK, int *INFO );
// DGELS solves overdetermined or underdetermined real linear systems
// involving an M-by-N matrix A, or its transpose, using a QR or LQ
// factorization of A.  It is assumed that A has full rank.
extern "C" void dgels_(char *TRANS, int *M,int *N, int *NRHS,
		double *A, int *LDA, double *B,int *LDB,
		double *WORK,int *LWORK,int *INFO);
