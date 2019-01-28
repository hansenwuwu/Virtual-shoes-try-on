#ifndef _LINEAR_H_
#define _LINEAR_H_

//
// linear algebra function
//


double det( double *, int, void *memory = 0x00 );
//int inverse( double *, double *, int n, void *memory = 0x00 );	// return 1 if it's a singular matrix
int inverse( double *, int n, void *memory = 0x00 );

// Gauss-Jordan Elimination
int GaussJordan( double *A, double *b, int n );

// Gram-Schmidt Process
void Gram_Schmidt( double *A, double *Q, int row, int col );

// Singular Value Decomposition	
void SVD( double *A, double *U, double *W, double *V, int row, int col, void *memory = 0x00 );
void svd( double *x, double *u, double *w, double *v, int row, int col, void *memory = 0x00 );

// QR factorization
void QRfact( double *A, double *Q, double *R, int row, int col );

// LU factorization
int LUfact( double *A, int *index, int n, double *mem = 0x00 );
void SolveLU( double *LU, int *index, double *b, int n );
void LUcomb( double *LU, int *index, int n );

// Least-Square Method
void LS_method( double *A, double *B, double *X, int row, int col, void *memory = 0x00 );

// eigensystem
void Jacobi( double *A, double *eigen, double *v, int n, int time = 30 );
double PowerIter( double *A, double *v, int n, int time = 30, double *mem = 0x00 );
double InvIter( double *A, double *v, double shift, int n, int time = 10, double *mem = 0x00 );
void ReductHessen( double *A, int n );	// Reduction to Hessenberg Form
void HessenQR( double *A, int n, double *w, int time = 30 );	// QR for real matrix only

//
// non-linear optimization
//

double levmar( void( *f )( double*, double*, int, int ),
			  double *x, double *p, int nx, int np, int time,
			  double tau, double, double, double, void *memory = 0x00 );
// f( x[], p[], nx, np )

double Newton( void( *f )( double *, int, double * ),
			 double x[], int n, double tolx, double tolf, 
			 int time = 30, void *memory = 0x00 );
// f( x[], n, fvector[] )


//
// Partial Differential Equations
//

void MultiGrid( double *u, int n, int cycle );
void Relax( double *, double *, double *, int, double );
int FMG( double *f, int n, int a1, int a2, int ncycle );
int MG( double *v, double *f, int n, int a1, int a2, int ncycle );

//
// C-style function
//



#endif