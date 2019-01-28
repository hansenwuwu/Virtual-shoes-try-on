
//
// Same as linear.h, but programming by c, not c++
//

#ifndef _LINEAR_2_H_
#define _LINEAR_2_H_

#ifdef __cplusplus
extern "C"{
#endif

//
// Type Macro
//

typedef int			Row, Col, Dim, Row_in, IterTime;
typedef float		Paramf;
typedef double		Paramd;
typedef float		*Matf, *Matf_in, *Matf_out;
typedef double		*Matd, *Matd_in, *Matd_out;
typedef int			*Veci, *Veci_in, *Veci_out;
typedef float		*Vecf, *Vecf_in, *Vecf_out;
typedef double		*Vecd, *Vecd_in, *Vecd_out;

typedef Row			Row_in, Row_out;
typedef Col			Col_in, Col_out;
typedef Vecf		*Mat2f, *Mat2f_in, *Mat2f_out;
typedef Vecd		*Mat2d, *Mat2d_in, *Mat2d_out;

typedef float		Errf;
typedef double		Errd;

// determinant
float detf( Matf, Dim );
double detd( Matd, Dim );

// inverse
int invf( Matf_out inv, Matf_in, Dim );
int invd( Matd_out inv, Matd_in, Dim );

// psudo inverse
int pinvf( Matf_out inv, Matf_in, Row, Col, Errf );
int pinvd( Matd_out inv, Matd_in, Row, Col, Errd );

// Gauss-Jordan Elimination
//int GaussJordan( double *A, double *b, int n );

// Gram-Schmidt Process
//void Gram_Schmidt( double *A, double *Q, int row, int col );

// Singular Value Decomposition	
void svdf( Matf_in x, Matf_out u, Vecf_out w, Matf_out vT, Row, Col, Errf, Matf xx );
void svdd( Matd_in x, Matd_out u, Vecd_out w, Matd_out vT, Row, Col, Errd, Matd xx );

// QR factorization
//void QRfact( double *A, double *Q, double *R, int row, int col );

// LU factorization
int LUfactf( Matf A, Veci i, Dim n, Vecf v );		// A will be upper + lower triangle matrix, return v as pivot vector, and index i
int LUfactd( Matd A, Veci i, Dim n, Vecd v );
void SolveLUf( Matf LU, Veci i, Vecf b, Dim n );	// solve linear equation by factorized matrix lu and row-interchange index
void SolveLUd( Matd LU, Veci i, Vecd b, Dim n );
void LUcombf( Matf LU, Veci i, Dim n );				// compute origin matrix from factorized matrix and row-interchange index
void LUcombd( Matd LU, Veci i, Dim n );

// Least-Square Method
//void LS_method( double *A, double *B, double *X, int row, int col, void *memory );

// eigensystem
void Jacobif( Matf a, Vecf eigen, Matf eig_vec, Dim n, int time );
void Jacobid( Matd a, Vecd eigne, Matd eig_vec, Dim n, int time );
float PowerIterf( Matf a, Vecf eigen, Dim n, int time, Errf, Vecf );
double PowerIterd( Matd a, Vecd eigen, Dim n, int time, Errd, Vecd );

// inverse iteration, input matrix A and shift value, return v as eigehvector, need pivot vector b1 and row-interchange index ind for LU factorization
float InvIterf( Matf A, Vecf v, float shift, Dim n, int time, Errf, Vecf b1, Veci i );	
double InvIterd( Matd A, Vecd v, double shift, Dim n, int time, Errd, Vecd b1, Veci i );

//void ReductHessen( double *A, int n );	// Reduction to Hessenberg Form
//void HessenQR( double *A, int n, double *w, int time );	// QR for real matrix only

#ifdef __cplusplus
}
#endif

#endif