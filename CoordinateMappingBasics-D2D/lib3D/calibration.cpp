
#include <windows.h>
#include "calibration.h"
#include "numerical\linear.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "include\version.h"

//
// Macro
//

#define ERR 0.000000000000001
#ifndef SQUARE
#define SQUARE(x) ( (x) * (x) )
#endif

namespace LevMar{

	// define
	#define alpha	p[0]
	#define gamma	p[1]
	#define u0		p[2]
	#define beta	p[3]
	#define v0		p[4]
	#define k1		p[5]
	#define k2		p[6]

	// variables
	double *M, A[9], Rt[12], x[4], temp[3];
	int point, i, j, k, n;

	// homography
	void hmg( double *m, double *p, int nm, int np ){
		for( i = 0; i < nm; i += 2 ){
			m[i] = ( p[0] * M[i] + p[1] * M[i + 1] + p[2] );
			m[i + 1] = ( p[3] * M[i] + p[4] * M[i + 1] + p[5] );
			m[i] /= ( p[6] * M[i] + p[7] * M[i + 1] + p[8] );
			m[i + 1] /= ( p[6] * M[i] + p[7] * M[i + 1] + p[8] );
		}
	}

	// calibration
	void clb( double *m, double *p, int nm, int np ){
		// initialize A
		A[0] = alpha;
		A[1] = gamma;
		A[2] = u0;
		A[4] = beta;
		A[5] = v0;
		A[3] = A[6] = A[7] = 0.0;
		A[8] = 1.0;

		for( n = 0, j = 5; n < np; n += 2 ){
			// set x
			x[0] = ( M + n )[0];
			x[1] = ( M + n )[1];
			x[2] = 0.0; 
			x[3] = 1.0;
			// next image
			if( ( n / 2 ) % point == 0 ){
				for( k = 0; k < 12; k++, j++ ) Rt[k] = p[j];
			}
			// x = A * Rt * x
			for( i = 0; i < 3; i++ ){
				temp[i] = 0.0;
				for( j = 0; j < 4; j++ ) temp[i] += Rt[i * 4 + j] * x[j];
			}
			for( i = 0; i < 3; i++ ){
				x[i] = 0.0;
				for( j = 0; j < 3; j++ ) x[i] += A[i * 3 + j] * temp[j];
			}
			// m[..., i, i+1, ...]
			m[i] = x[0] / x[2];
			m[i + 1] = x[1] / x[2];
		}
	}

	// Distortion
	void dist( double *m, double *p, int nm, int np ){
		for( k = 0, j = 7; k < np; k += 2 ){
			// set x
			x[0] = ( M + k )[0];
			x[1] = ( M + k )[1];
			x[2] = 0.0; 
			x[3] = 1.0;
			// next image, set Rt
			if( ( n / 2 ) % point == 0 ){
				for( i = 0; i < 12; i++, j++ ) Rt[i] = p[j];
			}
			// Normalize, temp = Rt * x
			for( i = 0; i < 3; i++ ){
				temp[i] = 0.0;
				for( j = 0; j < 4; j++ ) temp[i] += Rt[i * 4 + j] * x[j];
			}
			x[0] = temp[0] / temp[2];
			x[1] = temp[1] / temp[2];
			// distortion model
			x[2] = SQUARE( x[0] ) * SQUARE( x[1] );
			x[0] *= ( 1.0 + k1 * x[2] + k2 * SQUARE( x[2] ) );
			x[1] *= ( 1.0 + k1 * x[2] + k2 * SQUARE( x[2] ) );
			x[2] = 1.0;
			// intrinsic trnasform: m = A * x
			m[i] = x[0] * alpha + x[1] * gamma + u0;
			m[i + 1] = x[1] * beta + v0;
		}
	}

	#undef alpha
	#undef beta
	#undef gamma
	#undef u0
	#undef v0
	#undef k1
	#undef k2
};

// 
// CEMERA_PARAM
//

// Constructor

CAMERA_PARAM::CAMERA_PARAM( void ){
	_H = NULL;
	_n_img = _n_pt = 0;
	_m_size = 0;
	_mem = NULL;
}

CAMERA_PARAM::CAMERA_PARAM( int img, int pt, int w, int h ){
	_H = ( double* )malloc( img * 10 * sizeof( double ) );
	_m_size = 0;
	_mem = NULL;
	_n_img=img;
	_n_pt=pt;
	_h = h;
	_w = w;
}

CAMERA_PARAM::CAMERA_PARAM( char *filename ){
	_H = NULL;
	_n_img = _n_pt = 0;
	_m_size = 0;
	_mem = NULL;
	ReadParam( filename );
}

// Destructor

CAMERA_PARAM::~CAMERA_PARAM(void){
	Clear();
}

// Functions

bool CAMERA_PARAM::Calib( Point3d *m, Point3d *M ){
	// check is this object created for read parameter file
	if( _H == NULL ) return false;

	double *n = ( double* )malloc( _n_pt * _n_img * 4 * sizeof( double ) );
	double *N = n + _n_pt * _n_img * 2;
	int i;
	for( i = 0; i < _n_pt * _n_img; i++ ){
		n[i * 2 + 0] = m[i].x;
		n[i * 2 + 1] = m[i].y;
		N[i * 2 + 0] = M[i].x;
		N[i * 2 + 1] = M[i].y;
	}
	if( Calib( n, N ) ){
		free( n );
		return true;
	}
	free( n );
	return false;
}

bool CAMERA_PARAM::Calib( double *m, double *M ){
	// check is this object created for read parameter file
	if( _H == NULL ) return false;

	double *param;
	int i;

	if( !SolveHomo( m, M ) ) return false;
	if( !SolveCalib( m, M ) ) return false;
//	if( !SolveDistor( m, M ) ) return false;
	return true;

	// Complete Maximum Likelihood Estimation
	LevMar::M = M;				// set function
	LevMar::point = _n_pt;		// set point per image

	// memory pool
	i = ( 7 + 12 * _n_img ) * sizeof( double );
	if( ( int )_m_size < i ){
		_m_size = ( unsigned int )i;
		_mem = realloc( _mem, _m_size );
	}
	param = ( double* )_mem;					// initialize param
	param[0] = _alpha;							// the top 5 intrinsic parameters
	param[1] = _gamma;
	param[2] = _u0;
	param[3] = _beta;
	param[4] = _v0;
	param[5] = _k1;
	param[6] = _k2;
	for( i = 0; i < _n_img; i++ ) {					// initialize extrinsic parameters
		GetExtrinsic( param + 7 + 12 * i, i );		// for each image
	}

	// levmar algorithm
	levmar( LevMar::dist, m, param, _n_pt * _n_img * 2, 7 + _n_img * 12, 100, 0.001, ERR, ERR, ERR );

	// updata the parameter
	_alpha = param[0];
	_gamma = param[1];
	_u0 = param[2];
	_beta = param[3];
	_v0 = param[4];
	_k1 = param[5];
	_k2 = param[6];
	GetExtrinsic( _Rt, 0 );

	return true;
}

bool CAMERA_PARAM::SolveHomo( Point3d *m, Point3d *M){
	// check is this object created for read parameter file
	if( _H == NULL ) return false;
	
	double *n = ( double* )malloc( _n_pt * _n_img * 4 * sizeof( double ) );
	double *N = n + _n_pt * _n_img * 2;
	int i;
	for( i = 0; i < _n_pt * _n_img; i++ ){
		n[i * 2 + 0] = m[i].x;
		n[i * 2 + 1] = m[i].y;
		N[i * 2 + 0] = M[i].x;
		N[i * 2 + 1] = M[i].y;
	}
	if( !SolveHomo(n, N) ){
		free( n );
		return true;
	}
	free( n );
	return false;
}

#define l( i, j ) L[(i) * 9 + j]
bool CAMERA_PARAM::SolveHomo(double* m, double* M){
	// check is this object created for read parameter file
	if( _H == NULL ) return false;

	double min, *L, *A, *eig_vec, *eig_value, *x;
	int i, j, k, a, flag;

	// memory pool
	i = ( 18 * _n_pt + 180 ) * sizeof( double );
	if( ( int )_m_size < i ){
		_m_size = ( unsigned int )i;
		_mem = realloc( _mem, _m_size );
	}
	L = ( double* )_mem;			// 2_n_pt * 9
	eig_vec = L + 2 * _n_pt * 9;	// 9 * 9
	eig_value = eig_vec + 9 * 9;	// 9
	x = eig_value + 9;				// 9
	A = x + 9;						// 81
	
	for( a = 0; a < _n_img; a++ ){
		// [Mt	0	-uMt]
		// |			|x=0
		// [0	Mt	-vMt]
		for( i = 0; i < _n_pt; i++ ){
			flag = ( a * _n_pt + i ) * 2;			// the point's position in matrix
			l( i * 2, 3 ) = l( i * 2, 4 ) = l( i * 2, 5 ) =
			l( i * 2 + 1, 0 ) = l( i * 2 + 1, 1 ) = l( i * 2 + 1, 2 ) = 0.0;

			l( i * 2, 0 ) = l( i * 2 + 1, 3 ) = M[flag + 0];
			l( i * 2, 1 ) = l( i * 2 + 1, 4 ) = M[flag + 1];
			l( i * 2, 2 ) = l( i * 2 + 1, 5 ) = 1.0;

			// -uMt
			l( i * 2, 6 ) = -m[flag + 0] * M[flag + 0];
			l( i * 2, 7 ) = -m[flag + 0] * M[flag + 1];
			l( i * 2, 8 ) = -m[flag + 0];
			
			// -vMt
			l( i * 2 + 1, 6 ) = -m[flag + 1] * M[flag + 0];
			l( i * 2 + 1, 7 ) = -m[flag + 1] * M[flag + 1];
			l( i * 2 + 1, 8 ) = -m[flag + 1];
		}
		// generate A = LT * L
		flag = 2 * _n_pt;
		for( i = 0; i < 9; i++ ){
			for( j = i; j < 9; j++ ) {
				A[i * 9 + j] = 0.0;
				for( k = 0; k < flag; k++ ) {
					A[i * 9 + j] += L[k * 9 + i] * L[k * 9 + j];
				}
				A[j * 9 + i] = A[i * 9 + j];
			}
		}
		// Solve Lx=0
		Jacobi( A, eig_value, eig_vec, 9, 30 );
		for( i = 1, flag = 0, min = fabs( eig_value[0] ); i < 9; i++ ) {
			if( fabs( eig_value[i] ) < fabs( min ) ) {
				flag = i;
				min = eig_value[i];
			}
		}
		for( i = 0; i < 9; i++ ) x[i] = eig_vec[i * 9 + flag];

		// Maximum Likelihood Estimation
//		hmg::M=M;		
//		double *param;
//		int buf;
//		getvec(x, param, buf);
//		levmar(hmg::f, &m[j*_n_pt*2], param, _n_pt*2, 9, 10, 0.001, ERR, ERR, ERR);
//		x=vec(param, 9);
//		freevec(param);

		// Save the solution
		for( i = 0; i < 9; i++ ) _H[a * 9 + i] = x[i];

		// testing
	//	mat H(_H[j], 3, 3);
	//	vec MM(3); MM[0]=M[(j*_n_pt)*2+0]; MM[1]=M[(j*_n_pt)*2+1]; MM[2]=1.0;
	//	cout<< m[(j*_n_pt)*2+0]<< ' '<< m[(j*_n_pt)*2+1]<< endl
	//		<< H[0]*MM/(H[2]*MM)<< ' '<< H[1]*MM/(H[2]*MM)<< endl;
	}

	return true;
}
#undef l

bool CAMERA_PARAM::SolveCalib( double *m, double *M ){
	// check is this object created for read parameter file
	if( _H == NULL ) return false;

	int i, j, k;
	double *V, *v12, *v11, *v22, *eig_value, *eig_vec, *b, *A,
		*h, min, lambda;

	// memory pool
	i = 12 * _n_img + 96;
	if( ( int )_m_size < i ){
		_m_size = ( unsigned )i;
		_mem = realloc( _mem, _m_size );
	}
	V = ( double* )_mem;		// 2_n_img * 6
	v12 = V + 12 * _n_img;		// 6
	v11 = v12 + 6;				// 6
	v22 = v11 + 6;				// 6
	eig_value = v22 + 6;		// 6
	eig_vec = eig_value + 6;	// 6 * 6
	A = eig_vec + 36;			// 6 * 6
	b = v12;

	for( i = 0; i < _n_img; i++ ){
		// the i-th image
		h = _H + i * 9;

		#define h(i, j) h[( (j) - 1 ) * 3 + ( (i) - 1 )]
	
		// calculate v12
		v12[0] = h(1, 1) * h(2, 1);
		v12[1] = h(1, 1) * h(2, 2) + h(1, 2) * h(2, 1);
		v12[2] = h(1, 2) * h(2, 2);
		v12[3] = h(1, 3) * h(2, 1) + h(1, 1) * h(2, 3);
		v12[4] = h(1, 3) * h(2, 2) + h(1, 2) * h(2, 3);
		v12[5] = h(1, 3) * h(2, 3);
		
		// calculate v11
		v11[0] = h(1, 1) * h(1, 1);
		v11[1] = h(1, 1) * h(1, 2) + h(1, 2) * h(1, 1);
		v11[2] = h(1, 2) * h(1, 2);
		v11[3] = h(1, 3) * h(1, 1) + h(1, 1) * h(1, 3);
		v11[4] = h(1, 3) * h(1, 2) + h(1, 2) * h(1, 3);
		v11[5] = h(1, 3) * h(1, 3);
		
		// calculate v12
		v22[0] = h(2, 1) * h(2, 1);
		v22[1] = h(2, 1) * h(2, 2) + h(2, 2) * h(2, 1);
		v22[2] = h(2, 2) * h(2, 2);
		v22[3] = h(2, 3) * h(2, 1) + h(2, 1) * h(2, 3);
		v22[4] = h(2, 3) * h(2, 2) + h(2, 2) * h(2, 3);
		v22[5] = h(2, 3) * h(2, 3);
		
		#undef h 

		// generate V
		for( j = 0; j < 6; j++ ){
			V[i * 2 * 6 + j] = v12[j];
			V[( i * 2 + 1 ) * 6 + j] = v11[j] - v22[j];
		}
	}

	// A = VT * V
	for( i = 0; i < 6; i++ ){
		for( j = i; j < 6; j++ ){
			A[i * 6 + j] = 0.0;
			for( k = 0; k < 2 * _n_img; k++ ) A[i * 6 + j] += V[k * 6 + i] * V[k * 6 + j];
			A[j * 6 + i] = A[i * 6 + j];
		}
	}
	// calculate b of B=A(-t)*A(-1)
	Jacobi( A, eig_value, eig_vec, 6 );
	for( i = 1, j = 0, min = eig_value[0]; i < 6; i++ ) {
		if( eig_value[i] < min ) {
			j = i;
			min = eig_value[i];
		}
	}
	if( b[0] > 0.0 ) for( i = 0; i < 6; i++ ) b[i] = eig_vec[i * 6 + j];
	else for( i = 0; i < 6; i++ ) b[i] = -eig_vec[i * 6 + j];

	// calculate intrinsic parameters
	#define B11 b[0]
	#define B12 b[1]
	#define B22 b[2]
	#define B13 b[3]
	#define B23 b[4]
	#define B33 b[5]
	_v0 = ( B12 * B13 - B11 * B23) / ( B11 * B22 - B12 * B12 );
	lambda = B33 - ( B13 * B13 + _v0 * ( B12 * B13 - B11 * B23 ) ) / B11;
	_alpha = sqrt( lambda / B11 );
	_beta = sqrt( lambda * B11 / ( B11 * B22 - B12 * B12 ) );
	_gamma = -B12 * _alpha * _alpha * _beta / lambda;
	_u0 = _gamma * _v0 / _beta - B13 * _alpha * _alpha / lambda;

	// if the result is not correct
	if( _v0 * 2 == _v0 || _v0 * 2 == _v0 || _alpha * 2 == _alpha || _beta * 2 == _beta || _gamma * 2 == _gamma ) return false;

	#undef B11
	#undef B12
	#undef B22
	#undef B13
	#undef B23
	#undef B33
/*
	// Maximum Likelihood Estimation
	clb::M=M;				// set function
	clb::point=_n_pt;		// set point per image
	double *param;

	param=new double [5+12*_n_img];				// initialize param
	param[0]=_alpha;							// the top 5 intrinsic parameters
	param[1]=_gamma;
	param[2]=_u0;
	param[3]=_beta;
	param[4]=_v0;
	for(i=0; i<_n_img; i++) {					// initialize extrinsic parameters
		double temp[12];						// for each image
		GetExtrinsic(temp, i);
		for(j=0; j<12; j++) param[5+i*12+j]=temp[j];
	}

	try{
		// levmar algorithm
		levmar(clb::f, m, param, _n_pt*_n_img*2, 5+_n_img*12, 100, 0.001, ERR, ERR, ERR);
	}catch( ... ){
		delete param;
		return false;
	}
	// save the solution
	_alpha = param[0];
	_gamma = param[1];
	_u0 = param[2];
	_beta = param[3];
	_v0 = param[4];
/**/
//	delete param;
	return true;
}

bool CAMERA_PARAM::SolveDistor(double *m, double *M){
	double u, v, r2, *D, *d, Rt[12], xy[3], k[2];
	int i, j, n, row = 2 * _n_pt * _n_img;

	// memory pool
	i = 8 * SQUARE( _n_pt ) * SQUARE( _n_img );
	if( ( int )_m_size < i ){
		_m_size = ( unsigned )i;
		_mem = realloc( _mem, _m_size );
	}
	D = ( double* )_mem;			// 2_n_pt_*n_img * 2
	d = D + row * 2;				// 2*_n_pt*_n_img

	n = _n_pt * _n_img;
	for( i = 0; i < n; i++ ){
		if( i == 0 || ( i - 1 ) / _n_pt != i / _n_pt ) GetExtrinsic( Rt, i / _n_pt );
		// Set x, y = Rt * point
		for( j = 0; j < 3; j++ ){
			xy[j] = M[i * 2 + 0] * Rt[j * 4 + 0] + M[i * 2 + 1] * Rt[j * 4 + 1] + Rt[j * 4 + 3];
		}
		xy[0] /= xy[2];
		xy[1] /= xy[2];

		// Set u, v, r^2
		u = _alpha * xy[0] + _u0;
		v = _beta * xy[1] + _v0;
		r2 = xy[0] * xy[0] + xy[1] * xy[1];

		// Set matrix D
		D[( i * 2 + 0 ) * 2 + 0] = ( ( u - _u0 ) * r2 );
		D[( i * 2 + 0 ) * 2 + 1] = ( ( u - _u0 ) * SQUARE( r2 ) );
		D[( i * 2 + 1 ) * 2 + 0] = ( ( v - _v0 ) * r2 );
		D[( i * 2 + 1 ) * 2 + 1] = ( ( v - _v0 ) * SQUARE( r2 ) );

		// Set vector d
		d[i * 2 + 0] = m[i * 2 + 0] - u;
		d[i * 2 + 1] = m[i * 2 + 1] - v;
	}
	// Solve k by Least-Square Method
//	k = Eigen_method( D, d );
	LS_method( D, d, k, row, 2 );

	// Get Solution
	_k1 = k[0];
	_k2 = k[1];
	return true;
}

void CAMERA_PARAM::GetIntrinsic( double *A, int w, int h){
	if( h == 0 ) h = _h;
	if( w == 0 ) w = _w;

	A[0] = _alpha * w / _w;
	A[1] = _gamma * w / _w;
	A[2] = _u0 * w / _w;
	A[4] = _beta * h / _h;
	A[5] = _v0 * h / _h;
	A[3] = A[6] = A[7] = 0.0;
	A[8] = 1.0;
}

void CAMERA_PARAM::GetIntrinsicInv( double *A, int w, int h ){
	double rh = ( h == 0 ) ? 1.0 : ( double )h / ( double )_h;
	double rw = ( w == 0 ) ? 1.0 : ( double )w / ( double )_w;
	double a = _alpha * rw,
		r = _gamma * rw,
		u0 = _u0 * rw,
		b = _beta * rh,
		v0 = _v0 * rh;
	A[0] = 1.0 / a;
	A[1] = -r / ( a * b );
	A[2] = ( r * v0 - b * u0 ) / ( a * b );
	A[4] = 1.0 / b;
	A[5] = -v0 / b;
	A[3] = A[6] = A[7] = 0.0;
	A[8] = 1.0;
}

void CAMERA_PARAM::GetExtrinsic( double *Rt, int img_idx ){
	double A[9], r1[3], r2[3], r3[3], t[3], h1[3], h2[3], h3[3];
	double lambda;
	int i, j;
	
	// memory pool
	if( _m_size < 9 ){
		_m_size = 9;
		_mem = realloc( _mem, _m_size );
	}

	// initialize
	GetIntrinsic( A );
	inverse( A, 3, _mem );
	// H=[h1 h2 h3]
	h1[0] = _H[img_idx * 9 + 0]; 
	h1[1] = _H[img_idx * 9 + 3]; 
	h1[2] = _H[img_idx * 9 + 6]; 
	h2[0] = _H[img_idx * 9 + 1];  
	h2[1] = _H[img_idx * 9 + 4]; 
	h2[2] = _H[img_idx * 9 + 7]; 
	h3[0] = _H[img_idx * 9 + 2]; 
	h3[1] = _H[img_idx * 9 + 5];  
	h3[2] = _H[img_idx * 9 + 8]; 
	
	lambda = 0.0;
	for( i = 0; i < 3; i++ ){
		r1[i] = r2[i] = t[i] = 0.0;
		for( j = 0; j < 3; j++ ){
			r1[i] += A[i * 3 + j] * h1[j];	// r1
			r2[i] += A[i * 3 + j] * h2[j];	// r2
			t[i] += A[i * 3 + j] * h3[j];	// t
		}
		lambda += SQUARE( r1[i] );
	}
	lambda = sqrt( lambda );
	
	// normalize
	for( i = 0; i < 3; i++ ){
		r1[i] /= lambda;
		r2[i] /= lambda;
		t[i] /= lambda;
	}

	// r3
	r3[0] = r1[1] * r2[2] - r1[2] * r2[1];
	r3[1] = r1[2] * r2[0] - r1[0] * r2[2];
	r3[2] = r1[0] * r2[1] - r1[1] * r2[0];
	
	// generate a 3x4 matrix
	Rt[0] = r1[0];	Rt[1] = r2[0];	Rt[2] = r3[0];	Rt[3] = t[0];
	Rt[4] = r1[1];	Rt[5] = r2[1];	Rt[6] = r3[1];	Rt[7] = t[1];
	Rt[8] = r1[2];	Rt[9] = r2[2];	Rt[10] = r3[2];	Rt[11] = t[2];
}

void CAMERA_PARAM::GetExtrinsic( double *Rt ){
	for( int i = 0; i < 12; i++ ) Rt[i] = _Rt[i];
}

void CAMERA_PARAM::GetPPM( double *P, int index, int h, int w ){
	double A[9], Rt[12];
	int i, j, k;
	GetIntrinsic( A, h, w );
	GetExtrinsic( Rt, index );
	// P = A * Rt
	for( i = 0; i < 3; i++ ){
		for( j = 0; j < 4; j++ ){
			P[i * 4 + j] = 0.0;
			for( k = 0; k < 3; k++ ) P[i * 4 + j] += A[i * 3 + k] * Rt[k * 4 + j];
		}
	}
}

void CAMERA_PARAM::GetPPM( double *P, int h, int w ){
	double A[9], Rt[12];
	int i, j, k;
	GetIntrinsic( A, h, w );
	GetExtrinsic( Rt );
	// P = A * Rt
	for( i = 0; i < 3; i++ ){
		for( j = 0; j < 4; j++ ){
			P[i * 4 + j] = 0.0;
			for( k = 0; k < 3; k++ ) P[i * 4 + j] += A[i * 3 + k] * Rt[k * 4 + j];
		}
	}
}

void CAMERA_PARAM::SetIntrinsic( double a, double b, double r, double u0, double v0 ){
	_alpha = a;
	_beta = b;
	_gamma = r;
	_u0 = u0;
	_v0 = v0;
}

void CAMERA_PARAM::SetExtrinsic( double* m ){ 
	for( int i = 0; i < 12; i++ ) _Rt[i] = m[i];
}

void CAMERA_PARAM::SetSize( int w, int h ){
	_h = h;
	_w = w;
}

void CAMERA_PARAM::SetNum( int img, int pt ){
	_H = ( double* )realloc( _H, img * 10 * sizeof( double ) );
	_n_img = img;
	_n_pt = pt;
}

bool CAMERA_PARAM::WriteParam( char *filename ){
	FILE *outfile;
#if _IS_VC6_
	outfile = fopen( filename, "wb" );
#else
	fopen_s( &outfile, filename, "wb" );
#endif
	if( outfile == NULL ) return false;

	fwrite( &_h, sizeof( int ), 1, outfile );
	fwrite( &_w, sizeof( int ), 1, outfile );
	fwrite( &_alpha, sizeof( double ), 1, outfile );
	fwrite( &_beta, sizeof( double ), 1, outfile );
	fwrite( &_u0, sizeof( double ), 1, outfile );
	fwrite( &_v0, sizeof( double ), 1, outfile );
	fwrite( &_gamma, sizeof( double ), 1, outfile );
	fwrite( _Rt, sizeof( double ), 12, outfile );
	fclose( outfile );
	return true;
}

bool CAMERA_PARAM::ReadParam( char *filename ){
	FILE *infile;
#if _IS_VC6_
	infile = fopen( filename, "wb" );
#else
	fopen_s( &infile, filename, "wb" );
#endif
	if( infile == NULL ) return false;

	fread( &_h, sizeof( int ), 1, infile );
	fread( &_w, sizeof( int ), 1, infile );
	fread( &_alpha, sizeof( double ), 1, infile );
	fread( &_beta, sizeof( double ), 1, infile );
	fread( &_u0, sizeof( double ), 1, infile );
	fread( &_v0, sizeof( double ), 1, infile );
	fread( &_gamma, sizeof( double ), 1, infile );
	fread( _Rt, sizeof( double ), 12, infile );
	fclose( infile );
	return true;
}

void CAMERA_PARAM::Clear( void ){
	if( _H != NULL ){
		free( _H );
		_H = NULL;
	}
	if( _m_size ){
		free( _mem );
		_m_size = 0;
	}
}