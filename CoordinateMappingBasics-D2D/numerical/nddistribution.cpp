
#include "distribution.h"
#include "linear.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <list>
using namespace distribution;

//
// defination and macro

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif
#define SQR( x ) ( (x) * (x) )

//
// inline functions
//

inline void MultiMatrix( double *a, double *b, double *c, int n ){
	int i, j, k;
	for( i = 0; i < n; i++ ){
		for( j = 0; j < n; j++ ){
			c[i * n + j] = 0.0;
			for( k = 0; k < n; k++ ) c[i * n + j] += a[i * n + k] * b[k * n + j];
		}
	}
}

inline void MultiVector( double *a, double *b, double *c, int n ){
	int i, j;

	for( i = 0; i < n; i++ ){
		c[i] = 0.0;
		for( j = 0; j < n; j++ ) c[i] += a[i * n + j] * b[j];
	}
}

inline void Transpose( double *a, int n ){
	int i, j;
	double temp;

	for( i = 0; i < n; i++ ){
		for( j = i + 1; j < n; j++ ){
			temp = a[i * n + j];
			a[i * n + j] = a[j * n + i];
			a[j * n + i] = temp;
		}
	}
}

//
// N-D Distribution
//

Distribution_ND::Distribution_ND( void ){
	_d = _m_size = 0;
	_a = _b = NULL;
	_mem = NULL;
}

Distribution_ND::Distribution_ND( int d ){
	_d = d;
	_a = ( double* )malloc( d * sizeof( double ) );
	_b = ( double* )malloc( d * d * sizeof( double ) );
	_m_size = 0;
	_mem = NULL;
}

Distribution_ND::Distribution_ND( int d, double *a, double *b ){
	int i;

	_d = d;
	_a = ( double* )malloc( d * sizeof( double ) );
	for( i = 0; i < d; i++ ) _a[i] = a[i];
	d *= d;
	_b = ( double* )malloc( d * sizeof( double ) );
	for( i = 0; i < d; i++ ) _b[i] = b[i];
	_m_size = 0;
	_mem = NULL;
}

Distribution_ND::~Distribution_ND( void ){
	if( _a ) free( _a );
	if( _b ) free( _b );
	if( _m_size ) free( _mem );
}

int Distribution_ND::GetDimension( void ){
	return _d;
}

//
// N-D Normal Distribution
//

Normal_ND::Normal_ND( int d ): Distribution_ND( d ){
	int i;
	for( i = 0; i < d; i++ ) {
		_a[i] = 0.0;
		_b[i * d + i] = 1.0;
	}
}

Normal_ND::Normal_ND( int d, double *a, double *b ): Distribution_ND( d, a, b ){}

void Normal_ND::SetMean( double *a ){
	int i;
	for( i = 0; i < _d; i++ ) _a[i] = a[i];
}

void Normal_ND::SetVar( double *b ){
	int i, d = _d * _d;
	for( i = 0; i < d; i++ ) _b[i] = b[i];
}

void Normal_ND::Mean( double *a ){
	int i;
	for( i = 0; i < _d; i++ ) a[i] = _a[i];
}

void Normal_ND::Var( double *b ){
	int i, d = _d * _d;
	for( i = 0; i < d; i++ ) b[i] = _b[i];
}

double Normal_ND::f( double *x ){
	double ans = 0.0, deter, *sig_inv, *u;
	int i, d;

	// memory pool
	i = ( 2 * _d * _d + _d ) * sizeof( double );
	if( i > _m_size ){
		_m_size = i;
		_mem = realloc( _mem, _m_size );
	}
	u = ( double* )_mem + _d * _d;
	sig_inv = u + _d;

	// initialize
	for( i = 0; i < _d; i++ ) x[i] -= _a[i];		// x - u
	d = _d * _d;
	for( i = 0; i < d; i++ ) sig_inv[i] = _b[i];	// var
	deter = det( sig_inv, _d, _mem );				// det( var )
	inverse( sig_inv, _d, _mem );					// var^-1

	MultiVector( sig_inv, x, u, _d );				// var^-1 * ( x - u )
	for( i = 0; i < _d; i++ ) ans += x[i] * u[i];	// ( x - u )T * var^-1 * ( x - u )
	ans = exp( -ans / 2.0 );
	deter = sqrt( 2.0 * M_PI * deter );
	ans /= deter;

	return ans;
}

void Normal_ND::Rand( double *X ){
	double u1, u2, *Aw, *var, *eigen;
	int i, j;

	// 1D case
	if( _d == 1 ){
		u1 = (double)rand() / RAND_MAX;
		u2 = (double)rand() / RAND_MAX;
		X[0] = sqrt( -2.0 * log(u1) ) * sin( 2.0 * M_PI * u2 );
		X[0] = X[0] * sqrt( _b[0] ) + _a[0];
		return;
	}

	// memory pool
	i = ( 2 * _d * _d + _d ) * sizeof( double );
	if( _m_size < i ) {
		_m_size = i;
		_mem = realloc( _mem, _m_size );
	}

	// initialize
	Aw = ( double* )_mem;
	var = Aw + _d * _d;
	eigen = var + _d * _d;

	// Box Muller Method
	for( i = 0; i < _d; i += 2 ){
		u1 = (double)rand() / RAND_MAX;
		u2 = (double)rand() / RAND_MAX;
		X[i] = sqrt( -2.0 * log(u1) ) * sin( 2.0 * M_PI * u2 );
		if( i + 1 < _d )
			X[i + 1] = sqrt( -2.0 * log(u1) ) * cos( 2.0 * M_PI * u2 );
	}

	// transform
	for( i = 0; i < _d * _d; i++ ) var[i] = _b[i];
	Jacobi( var, eigen, Aw, _d, 30 );
	for( i = 0; i < _d; i++ ){
		if( eigen[i] > 0 ){
			u1 = sqrt( eigen[i] );
			for( j = 0; j < _d; j++ ) Aw[j * _d + i] *= u1;
		}
		else{
			u1 = sqrt( -eigen[i] );
			for( j = 0; j < _d; j++ ) Aw[j * _d + i] *= -u1;
		}
	}
	MultiVector( Aw, X, eigen, _d );
	for( i = 0; i < _d; i++ ) X[i] = _a[i] + eigen[i];
}

//
// Mixture Model
//

Mixture::Mixture( int n, Distribution_ND *x ): Distribution_ND(){
	int i;
	double p = 1.0 / ( double )n;

	_d = n;													// number of mixture
	_a = ( double * )malloc( n * sizeof( double ) );		// probability
	_b = ( double * )malloc( n * n * sizeof( double ) );	// temp
	_X = ( Distribution_ND ** )calloc( n, sizeof( Distribution_ND * ) );

	// initialize
	for( i = 0; i < _d; i++ ) {
		_a[i] = p;
		_X[i] = x + i;
	}
}

Mixture::~Mixture( void ){
	if( _X ) free( _X );
}

void Mixture::SetProb( double *p ){
	for( int i = 0; i < _d; i++ ) _a[i] = p[i];
}

double Mixture::Prob( int i ){
	return _a[i];
}

void Mixture::SetMean( int i, double *x ){
	_X[i]->SetMean( x );
}

void Mixture::SetVar( int i, double *x ){
	_X[i]->SetVar( x );
}

void Mixture::Mean( int i, double *x ){
	_X[i]->Mean( x );
}

void Mixture::Var( int i, double *x ){
	_X[i]->Var( x );
}

double Mixture::f( int i, double *x ){
	return _X[i]->f( x );
}

void Mixture::SetMean( double *u ){
	int i, d = _X[0]->GetDimension();
	for( i = 0; i < _d; i++ ) _X[i]->SetMean( u + d * i );
}

void Mixture::SetVar( double *v ){
	int i, d = _X[0]->GetDimension();
	d = SQR( d );	// d * d
	for( i = 0; i < _d; i++ ) _X[i]->SetVar( v + d * i );
}

void Mixture::Mean( double *u ){
	int i, j, d = _X[0]->GetDimension();
	d *= d;

	// initialize
	for( i = 0; i < d; i++ ) _b[i] = 0.0;

	// mean of mixed model
	for( i = 0; i < _d; i++ ){
		_X[i]->Var( _b );
		for( j = 0; j < d; j++ ) u[j] += _a[i] * _a[i] * _b[j];
	}
}

void Mixture::Var( double *v ){
	/*
	do something
	*/
}

double Mixture::f( double *x ){
	double ans = 0.0;
	int i;

	// pi * fi(x)
	for( i = 0; i < _d; i++ ) ans += _a[i] * _X[i]->f( x );
	return ans;
}

void Mixture::Rand( double *x ){
	double r = ( double )rand() / ( double )RAND_MAX;
	int i;
	
	// select model
	for( i = 0; i < _d; i++ ){
		if( r <= _a[i] ){
			_X[i]->Rand( x );	// random
			return;
		}
		r -= _a[i];
	}
}
