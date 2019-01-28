//#pragma comment( lib, "..\\matrixLib\\Release\\matrixLib.lib" )

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <typeinfo>
#include <exception>
#include "surface.h"
//#include "depthmap.h"
#include "numerical\linear.h"

// macro

#define SQR(x) ( (x) * (x) )
#define CUBE(x) ( (x) * (x) * (x) )
//#define min( a, b ) ( ( (a) < (b) ) ? (a) : (b) )
//#define max( a, b ) ( ( (a) > (b) ) ? (a) : (b) )
#define LEN2( a ) ( SQR( (a)[0] ) + SQR( (a)[1] ) + SQR( (a)[2] ) )
#define LENG( a ) sqrt( LEN2(a) )

#define ERR_BOUND			1.5
#define SMALL				0.000001
#ifndef M_PI
#define M_PI				3.1415926535897932384626433832795
#endif

// template functions

template< typename T > T powi( T &x, int a ){
	T ans = 1;
	for( int i = 0; i < a; i++ ) ans *= x;
	return ans;
}

// inline functions

inline double C( int a, int b ){
	double ans = 1;
	int i;
	if( b + b > a ) b = a - b;
	for( i = 0; i < b; ){
		ans *= a;
		ans /= ++i;
	}
	return ans;
}

inline void Cross( double x[], double y[], double z[] ){
	z[0] = x[1] * y[2] - x[2] * y[1];
	z[1] = x[2] * y[0] - x[0] * y[2];
	z[2] = x[0] * y[1] - x[1] * y[0];
}

inline double Dot( double x[], double y[] ){
	return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

//
// Surface
//

Surface::Surface( void ){
	_p = NULL;
	_mem = NULL;
	_n = _m = _m_size = 0;;
}

Surface::~Surface( void ){
	Clean();
}

int Surface::SetDegree( int m, int n ){
	return SetCtrPointNum( m + 1, n + 1 );
}

int Surface::SetCtrPointNum( int m, int n ){
	if( n >= 0 && m >= 0 && _m * _n < m * n ){
		double *p;
		int i, j = _n * _m * 3;

		p = ( double* )malloc( m * n * 3 * sizeof( double ) );	// allocate
		if( _p != NULL ){
			for( i = 0; i < j; i++ ) p[i] = _p[i];
			free( _p );
		}
		_p = p;
		_n = n;
		_m = m;
		return 1;
	}
	else return 0;
}

void Surface::SetCtrPoint( double p[] ){
	int i, n = 3 * _m * _n;

	for( i = 0; i < n; i++ ) _p[i] = p[i];
}

int Surface::SetCtrPoint( int m, int n, double p[3] ){
	if( n >= 0 && n < _n && m >= 0 && m < _m ){
		double *pp = _p + ( m * _n + n ) * 3;
		pp[0] = p[0];
		pp[1] = p[1];
		pp[2] = p[2];
		return 1;
	}
	else return 0;
}

int Surface::SetCtrPoint( int m, int n, ... ){
	if( n >= 0 && n < _n && m >= 0 && m < _m ){
		va_list va;
		double *p =  _p + ( m * _n + n ) * 3;

		va_start( va, n );
		p[0] = va_arg( va, double );
		p[1] = va_arg( va, double );
		p[2] = va_arg( va, double );
		va_end( va );
		return 1;
	}
	else return 0;
}

int Surface::GetCtrPoint( int m, int n, double p[3] ){
	if( n >= 0 && n < _n ){
		double *pp = _p + ( m * _n + n ) * 3;
		p[0] = pp[0];
		p[1] = pp[1];
		p[2] = pp[2];
		return 1;
	}
	else return 0;
}

void Surface::GetDegree( int *m, int *n ){
	*m = _m - 1;
	*n = _n - 1;
}

void Surface::Clean( void ){
	if( _p ) { 
		free( _p ); 
		_n = _m = 0;
		_p = NULL;
	}
	if( _mem ) { 
		free( _mem ); 
		_m_size = 0;
		_mem = NULL;
	}
}


double Surface::Nearest( double p0[3], double *pu, double *pv ){
	const double e = 1.0 / SMALL;
	double du, dv, p1[3], p2[3], p3[3], p4[3], p[3],
		unit, u, v, min = DBL_MAX, d, temp;
	int i, j, k, m, n;

	// initialize
	GetDegree( &m, &n );
	n = ( m > n ) ? m : n;
	unit = 1.0 / ( double )n;

	// initial guess
	u = 0.0;
	for( i = 0; i < n; i++, u += unit ){
		v = 0.0;
		for( j = 0; j < n; j++, v += unit ){
			// compute distance
			Point( u, v, p );
			d = 0.0;
			for( k = 0; k < 3; k++ ) {
				temp = p[k] - p0[k];
				d += SQR( temp );
			}
			// minimize
			if( d < min ) {
				min = d;
				*pu = u;
				*pv = v;
			}
		}
	}
	// initial guess
	u = *pu;
	v = *pv;
	// solve by gradient descent
	do{
		// compute gradient = 2 * [ Pu * ( P - P0 ) Pv * ( P - P0 ) ]
		Point( u, v, p );
		Point( u + SMALL, v, p1 );
		Point( u - SMALL, v, p2 );
		Point( u, v + SMALL, p3 );
		Point( u, v - SMALL, p4 );
		// gradient
		du = dv = 0.0;
		for( k = 0; k < 3; k++ ) {
			temp = ( p[k] - p0[k] ) * e;
			du += ( p1[k] - p2[k] ) * temp;
			dv += ( p3[k] - p4[k] ) * temp;
		}
		// update
		u -= 0.01 * du;
		v -= 0.01 * dv;
		if( u > 1.0 ) { du -= u - 1.0; u = 1.0; }
		else if( u < 0.0 ) { du -= u; u = 0.0; }
		if( v > 1.0 ) { dv -= v - 1.0; v = 1.0; }
		else if( v < 0.0 ) { dv -= v; v = 0.0; }
	} while( SQR( du ) + SQR( dv ) > SMALL );
	// compute distance
	Point( u, v, p );
	d = 0.0;
	for( k = 0; k < 3; k++ ) {
		temp = p[k] - p0[k];
		d += SQR( temp );
	}
	min = d;
	*pu = u;
	*pv = v;
	return sqrt( min );
}

//
// Bazier Surface
//

BazierSurface::BazierSurface( void ): Surface(){}

void BazierSurface::Point( double u, double v, double p[3] ){
	const double u2 = 1.0 - u, v2 = 1.0 - v;
	double *pp, bi, bj;
	int i, j;

	memset( p, 0x00, 3 * sizeof( double ) );
	for( i = 0; i <= _n; i++ ){
		bi = C( i, _n ) * powi( u, i ) * ( u2, _n - i );
		for( j = 0; j <= _m; j++ ){
			pp = _p + ( j * _n + i ) * 3;
			bj = C( j, _m ) * powi( v, j ) * ( v2, _m - j );
			bj *= bi;
			p[0] += pp[0] * bj;
			p[1] += pp[1] * bj;
			p[2] += pp[2] * bj;
		}
	}
}
