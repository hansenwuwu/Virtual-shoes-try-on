
//#pragma warning( disable: 4793 )

#include <stdlib.h>
#include <math.h>
#include "numerical\linear.h"
#include <typeinfo>
#include "curve.h"

//
// macor and define
//

#define ERR		0.000001
#define SQR(x)	( (x) * (x) )

//
// inline functions
//

inline unsigned int Fact( int x, int n ){
	int i;
	unsigned int ans = 1;
	for( i = 0; i < n; i++ ) ans *= x--;
	return ans;
}

inline unsigned int C( int a, int b ){
	int i = 1, j;
	unsigned int ans = 1;
	
	if( b + b > a ) b = a - b;				// check b > a/2
	j = a - b;
	while( i <= b ){
		ans *= a--;
		if( ans % i == 0 ) ans /= i++;
	}
	for( i = a - b; a > j; a-- ) ans *= a; 
	return ans;
}

template < typename T > T powi( T &x, int n ){
	T ans = 1;
	int i;
	for( i = 0; i < n; i++ ) ans *= x;
	return ans;
}

//
// Curve
//

Curve::Curve( void ){
	_p = NULL;
	_mem = NULL;
	_n = _m_size = NULL;
}

Curve::~Curve( void ){
	Clean();
}

int Curve::SetDegree( int n ){
	return SetCtrPointNum( n + 1 );
}

int Curve::SetCtrPointNum( int n ){
	if( n >= 0 && _n < n ){
		double *p;
		int i, j = _n * 3;

		p = ( double* )malloc( n * 3 * sizeof( double ) );	// allocate
		if( _p != NULL ){
			for( i = 0; i < j; i++ ) p[i] = _p[i];			// copy
			free( _p );										// free
		}
		_p = p;
		_n = n;
		return 1;
	}
	else return 0;
}

void Curve::SetCtrPoint( double p[] ){
	int i, n = 3 * _n;

	for( i = 0; i < n; i++ ) _p[i] = p[i];
}

int Curve::SetCtrPoint( int n, double p[3] ){
	if( n >= 0 && n < _n ){
		n *= 3;
		_p[ n + 0 ] = p[0];
		_p[ n + 1 ] = p[1];
		_p[ n + 2 ] = p[2];
		return 1;
	}
	else return 0;
}

int Curve::SetCtrPoint( int i, double x, double y, double z ){
	if( i >= 0 && i < _n ){
		double *p = _p + i * 3;
		p[0] = x;
		p[1] = y;
		p[2] = z;
		return 1;
	}
	else return 0;
}

int Curve::GetDegree( void ){
	return _n - 1;
}

int Curve::GetCtrPoint( int n, double p[3] ){
	if( n >= 0 && n < _n ){
		n *= 3;
		p[0] = _p[ n + 0 ];
		p[1] = _p[ n + 1 ];
		p[2] = _p[ n + 2 ];
		return 1;
	}
	else return 0;
}

//
// virtual functions
//

void Curve::Clean( void ){
	if( _p ) { free( _p ); _n = 0; }
	if( _mem ) { free( _mem ); _m_size = 0; }
}

void Curve::TangentVec( double u, double p[3] ){
	double p1[3], p2[3], norm = 0.0;
	int i;
	
	Point( u + ERR, p1 );
	Point( u - ERR, p2 );
	for( i = 0; i < 3; i++ ) {
		p[i] = ( p1[i] - p2[i] ) * 0.5;
		norm += SQR( p[i] );
	}
	norm = 1.0 / sqrt( norm );
	for( i = 0; i < 3; i++ ) p[i] *= norm;
}

void Curve::NormalVec( double u, double p[3] ){
	double p0[3], p1[3], p2[3], pu[3], puu[3], 
		dot = 0.0, norm2 = 0.0;
	int i;
	
	// compute pu and puu
	Point( u, p0 );
	Point( u + ERR, p1 );
	Point( u - ERR, p2 );
	for( i = 0; i < 3; i++ ) {
		pu[i] = ( p1[i] - p2[i] ) * 0.5;
		puu[i] = p1[i] + p2[i] - p0[i] - p0[i];
		dot += pu[i] * puu[i];		// pu * puu
		norm2 += SQR( pu[i] );		// || pu ||2
	}

	// compue ( puu * pu ) / || pu ||2
	dot /= norm2;
	norm2 = 0.0;
	for( i = 0; i < 3; i++ ) {
		p[i] = puu[i] - dot * pu[i];
		norm2 += SQR( p[i] );
	}
	// normalize
	norm2 = 1.0 / sqrt( norm2 );
	for( i = 0; i < 3; i++ ) p[i] *= norm2;
}

void Curve::BiNormalVec( double u, double p[3] ){
	double p0[3], p1[3], p2[3], pu[3], puu[3], 
		dot = 0.0, norm2 = 0.0;
	int i;
	
	// compute pu and puu
	Point( u, p0 );
	Point( u + ERR, p1 );
	Point( u - ERR, p2 );
	for( i = 0; i < 3; i++ ) {
		pu[i] = ( p1[i] - p2[i] ) * 0.5;
		puu[i] = p1[i] + p2[i] - p0[i] - p0[i];
		dot += pu[i] * puu[i];		// pu * puu
		norm2 += SQR( pu[i] );		// || pu ||2
	}

	// compue ( puu * pu ) / || pu ||2
	dot /= norm2;
	norm2 = 0.0;
	for( i = 0; i < 3; i++ ) {
		p[i] = puu[i] - dot * pu[i];
		norm2 += SQR( p[i] );
	}
	// p1: principle normal vector form normalize
	norm2 = 1.0 / sqrt( norm2 );
	for( i = 0; i < 3; i++ ) p1[i] = p[i] * norm2;
	// p2: tangent vector from normalize
	norm2 = 0.0;
	for( i = 0; i < 3; i++ ) norm2 += SQR( pu[i] );
	norm2 = 1.0 / sqrt( norm2 );
	for( i = 0; i < 3; i++ ) p2[i] = pu[i] * norm2;

	// p = p1 x p2;
	p[0] = p2[1] * p1[2] - p2[2] * p1[1];
	p[1] = p2[2] * p1[0] - p2[0] * p1[2];
	p[2] = p2[0] * p1[1] - p2[1] * p1[0];
}

//
// Bazier Curve
//

Bazier::Bazier( void ): Curve(){}

int Bazier::DegEvaluate( void ){
	if( _n == 0 ) return 0;

	double *p, *p_new, *p1, *p2, u, v, unit;
	int i, j, n = _n;

	// initialize
	unit = 1.0 / ( double )_n;	// new degree ^-1
	u = unit;					// i / ( n+1 )
	v = 1.0 - u;				// 1 - v
	p = ( double* )malloc( ++_n * 3 * sizeof( double ) );
	j = n * 3;
	for( i = 0; i < j; i++ ) p[i] = _p[i];

	// evaluation
	p1 = p + n * 3;
	p2 = p1 - 3;
	p1[0] = p2[0];	p1[1] = p2[1];	p1[2] = p2[2];	// pn+1' =pn 
	p_new = p + 3;
	p2 = _p + 3;
	p1 = _p;
	for( i = 1; i < n; i++, u += unit, v -= unit ){
		p_new[0] = u * p1[0] + v * p2[0];
		p_new[1] = u * p1[1] + v * p2[1];
		p_new[2] = u * p1[2] + v * p2[2];
		p_new += 3;
		p1 = p2;
		p2 += 3;
	}
	free( _p );					// free
	_p = p;						// assign
	return n;
}

void Bazier::DeCasteljau( double u, double p[3] ){
	double v = 1.0 - u, *pj, *p1, *p2;
	int i, j, n = _n - 1;
	
	// memory pool
	j = _n * 3;
	i = j * sizeof( double );
	if( _m_size < i ){
		_m_size = i;
		_mem = realloc( _mem, _m_size );
	}
	pj = ( double* )_mem;

	for( i = 0; i < j; i++ ) pj[i] = _p[i];	// initialize
	
	for( i = 0; i < _n - 1; i++ ){
		for( j = 0; j < n; j++ ){
			p1 = pj + j * 3;
			p2 = p1 + 3;
			p1[0] = v * p1[0] + u * p2[0];
			p1[1] = v * p1[1] + u * p2[1];
			p1[2] = v * p1[2] + u * p2[2];
		}
		n--;
	}
	p[0] = pj[0];	p[1] = pj[1];	p[2] = pj[2];
}

//
// virtual functions
//

void Bazier::Point( double u, double p[3] ){
	double b, v = 1.0 - u;
	int i, j, n = _n - 1;

	p[0] = p[1] = p[2] = 0.0;			// initialize
	for( i = 0; i <= n; i++ ){
		j = i * 3;
		b = ( double )C( n, i );
		b *= powi( u, i ) * powi( v, n - i );
		p[0] += b * _p[j + 0];
		p[1] += b * _p[j + 1];
		p[2] += b * _p[j + 2];
	}
}

int Bazier::Derivative( Curve *bz ){
	// check curve type
	if( typeid( *bz ) != typeid( Bazier ) ) return 0;

	double *p, n = ( double )_n;
	int i, n1 = _n - 1;

	bz->SetCtrPointNum( _n - 1 );
	p = _p;
	for( i = 0; i < n1; i++, p += 3 ){
		bz->SetCtrPoint( i, 
			n * ( p[3] - p[0] ),
			n * ( p[4] - p[1] ),
			n * ( p[5] - p[2] ) );
	}
	return 1;
}

//
// B-Spline
//

BSpline::BSpline( void ){
	_t = NULL;
	_k = 0;
}

double BSpline::_N( int i, int k, double u ){
	double ans = 0.0, w;

	if( k == 1 ){
		if( _t[i] == _t[ i + 1 ] ) return 0.0;
		else if( _t[i] > u || _t[ i + 1 ] < u ) return 0.0;
		else return 1.0;
	}
	if( ( w = _N( i, k - 1, u ) ) != 0.0 )
		ans += ( u - _t[i] ) * w / ( _t[ i + k - 1 ] - _t[i] );
	if( ( w = _N( i + 1, k - 1, u ) ) != 0.0 )
		ans += ( _t[ i + k ] - u ) * w / ( _t[ i + k ] - _t[ i + 1 ] );
	return ans;
}

int BSpline::SetOrder( int k ){
	if( k >= 0 && _k < k ){
		double *t;
		int i, j = _n + k;

		t = ( double* )malloc( j * sizeof( double ) );	// allocate
		if( _t != NULL ){
			for( i = _n + _k; i >= 0; i-- ) t[i] = _t[i];// copy
			free( _t );									// free
		}
		_t = t;
		_k = k;

		return 1;
	}
	else return 0;
}

void BSpline::SetKnots( double t[] ){
	int i;
	for( i = _n + _k; i >= 0; i-- ) _t[i] = t[i];
}

int BSpline::SetKnots( int i, double t ){
	if( i >= 0 && i < _n + _k ) {
		_t[i] = t;
		return 1;
	}
	else return 0;
}
/**/
int BSpline::SetKnots( int n, ... ){
	// maximum memory
	if( n > _n + _k ) return 0;
	// set
	int i;
	va_list va;
	va_start( va, n );
	for( i = 0; i < n; i++ ) _t[i] = va_arg( va, double );
	va_end( va );

	return 1;
}
/**/
double BSpline::GetKnots( int i ){
	if( i >= 0 && i < _n + _k ) return _t[i];
	else return 0x00;
}

double BSpline::ParamStart( void ){
	return _t[ _n - 1 ];
}

double BSpline::ParamEnd( void ){
	return _t[_k];
}

void BSpline::Clean( void ){
	if( _n > 0 ) free( _p );
	if( _t != NULL ) free( _t );
	if( _m_size > 0 ) free( _mem );
	_n = _k = _m_size = 0;
}

void BSpline::Point( double u, double p[3] ){
	int i, n = _n - 1;
	double *pp, w;

	pp = _p;
	p[0] = p[1] = p[2] = 0.0;
	for( i = 0; i <= n ; i++, pp += 3 ){
		w = _N( i, _k, u );
		p[0] += w * pp[0];
		p[1] += w * pp[1];
		p[2] += w * pp[2];
	}
}

int BSpline::Derivative( Curve *c ){
	// check curve type
	if( typeid( *c ) != typeid( BSpline ) ) return 0;

	double *p, k1 = ( double )( _k - 1 ), temp;
	BSpline *bs = dynamic_cast< BSpline * >( c );
	int i, n1 = _n - 1;

	// Order - 1 and degree - 1
	bs->SetCtrPointNum( _n - 1 );
	bs->SetOrder( _k - 1 );
	p = _p;
	for( i = 1; i <= _n; i++, p += 3 ){
		temp = k1 / ( _t[ i + _k - 1 ] - _t[i] );
		bs->SetCtrPoint( i - 1,
			temp * p[3] - p[0],
			temp * p[4] - p[1],
			temp * p[5] - p[2] );
	}

	// Set knots value
	n1 = _n + _k - 1;
	for( i = 0; i < n1; i++ ) bs->SetKnots( i, _t[ i + 1 ] );

	return 1;
}