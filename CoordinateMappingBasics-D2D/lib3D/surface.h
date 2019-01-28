#ifndef _SURFACE_H_
#define _SURFACE_H_

#include "curve.h"
#include <stdarg.h>
#include <exception>

class DepthImage;

//
// structures
//

typedef struct _curvatures{
	double gauss, mean;
} Curvature2;

//
// Macro
//

#define SURF_TYPE			int

// Surface Types: 0bKHKHKH
#define SF_SADDLE_RIDGE		0x30		// 32 + 16 = 0b110000
#define SF_MINIMAL			0x24		// 32 + 4 = 0b100100
#define SF_SADDLE_VALLEY	0x21		// 32 + 1 = 0b100001
#define SF_RIDGE			0x18		// 16 + 8 = 0b011000
#define SF_FLAT				0x0c		// 8 + 4 = 0b001100
#define SF_VALLEY			0x09		// 8 + 1 = 0b001001
#define SF_PEAK				0x12		// 16 + 2 = 0b010010
#define SF_PIT				0x03		// 16 + 1 = 0b000011
#define SF_NONE				0x00

//
// Basic Surface
//

class Surface{
protected:
	double *_p;
	int _m, _n;		// n for u, m for v

	void *_mem;
	int _m_size;

public:
	Surface( void );
	~Surface( void );

	int SetDegree( int m, int n );
	int SetCtrPointNum( int m, int n );
	void SetCtrPoint( double [] );
	int SetCtrPoint( int m, int n, double [3] );
	int SetCtrPoint( int m, int n, ... );
	int GetCtrPoint( int m, int n, double [3] );

	virtual void GetDegree( int *m, int *n );
	virtual void Clean( void );
	virtual void Point( double, double, double [3] ) = 0;
	virtual double Nearest( double p[3], double *u, double *v );
};

class BazierSurface: public Surface{
public:
	BazierSurface( void );
	virtual void Point( double, double, double [3] );
};

template< class CV1, class CV2 >
class RuledSurface: public Surface{
protected:
	CV1 _a;
	CV2 _b;

public:
	RuledSurface( void ): Surface(){}

	Curve& GetCurve( int i );
	virtual void GetDegree( int *m, int *n );
	virtual void Clean( void );
	virtual void Point( double, double, double [3] );
	virtual double Nearest( double p[3], double *u, double *v );
};

//
// template class functions
//

template < class CV1, class CV2 >
Curve& RuledSurface< CV1, CV2 >::GetCurve( int i ){
	switch( i ){
	case 1: return _a;
	case 2: return _b;
	default: throw 0;
	}
}

template < class CV1, class CV2 >
void RuledSurface< CV1, CV2 >::GetDegree( int *m, int *n ){
	*m = _a.GetDegree();
	*n = _b.GetDegree();
}

template < class CV1, class CV2 >
void RuledSurface< CV1, CV2 >::Clean( void ){
	_a.Clean();
	_b.Clean();
}

template < class CV1, class CV2 >
void RuledSurface< CV1, CV2 >::Point( double u, double v, double p[3] ){
	double p1[3], p2[3];
	int i;

	_a.Point( u, p1 );
	_b.Point( u, p2 );
	u = 1.0 - v;
	for( i = 0; i < 3; i++ ) p[i] = u * p1[i] + v * p2[i];
}

template < class CV1, class CV2 >
double  RuledSurface< CV1, CV2 >::Nearest( double p[], double *pu, double *pv ){
	return this->Surface::Nearest( p, pu, pv );
}

#endif