#ifndef _CURVE_H_
#define _CURVE_H_

#include <stdarg.h>

class Curve{
protected:
	double *_p;
	int _n;		// number of control points

	void *_mem;
	int _m_size;

public:
	Curve( void );
	~Curve( void );

	int SetDegree( int n );	// # of control points = n + 1
	int SetCtrPointNum( int n );
	void SetCtrPoint( double [] );
	int SetCtrPoint( int, double [3] );
	int SetCtrPoint( int, double, double, double );
	int GetDegree( void );
	int GetCtrPoint( int, double [3] );
	
	virtual void Clean( void );
	virtual void TangentVec( double, double [3] );
	virtual void NormalVec( double, double [3] );
	virtual void BiNormalVec( double, double [3] );
	virtual void Point( double, double [3] ) = 0;
	virtual int Derivative( Curve * ) = 0;

	friend class Surface;
//	template< class C1, class C2 > friend class RuledSurface;
};

//
// Bazier Curve
//

class Bazier: public Curve{
public:
	Bazier( void );

	int DegEvaluate( void );	// return the degree evaluated
	void DeCasteljau( double, double [3] );
	virtual void Point( double, double [3] );
	virtual int Derivative( Curve * );
};

//
// B-Spline Curve
//

class BSpline: public Curve{
protected:
	double *_t;		// knot values
	int _k;			// order = degree + 1

	double _N( int, int, double );
public:
	BSpline( void );

	int SetOrder( int );
	void SetKnots( double [] );
	int SetKnots( int, double );
	int SetKnots( int, ... );		// set all
	double GetKnots( int );
	double ParamStart( void );
	double ParamEnd( void );
	int GetOrder( void );

	virtual void Clean( void );
	virtual void Point( double, double[3] );
	virtual int Derivative( Curve * );
};

#endif