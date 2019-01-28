
#ifndef __MATHBOX_REGRESSION__
#define __MATHBOX_REGRESSION__

#include <math.h>
#include <scene.h>
#include "mathbox_matrix.hpp"


// multiple regression 3D
bool	Regression3D			( point3f *data, int member_size, float answer[3] );

// multiple regression 2D
bool	Regression2D			( point2f *data, int member_size, float answer[2] );

// multiple regression of ellipse 2D
bool	RegressionEllipse2D		( point2f *data, int member_size, ellipsef &ellipse );


/*	Do
________________________________________________________________*/

// multiple regression 3D

bool Regression3D( point3f *_s, int _mn, float _a[3] )
{
	static float m[9];
	static float c[3];
	static float temp;
	static int i;

	memset( m, 0, 9*sizeof(float) );
	memset( c, 0, 3*sizeof(float) );

	for( i = 0; i < _mn; i++, _s++ )
	{
		m[0] += _s[0].x*_s[0].x;
		m[4] += _s[0].y*_s[0].y;
		m[2] += _s[0].x;
		m[6] += _s[0].x;
		m[5] += _s[0].y;
		m[7] += _s[0].y;

		temp = _s[0].x*_s[0].y;
		m[1] += temp;
		m[3] += temp;

		c[0] += _s[0].x*_s[0].z;
		c[1] += _s[0].y*_s[0].z;
		c[2] += _s[0].z;
	}
	m[8] = (float)_mn;

	static float inv[9];
	if( !mInverse( m, inv, 3 ) )
		return false;

	_a[0] = inv[0]*c[0] + inv[1]*c[1] + inv[2]*c[2];
	_a[1] = inv[3]*c[0] + inv[4]*c[1] + inv[5]*c[2];
	_a[2] = inv[6]*c[0] + inv[7]*c[1] + inv[8]*c[2];

	return true;
}

// multiple regression 2D

bool Regression2D( point2f *_s, int _mn, float _a[2] )
{
	static float m[4];
	static float c[2];
	int i;

	memset( m, 0, 4*sizeof(float) );
	memset( c, 0, 2*sizeof(float) );

	for( i = 0; i < _mn; i++, _s++ )
	{
		m[0] += _s[0].x*_s[0].x;
		m[1] += _s[0].x;
		m[2] += _s[0].x;
	}
	m[3] = (float)_mn;

	static float inv[4];
	if( !mInverse( m, inv, 2 ) )
		return false;

	_a[0] = inv[0]*c[0] + inv[1]*c[1];
	_a[1] = inv[2]*c[0] + inv[3]*c[1];

	return true;
}

// multiple regression of ellipse 2D

bool RegressionEllipse2D( point2f *_s, int _mn, ellipsef &_e )
{
	/*
		先求出橢圓中心( xbar, ybar )，旋轉角theta，將input的點反轉theta角後平移到原點( 0, 0 )，再求橢圓的長軸( alfa )短軸( beta )
	*/

	static int i;
	static float temp;
	static float sumx, sumy, sumxx, sumyy, sumxy;
	static float n;
	static point2f *ptr;

	sumx = sumy = sumxx = sumyy = sumxy = 0.f;
	n = (float)_mn;
	ptr = _s;

	for( i = 0; i < n; i++, ptr++ )
	{
		sumx += ptr[0].x;
		sumy += ptr[0].y;
		sumxx += ptr[0].x*_s[0].x;
		sumyy += ptr[0].y*_s[0].y;
		sumxy += ptr[0].x*_s[0].y;
	}

	static float xbar, ybar, varx, vary, covarxy;
	static float sumdxdx, sumdydy, sumdxdy;
	static float dx, dy;

	xbar = sumx/n;
	ybar = sumy/n;
	varx = sumxx/n;
	vary = sumyy/n;
	covarxy = sumxy/n;

	sumdxdx = sumdydy = sumdxdy = 0.f;

	ptr = _s;
	for( i = 0; i < n; i++, ptr++ )
	{
		dx = ptr[0].x - xbar;
		dy = ptr[0].y - ybar;
		sumdxdx += dx*dx;
		sumdydy += dy*dy;
		sumdxdy += dx*dy;
	}

	static float theta, c, s;
	if( ( temp = sumdydy - sumdxdx ) == 0.f )
		theta = .5f*3.1415926535f;
	else
		theta = .5f*atanf( 2.f*sumdxdy/temp );
	c = cosf( theta );
	s = sinf( theta );

	static float sumXX, sumYY, varX, varY;
	static float X, Y;

	sumXX = sumYY = varX = varY = 0.f;
	ptr = _s;
	for( i = 0; i < n; i++, ptr++ )
	{
		dx = ptr[0].x - xbar;
		dy = ptr[0].y - ybar;
		X = c*dx - s*dy;
		Y = s*dx + c*dy;
		sumXX += X*X;
		sumYY += Y*Y;
	}
	varX = sumXX/n;
	varY = sumYY/n;

	static float a, b;
	a = sqrtf( varX );
	b = sqrtf( varY );

	_e.alfa = a;
	_e.beta = b;
	_e.center = point2f( xbar, ybar );
	_e.theta = theta;

	return true;
}


#endif