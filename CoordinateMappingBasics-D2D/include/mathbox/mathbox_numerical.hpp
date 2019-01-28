
#ifndef __MATHBOX_NUMERICAL__
#define __MATHBOX_NUMERICAL__

#include "mathbox_matrix.hpp"


bool Linear3Element1Time( const float &f1, const float &f2, const float &f3, float *xyz );

/*	Do
________________________________________________________________*/

bool Linear3Element1Time( const float *f1, const float *f2, const float *f3, float *xyz )
{
	static float d;
	static float m[9];
	/*
	float m[9] = {
		f1[0], f1[1], f1[2],
		f2[0], f2[1], f2[2],
		f3[0], f3[1], f3[2]
	};
	*/
	memcpy( m, f1, 3*sizeof(float) );
	memcpy( m + 3, f2, 3*sizeof(float) );
	memcpy( m + 6, f3, 3*sizeof(float) );

	if( ( d = mDet( m, 3 ) ) == 0 )
		return false;

	static float c[9];
	d = 1.f/d;
	memcpy( c, m, 9*sizeof(float) );
	c[0] = -f1[3];
	c[3] = -f2[3];
	c[6] = -f3[3];
	xyz[0] = mDet( c, 3 )*d;

	memcpy( c, m, 9*sizeof(float) );
	c[1] = -f1[3];
	c[4] = -f2[3];
	c[7] = -f3[3];
	xyz[1] = mDet( c, 3 )*d;
	
	memcpy( c, m, 9*sizeof(float) );
	c[2] = -f1[3];
	c[5] = -f2[3];
	c[8] = -f3[3];
	xyz[2] = mDet( c, 3 )*d;
	return true;
}

#endif