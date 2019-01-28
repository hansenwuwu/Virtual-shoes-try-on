
#ifndef __MATHBOX_VECTOR__
#define __MATHBOX_VECTOR__

#include <math.h>
#include <scene.h>


float	vLength				( const vector3f &vector );
float	vDot				( const vector3f &vector_a, const vector3f &vector_b );
void	vCross				( const vector3f &vector_a, const vector3f &vector_b, vector3f &result_vector );
void	vNormalize			( vector3f &vector );
float	vCos				( const vector3f &vector_fix, const vector3f &vector );


/*	Do
________________________________________________________________*/

float vLength( const vector3f &v )
{
	return sqrtf( vDot( v, v ) );
}

float vDot( const vector3f &va, const vector3f &vb )
{
	return va.x*vb.x + va.y*vb.y + va.z*vb.z;
}

void vCross( const vector3f &va, const vector3f &vb, vector3f &v )
{
	v = vector3f( va.y*vb.z - vb.y*va.z, -va.x*vb.z + vb.x*va.z, va.x*vb.y - vb.x*va.y );
}

void vNormalize( vector3f &v )
{
	v *= 1.f/vLength( v );
}

float vCos( const vector3f &va, const vector3f &vb )
{
	return vDot( va, vb )/( vLength( va )*vLength( vb ) );
}


#endif