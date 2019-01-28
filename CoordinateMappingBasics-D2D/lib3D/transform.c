
#include "track.h"
#include <math.h>
#include <string.h>
#include <float.h>
#include "numerical\linear2.h"

#include <windows.h>
#include <stdio.h>

//
// thread global variables
//

__declspec( thread )  extern double gxd, gyd, gzd, gcad, gcbd, gccd, gsad, gsbd, gscd;
__declspec( thread )  extern float gxf, gyf, gzf, gcaf, gcbf, gccf, gsaf, gsbf, gscf;

//
// macro
//

#define SQR(x)			( (x) * (x) )
#define SMALL			0.000000000000001
#define SMALL_ANGLE		0.2
#define M_PI			3.1415926535897932384626433832795
#define DegToRec(x)		( 0.01745329251994329576923690768489 * (x) )
#define RecToDeg(x)		( 57.295779513082320876798154814105 * (x) )
#define MAX2( a, b )	( ( (a) > (b) ) ? (a) : (b) )
#define MIN2( a, b )	( ( (a) < (b) ) ? (a) : (b) )
#define MAX3( a, b, c )	( ( (a) > (b) ) ? MAX2( (a), (c) ) : MAX2( (b), (c) ) )
#define MIN3( a, b, c )	( ( (a) < (b) ) ? MIN2( (a), (c) ) : MIN2( (b), (c) ) )
#define CompDist2f( d2, m, p1, p2 ) \
	gxf = (m)[0][0] * (p1).x + (m)[0][1] * (p1).y + (m)[0][2] * (p1).z + (m)[0][3] - (p2).x;\
	gyf = (m)[1][0] * (p1).x + (m)[1][1] * (p1).y + (m)[1][2] * (p1).z + (m)[1][3] - (p2).y;\
	gzf = (m)[2][0] * (p1).x + (m)[2][1] * (p1).y + (m)[2][2] * (p1).z + (m)[2][3] - (p2).z;\
	(d2) = SQR( gxf ) + SQR( gyf ) + SQR( gzf );

#define CompDist2d( d2, m, p1, p2 ) \
	gxd = (m)[0][0] * (p1).x + (m)[0][1] * (p1).y + (m)[0][2] * (p1).z + (m)[0][3] - (p2).x;\
	gyd = (m)[1][0] * (p1).x + (m)[1][1] * (p1).y + (m)[1][2] * (p1).z + (m)[1][3] - (p2).y;\
	gzd = (m)[2][0] * (p1).x + (m)[2][1] * (p1).y + (m)[2][2] * (p1).z + (m)[2][3] - (p2).z;\
	(d2) = SQR( gxd ) + SQR( gyd ) + SQR( gzd );

#define CompDist2fv( d2, m, p1, p2 ) \
	gxf = (m)[0][0] * (p1).ptr[0] + (m)[0][1] * (p1).ptr[1] + (m)[0][2] * (p1).ptr[2] + (m)[0][3] - (p2).ptr[0];\
	gyf = (m)[1][0] * (p1).ptr[0] + (m)[1][1] * (p1).ptr[1] + (m)[1][2] * (p1).ptr[2] + (m)[1][3] - (p2).ptr[1];\
	gzf = (m)[2][0] * (p1).ptr[0] + (m)[2][1] * (p1).ptr[1] + (m)[2][2] * (p1).ptr[2] + (m)[2][3] - (p2).ptr[2];\
	(d2) = SQR( gxf ) + SQR( gyf ) + SQR( gzf );

#define CompDist2dv( d2, m, p1, p2 ) \
	gxd = (m)[0][0] * (p1).ptr[0] + (m)[0][1] * (p1).ptr[1] + (m)[0][2] * (p1).ptr[2] + (m)[0][3] - (p2).ptr[0];\
	gyd = (m)[1][0] * (p1).ptr[0] + (m)[1][1] * (p1).ptr[1] + (m)[1][2] * (p1).ptr[2] + (m)[1][3] - (p2).ptr[1];\
	gzd = (m)[2][0] * (p1).ptr[0] + (m)[2][1] * (p1).ptr[1] + (m)[2][2] * (p1).ptr[2] + (m)[2][3] - (p2).ptr[2];\
	(d2) = SQR( gxd ) + SQR( gyd ) + SQR( gzd );

#define FastCompDist2f( d2, rt, p1, p2 ) \
	gxf = (p1).x - (rt).c * (p1).y + (rt).b * (p1).z + (rt).x - (p2).x;\
	gyf = (rt).c * (p1).x + (p1).y - (rt).a * (p1).z + (rt).y - (p2).y;\
	gzf = -(rt).b * (p1).x + (rt).a * (p1).y + (p1).z + (rt).z - (p2).z;\
	(d2) = SQR( gxf ) + SQR( gyf ) + SQR( gzf );

#define FastCompDist2d( d2, rt, p1, p2 ) \
	gxd = (p1).x - (rt).c * (p1).y + (rt).b * (p1).z + (rt).x - (p2).x;\
	gyd = (rt).c * (p1).x + (p1).y - (rt).a * (p1).z + (rt).y - (p2).y;\
	gzd = -(rt).b * (p1).x + (rt).a * (p1).y + (p1).z + (rt).z - (p2).z;\
	(d2) = SQR( gxd ) + SQR( gyd ) + SQR( gzd );

//
// Data Translation by Paired Points
//

float RigidTransformMatf( Point3f_Pair pa, float m[3][4], const Transform_Paramf param ){
	RT_Paramf theta;
	float ret;
	CompRTFromMatf( &theta, m );
	ret = RigidTransformf( pa, &theta, param );
	CompRTMatf( m, theta );
	return ret;
}

float RigidTransformf( Point3f_Pair pa, RT_Paramf *theta, const Transform_Paramf param ){
	const float ni = 1.0f / ( float )pa.n;
	RT_Paramf buf;
	float M[3][4], Ma[3][3], Mb[3][3], Mc[3][3];
	float e[3], grad[6], temp[3], a1, a2;
	struct _point_3d_float *p1, *p2;
	int i, k;

	for( k = 0; k < param.iter; k++ ){
		// compute cos ans sin
		gsaf = ( float )sin( theta->a );
		gsbf = ( float )sin( theta->b );
		gscf = ( float )sin( theta->c );
		gcaf = ( float )cos( theta->a );
		gcbf = ( float )cos( theta->b );
		gccf = ( float )cos( theta->c );
		
		// update matrix M
		M[0][0] = gccf * gcbf;	M[0][1] = gccf * gsbf * gsaf - gscf * gcaf;	M[0][2] = gccf * gsbf * gcaf + gscf * gsaf;	M[0][3] = theta->x;
		M[1][0] = gscf * gcbf;	M[1][1] = gscf * gsbf * gsaf + gccf * gcaf;	M[1][2] = gscf * gsbf * gcaf - gccf * gsaf;	M[1][3] = theta->y;
		M[2][0] = -gsbf;	M[2][1] = gcbf * gsaf;	M[2][2] = gcbf * gcaf;	M[2][3] = theta->z;

		// update matrix dM/da
//		Ma[0][0] = M[1][0] = M[2][0] = 0.0;
		Ma[0][1] = gccf * gsbf * gcaf + gscf * gsaf;	Ma[0][2] = -gccf * gsbf * gsaf + gscf * gcaf;
		Ma[1][1] = gscf * gsbf * gcaf - gccf * gsaf;	Ma[1][2] = -gscf * gsbf * gsaf - gccf * gcaf;
		Ma[2][1] = gcbf * gcaf;	Ma[2][2] = -gcbf * gsaf;

		// update matrix dM/db
		Mb[0][0] = -gccf * gsbf;	Mb[0][1] = gccf * gcbf * gsaf;	Mb[0][2] = gccf * gcbf * gcaf;
		Mb[1][0] = -gscf * gsbf;	Mb[1][1] = gscf * gcbf * gsaf;	Mb[1][2] = gscf * gcbf * gcaf;
		Mb[2][0] = -gcbf;	Mb[2][1] = -gsbf * gsaf;	Mb[2][2] = -gsbf * gcaf;

		// update matrix dM/dc
		Mc[0][0] = -gscf * gcbf;	Mc[0][1] = -gscf * gsbf * gsaf - gccf * gcaf;	Mc[0][2] = -gscf * gsbf * gcaf + gccf * gsaf;
		Mc[1][0] = gccf * gcbf;	Mc[1][1] = gccf * gsbf * gsaf - gscf * gcaf;	Mc[1][2] = gccf * gsbf * gcaf + gscf * gsaf;
//		Mc[2][0] = Mc[2][1] = Mc[2][2] = 0.0;

		e[0] = e[1] = e[2] = 0.0;		// zerolize error

		// compute gradient
		memset( grad, 0x00, 6 * sizeof( float ) );
		p1 = pa.pt1;
		p2 = pa.pt2;
		for( i = 0; i < pa.n; i++, p1++, p2++ ){
			// temp = ( M * p1 - p2 )T
			temp[0] = ( float )( M[0][0] * p1->x + M[0][1] * p1->y + M[0][2] * p1->z + M[0][3] - p2->x );
			temp[1] = ( float )( M[1][0] * p1->x + M[1][1] * p1->y + M[1][2] * p1->z + M[1][3] - p2->y );
			temp[2] = ( float )( M[2][0] * p1->x + M[2][1] * p1->y + M[2][2] * p1->z + M[2][3] - p2->z );
			e[0] += SQR( temp[0] ) + SQR( temp[1] ) + SQR( temp[2] );

			// gradient = tempT * [ e1 e2 e3 | Ma * p | Mb * p | Mc * p ]
			grad[0] += temp[0];
			grad[1] += temp[1];
			grad[2] += temp[2];
			grad[3] += ( float )(
				temp[0] * ( Ma[0][1] * p1->y + Ma[0][2] * p1->z ) +
				temp[1] * ( Ma[1][1] * p1->y + Ma[1][2] * p1->z ) +
				temp[2] * ( Ma[2][1] * p1->y + Ma[2][2] * p1->z ) );
			grad[4] += ( float )(
				temp[0] * ( Mb[0][0] * p1->x + Mb[0][1] * p1->y + Mb[0][2] * p1->z ) +
				temp[1] * ( Mb[1][0] * p1->x + Mb[1][1] * p1->y + Mb[1][2] * p1->z ) +
				temp[2] * ( Mb[2][0] * p1->x + Mb[2][1] * p1->y + Mb[2][2] * p1->z ) );
			grad[5] += ( float )(
				temp[0] * ( Mc[0][0] * p1->x + Mc[0][1] * p1->y + Mc[0][2] * p1->z ) +
				temp[1] * ( Mc[1][0] * p1->x + Mc[1][1] * p1->y + Mc[1][2] * p1->z ) );
		}
		e[0] *= ni;					// normalize error
		if( e[0] < param.err ) break;		// break criteria

		gxf = MAX3( grad[3], grad[4], grad[5] );
		gyf = MIN3( grad[3], grad[4], grad[5] );
		a1 = ( float )MAX2( fabs( gxf ), fabs( gyf ) );
		a1 = ( a1 < M_PI ) ? param.lr1 * ni : param.lr1 * ( float )M_PI / a1;
		a2 = param.lr2 * a1;

		// compute error at a1
		buf.x = theta->x - a1 * grad[0];	buf.y = theta->y - a1 * grad[1];	buf.z = theta->z - a1 * grad[2]; 
		buf.a = theta->a - a1 * grad[3];	buf.b = theta->b - a1 * grad[4];	buf.c = theta->c - a1 * grad[5]; 
		CompRTMatf( M, buf );
		p1 = pa.pt1;
		p2 = pa.pt2;
		for( i = 0; i < pa.n; i++, p1++, p2++ ) {
			CompDist2f( temp[0], M, *p1, *p2 );
			e[1] += temp[0];
		}
		e[1] *= ni;

		// compute error at a2 = 0.5
		buf.x = theta->x - a2 * grad[0];	buf.y = theta->y - a2 * grad[1];	buf.z = theta->z - a2 * grad[2]; 
		buf.a = theta->a - a2 * grad[3];	buf.b = theta->b - a2 * grad[4];	buf.c = theta->c - a2 * grad[5]; 
		CompRTMatf( M, buf );
		p1 = pa.pt1;
		p2 = pa.pt2;
		for( i = 0; i < pa.n; i++, p1++, p2++ ) {
			CompDist2f( temp[0], M, *p1, *p2 );
			e[2] += temp[0];
		}
		e[2] *= ni;

		// find learning rate that minimize object function
		gxf = a1 * ( e[2] - e[0] ) - a2 * ( e[1] - e[0] );
		if( gxf == 0.0 ) gzf = a2;		// if solution is max, then set to solution, else set to min step
		else{
			gyf = SQR( a2 ) * ( e[1] - e[0] ) - SQR( a1 ) * ( e[2] - e[0] );
			gzf = -0.5f * gyf / gxf;
		}
		if( gzf < SMALL ) {
			if( e[0] < e[1] && e[0] < e[2] ) break;
			if( e[1] < e[2] ) gzf = a1;
			else gzf = a2;
		}
		// update agccording gradient
		theta->x -= gzf * grad[0];
		theta->y -= gzf * grad[1];
		theta->z -= gzf * grad[2];
		theta->a -= gzf * grad[3];
		theta->b -= gzf * grad[4];
		theta->c -= gzf * grad[5];
	}
	return e[0];
}

float RigidTransformMatfv( Point3fv_Pair pa, float m[3][4], const Transform_Paramf param ){
	RT_Paramf theta;
	float ret;
	CompRTFromMatf( &theta, m );
	ret = RigidTransformfv( pa, &theta, param );
	CompRTMatf( m, theta );
	return ret;
}

float RigidTransformfv( Point3fv_Pair pa, RT_Paramf *theta, const Transform_Paramf param ){
	const float ni = 1.0f / ( float )pa.n;
	RT_Paramf buf;
	float M[3][4], Ma[3][3], Mb[3][3], Mc[3][3];
	float e[3], grad[6], temp[3], a1, a2;
	struct _point_3d_float_ptr *p1, *p2;
	int i, k;

	for( k = 0; k < param.iter; k++ ){
		// compute cos ans sin
		gsaf = ( float )sin( theta->a );
		gsbf = ( float )sin( theta->b );
		gscf = ( float )sin( theta->c );
		gcaf = ( float )cos( theta->a );
		gcbf = ( float )cos( theta->b );
		gccf = ( float )cos( theta->c );
		
		// update matrix M
		M[0][0] = gccf * gcbf;	M[0][1] = gccf * gsbf * gsaf - gscf * gcaf;	M[0][2] = gccf * gsbf * gcaf + gscf * gsaf;	M[0][3] = theta->x;
		M[1][0] = gscf * gcbf;	M[1][1] = gscf * gsbf * gsaf + gccf * gcaf;	M[1][2] = gscf * gsbf * gcaf - gccf * gsaf;	M[1][3] = theta->y;
		M[2][0] = -gsbf;	M[2][1] = gcbf * gsaf;	M[2][2] = gcbf * gcaf;	M[2][3] = theta->z;

		// update matrix dM/da
//		Ma[0][0] = M[1][0] = M[2][0] = 0.0;
		Ma[0][1] = gccf * gsbf * gcaf + gscf * gsaf;	Ma[0][2] = -gccf * gsbf * gsaf + gscf * gcaf;
		Ma[1][1] = gscf * gsbf * gcaf - gccf * gsaf;	Ma[1][2] = -gscf * gsbf * gsaf - gccf * gcaf;
		Ma[2][1] = gcbf * gcaf;	Ma[2][2] = -gcbf * gsaf;

		// update matrix dM/db
		Mb[0][0] = -gccf * gsbf;	Mb[0][1] = gccf * gcbf * gsaf;	Mb[0][2] = gccf * gcbf * gcaf;
		Mb[1][0] = -gscf * gsbf;	Mb[1][1] = gscf * gcbf * gsaf;	Mb[1][2] = gscf * gcbf * gcaf;
		Mb[2][0] = -gcbf;	Mb[2][1] = -gsbf * gsaf;	Mb[2][2] = -gsbf * gcaf;

		// update matrix dM/dc
		Mc[0][0] = -gscf * gcbf;	Mc[0][1] = -gscf * gsbf * gsaf - gccf * gcaf;	Mc[0][2] = -gscf * gsbf * gcaf + gccf * gsaf;
		Mc[1][0] = gccf * gcbf;	Mc[1][1] = gccf * gsbf * gsaf - gscf * gcaf;	Mc[1][2] = gccf * gsbf * gcaf + gscf * gsaf;
//		Mc[2][0] = Mc[2][1] = Mc[2][2] = 0.0;

		e[0] = e[1] = e[2] = 0.0;		// zerolize error

		// compute gradient
		memset( grad, 0x00, 6 * sizeof( float ) );
		p1 = pa.pt1;
		p2 = pa.pt2;
		for( i = 0; i < pa.n; i++, p1++, p2++ ){
			// temp = ( M * p1 - p2 )T
			temp[0] = ( float )( M[0][0] * p1->ptr[0] + M[0][1] * p1->ptr[1] + M[0][2] * p1->ptr[2] + M[0][3] - p2->ptr[0] );
			temp[1] = ( float )( M[1][0] * p1->ptr[0] + M[1][1] * p1->ptr[1] + M[1][2] * p1->ptr[2] + M[1][3] - p2->ptr[1] );
			temp[2] = ( float )( M[2][0] * p1->ptr[0] + M[2][1] * p1->ptr[1] + M[2][2] * p1->ptr[2] + M[2][3] - p2->ptr[2] );
			e[0] += SQR( temp[0] ) + SQR( temp[1] ) + SQR( temp[2] );

			// gradient = tempT * [ e1 e2 e3 | Ma * p | Mb * p | Mc * p ]
			grad[0] += temp[0];
			grad[1] += temp[1];
			grad[2] += temp[2];
			grad[3] += ( float )(
				temp[0] * ( Ma[0][1] * p1->ptr[1] + Ma[0][2] * p1->ptr[2] ) +
				temp[1] * ( Ma[1][1] * p1->ptr[1] + Ma[1][2] * p1->ptr[2] ) +
				temp[2] * ( Ma[2][1] * p1->ptr[1] + Ma[2][2] * p1->ptr[2] ) );
			grad[4] += ( float )(
				temp[0] * ( Mb[0][0] * p1->ptr[0] + Mb[0][1] * p1->ptr[1] + Mb[0][2] * p1->ptr[2] ) +
				temp[1] * ( Mb[1][0] * p1->ptr[0] + Mb[1][1] * p1->ptr[1] + Mb[1][2] * p1->ptr[2] ) +
				temp[2] * ( Mb[2][0] * p1->ptr[0] + Mb[2][1] * p1->ptr[1] + Mb[2][2] * p1->ptr[2] ) );
			grad[5] += ( float )(
				temp[0] * ( Mc[0][0] * p1->ptr[0] + Mc[0][1] * p1->ptr[1] + Mc[0][2] * p1->ptr[2] ) +
				temp[1] * ( Mc[1][0] * p1->ptr[0] + Mc[1][1] * p1->ptr[1] + Mc[1][2] * p1->ptr[2] ) );
		}
		e[0] *= ni;					// normalize error
		if( e[0] < param.err ) break;		// break criteria

		gxf = MAX3( grad[3], grad[4], grad[5] );
		gyf = MIN3( grad[3], grad[4], grad[5] );
		a1 = ( float )MAX2( fabs( gxf ), fabs( gyf ) );
		a1 = ( a1 < M_PI ) ? param.lr1 * ni : param.lr1 * ( float )M_PI / a1;
		a2 = param.lr2 * a1;

		// compute error at a1
		buf.x = theta->x - a1 * grad[0];	buf.y = theta->y - a1 * grad[1];	buf.z = theta->z - a1 * grad[2]; 
		buf.a = theta->a - a1 * grad[3];	buf.b = theta->b - a1 * grad[4];	buf.c = theta->c - a1 * grad[5]; 
		CompRTMatf( M, buf );
		p1 = pa.pt1;
		p2 = pa.pt2;
		for( i = 0; i < pa.n; i++, p1++, p2++ ) {
			CompDist2fv( temp[0], M, *p1, *p2 );
			e[1] += temp[0];
		}
		e[1] *= ni;

		// compute error at a2 = 0.5
		buf.x = theta->x - a2 * grad[0];	buf.y = theta->y - a2 * grad[1];	buf.z = theta->z - a2 * grad[2]; 
		buf.a = theta->a - a2 * grad[3];	buf.b = theta->b - a2 * grad[4];	buf.c = theta->c - a2 * grad[5]; 
		CompRTMatf( M, buf );
		p1 = pa.pt1;
		p2 = pa.pt2;
		for( i = 0; i < pa.n; i++, p1++, p2++ ) {
			CompDist2fv( temp[0], M, *p1, *p2 );
			e[2] += temp[0];
		}
		e[2] *= ni;

		// find learning rate that minimize object function
		gxf = a1 * ( e[2] - e[0] ) - a2 * ( e[1] - e[0] );
		if( gxf == 0.0 ) gzf = a2;		// if solution is max, then set to solution, else set to min step
		else{
			gyf = SQR( a2 ) * ( e[1] - e[0] ) - SQR( a1 ) * ( e[2] - e[0] );
			gzf = -0.5f * gyf / gxf;
		}
		if( gzf < SMALL ) {
			if( e[0] < e[1] && e[0] < e[2] ) break;
			if( e[1] < e[2] ) gzf = a1;
			else gzf = a2;
		}
		// update agccording gradient
		theta->x -= gzf * grad[0];
		theta->y -= gzf * grad[1];
		theta->z -= gzf * grad[2];
		theta->a -= gzf * grad[3];
		theta->b -= gzf * grad[4];
		theta->c -= gzf * grad[5];
	}
	return e[0];
}

//
// PCA
//

void CompPCACoordf( Point3f_Array pa, float m[4][4] ){
	const int n = pa.n;
	Point3f *p3d;
	float p, u[3];
	float r[9], eigen[3], eig_vec[9];		// r: covariance matrix
	int i;

	// initialize
	memset( u, 0, 3 * sizeof( float ) );
	memset( r, 0, 9 * sizeof( float ) );

	// compute mean of points
	p3d = pa.pt;
	for( i = 0; i < n; i++, p3d++ ){
		u[0] += ( *p3d ).x;
		u[1] += ( *p3d ).y;
		u[2] += ( *p3d ).z;
	}
	p = 1.0f / ( float )n;
	u[0] *= p;
	u[1] *= p;
	u[2] *= p;

	// compute covariance matrix
	p3d = pa.pt;
	for( i = 0; i < n; i++, p3d++ ){
		// normalize
		gxf = ( *p3d ).x - u[0];
		gyf = ( *p3d ).y - u[1];
		gzf = ( *p3d ).z - u[2];
		p = 1.0f / ( float )sqrt( SQR( gxf ) + SQR( gyf ) + SQR( gzf ) );
		gxf *= p;
		gyf *= p;
		gzf *= p;
		// covariance matrix
		r[0] += SQR( gxf );
		r[1] += gxf * gyf;
		r[2] += gxf * gzf;
		r[4] += SQR( gyf );
		r[5] += gyf * gzf;
		r[8] += SQR( gzf );
	}
	p = 1.0f / ( float )n;
	r[0] *= p;
	r[3] = r[1] = r[1] * p;
	r[6] = r[2] = r[2] * p;
	r[4] *= p;
	r[7] = r[5] = r[5] * p;
	r[8] *= p;

	// solve eigen problem
	Jacobif( r, eigen, eig_vec, 3, 30 );
	// set x, y, z coordinate, x: max eigenvalue, y: second eigenvalue, z: min eigenvalue
	if( eigen[0] > eigen[1] ){			// 0 > 1
		if( eigen[0] > eigen[2] ){		// 0 > 2
			m[0][0] = eig_vec[0]; m[1][0] = eig_vec[3]; m[2][0] = eig_vec[6];
			if( eigen[1] > eigen[2] ){	// 0 > 1 > 2
				m[0][1] = eig_vec[1]; m[1][1] = eig_vec[4]; m[2][1] = eig_vec[7];
				m[0][2] = eig_vec[2]; m[1][2] = eig_vec[5]; m[2][2] = eig_vec[8];
			}
			else{						// 0 > 2 > 1
				m[0][1] = eig_vec[2]; m[1][1] = eig_vec[5]; m[2][1] = eig_vec[8];
				m[0][2] = eig_vec[1]; m[1][2] = eig_vec[4]; m[2][2] = eig_vec[7];
			}
		}
		else{							// 2 > 0 > 1
			m[0][0] = eig_vec[2]; m[1][0] = eig_vec[5]; m[2][0] = eig_vec[8];
			m[0][1] = eig_vec[0]; m[1][1] = eig_vec[3]; m[2][1] = eig_vec[6];
			m[0][2] = eig_vec[1]; m[1][2] = eig_vec[4]; m[2][2] = eig_vec[7];
		}
	}
	else{								// 1 > 0
		if( eigen[1] > eigen[2] ){		// 1 > 2
			m[0][0] = eig_vec[1]; m[1][0] = eig_vec[4]; m[2][0] = eig_vec[7];
			if( eigen[0] > eigen[2] ){	// 1 > 0 > 2
				m[0][1] = eig_vec[0]; m[1][1] = eig_vec[3]; m[2][1] = eig_vec[6];
				m[0][2] = eig_vec[2]; m[1][2] = eig_vec[5]; m[2][2] = eig_vec[8];
			}
			else{						// 1 > 2 > 0
				m[0][1] = eig_vec[2]; m[1][1] = eig_vec[5]; m[2][1] = eig_vec[8];
				m[0][2] = eig_vec[0]; m[1][2] = eig_vec[3]; m[2][2] = eig_vec[6];
			}
		}
		else{							// 2 > 1 > 0
			m[0][0] = eig_vec[2]; m[1][0] = eig_vec[5]; m[2][0] = eig_vec[8];
			m[0][1] = eig_vec[1]; m[1][1] = eig_vec[4]; m[2][1] = eig_vec[7];
			m[0][2] = eig_vec[0]; m[1][2] = eig_vec[3]; m[2][2] = eig_vec[6];
		}
	}
	if( ( m[1][0] * m[2][1] - m[2][0] * m[1][1] ) * m[0][2] < 0.0f ){	// not Euclidean coordinate
		m[0][2] *= -1.0f;	m[1][2] *= -1.0f;	m[2][2] *= -1.0f;
	}

	m[0][3] = u[0];	m[1][3] = u[1];	m[2][3] = u[2];
	m[3][0] = m[3][1] = m[3][2] = 0.0;
	m[3][3] = 1.0;
}

void CompPCACoordd( Point3d_Array pa, double m[4][4] ){
	const int n = pa.n;
	Point3d *p3d;
	double p, u[3];
	double r[9], eigen[3], eig_vec[9];		// r: covariance matrix
	int i;

	// initialize
	memset( u, 0, 3 * sizeof( double ) );
	memset( r, 0, 9 * sizeof( double ) );

	// compute mean of points
	p3d = pa.pt;
	for( i = 0; i < n; i++, p3d++ ){
		u[0] += ( *p3d ).x;
		u[1] += ( *p3d ).y;
		u[2] += ( *p3d ).z;
	}
	p = 1.0 / ( double )n;
	u[0] *= p;
	u[1] *= p;
	u[2] *= p;

	// compute covariance matrix
	p3d = pa.pt;
	for( i = 0; i < n; i++, p3d++ ){
		// normalize
		gxd = ( *p3d ).x - u[0];
		gyd = ( *p3d ).y - u[1];
		gzd = ( *p3d ).z - u[2];
		p = 1.0 / sqrt( SQR( gxd ) + SQR( gyd ) + SQR( gzd ) );
		gxd *= p;
		gyd *= p;
		gzd *= p;
		// covariance matrix
		r[0] += SQR( gxd );
		r[1] += gxd * gyd;
		r[2] += gxd * gzd;
		r[4] += SQR( gyd );
		r[5] += gyd * gzd;
		r[8] += SQR( gzd );
	}
	p = 1.0 / ( double )n;
	r[0] *= p;
	r[3] = r[1] = r[1] * p;
	r[6] = r[2] = r[2] * p;
	r[4] *= p;
	r[7] = r[5] = r[5] * p;
	r[8] *= p;

	// solve eigen problem
	Jacobid( r, eigen, eig_vec, 3, 30 );
	// set x, y, z coordinate, x: max eigenvalue, y: second eigenvalue, z: min eigenvalue
	if( eigen[0] > eigen[1] ){			// 0 > 1
		if( eigen[0] > eigen[2] ){		// 0 > 2
			m[0][0] = eig_vec[0]; m[1][0] = eig_vec[3]; m[2][0] = eig_vec[6];
			if( eigen[1] > eigen[2] ){	// 0 > 1 > 2
				m[0][1] = eig_vec[1]; m[1][1] = eig_vec[4]; m[2][1] = eig_vec[7];
				m[0][2] = eig_vec[2]; m[1][2] = eig_vec[5]; m[2][2] = eig_vec[8];
			}
			else{						// 0 > 2 > 1
				m[0][1] = eig_vec[2]; m[1][1] = eig_vec[5]; m[2][1] = eig_vec[8];
				m[0][2] = eig_vec[1]; m[1][2] = eig_vec[4]; m[2][2] = eig_vec[7];
			}
		}
		else{							// 2 > 0 > 1
			m[0][0] = eig_vec[2]; m[1][0] = eig_vec[5]; m[2][0] = eig_vec[8];
			m[0][1] = eig_vec[0]; m[1][1] = eig_vec[3]; m[2][1] = eig_vec[6];
			m[0][2] = eig_vec[1]; m[1][2] = eig_vec[4]; m[2][2] = eig_vec[7];
		}
	}
	else{								// 1 > 0
		if( eigen[1] > eigen[2] ){		// 1 > 2
			m[0][0] = eig_vec[1]; m[1][0] = eig_vec[4]; m[2][0] = eig_vec[7];
			if( eigen[0] > eigen[2] ){	// 1 > 0 > 2
				m[0][1] = eig_vec[0]; m[1][1] = eig_vec[3]; m[2][1] = eig_vec[6];
				m[0][2] = eig_vec[2]; m[1][2] = eig_vec[5]; m[2][2] = eig_vec[8];
			}
			else{						// 1 > 2 > 0
				m[0][1] = eig_vec[2]; m[1][1] = eig_vec[5]; m[2][1] = eig_vec[8];
				m[0][2] = eig_vec[0]; m[1][2] = eig_vec[3]; m[2][2] = eig_vec[6];
			}
		}
		else{							// 2 > 1 > 0
			m[0][0] = eig_vec[2]; m[1][0] = eig_vec[5]; m[2][0] = eig_vec[8];
			m[0][1] = eig_vec[1]; m[1][1] = eig_vec[4]; m[2][1] = eig_vec[7];
			m[0][2] = eig_vec[0]; m[1][2] = eig_vec[3]; m[2][2] = eig_vec[6];
		}
	}
	if( ( m[1][0] * m[2][1] - m[2][0] * m[1][1] ) * m[0][2] < 0.0 ){	// not Euclidean coordinate
		m[0][2] *= -1.0;	m[1][2] *= -1.0;	m[2][2] *= -1.0;
	}

	m[0][3] = u[0];	m[1][3] = u[1];	m[2][3] = u[2];
	m[3][0] = m[3][1] = m[3][2] = 0.0;
	m[3][3] = 1.0;
}

//
// Rigid Transform by SVD
//

void RigidTransform_SVDMatf( Point3f_Pair pa, float m[3][4] ){
	const float ni = 1.0f / ( float )pa.n;
	Point3f *p1, *p2, p1c, p2c;
	float ave[2][3], h[3][3], u[3][3], vT[3][3], w[3], hh[3][3];
	int i;

	// calculate mean
	memset( ave[0], 0x00, 6 * sizeof( float ) );
	p1 = pa.pt1;
	p2 = pa.pt2;
	for( i = 0; i < pa.n; i++, p1++, p2++ ){
		ave[0][0] += p1->x;	ave[0][1] += p1->y;	ave[0][2] += p1->z;
		ave[1][0] += p2->x;	ave[1][1] += p2->y;	ave[1][2] += p2->z;	
	}
	ave[0][0] *= ni;	ave[0][1] *= ni;	ave[0][2] *= ni;
	ave[1][0] *= ni;	ave[1][1] *= ni;	ave[1][2] *= ni;	

	// calculate matrix H = ( mi - m ) * ( di - d )
	p1 = pa.pt1;
	p2 = pa.pt2;
	memset( h[0], 0x00, 9 * sizeof( float ) );
	for( i = 0; i < pa.n; i++, p1++, p2++ ){
		p1c.x = p1->x - ave[0][0];
		p1c.y = p1->y - ave[0][1];
		p1c.z = p1->z - ave[0][2];
		p2c.x = p2->x - ave[1][0];
		p2c.y = p2->y - ave[1][1];
		p2c.z = p2->z - ave[1][2];
		h[0][0] += p1c.x * p2c.x;	h[0][1] += p1c.x * p2c.y;	h[0][2] += p1c.x * p2c.z;
		h[1][0] += p1c.y * p2c.x;	h[1][1] += p1c.y * p2c.y;	h[1][2] += p1c.y * p2c.z;
		h[2][0] += p1c.z * p2c.x;	h[2][1] += p1c.z * p2c.y;	h[2][2] += p1c.z * p2c.z;
	}
	// H = U * A * VT
	svdf( h[0], u[0], w, vT[0], 3, 3, 0.00001f, hh[0] );
	// desired R = V * UT
	h[0][0] = u[0][0] * vT[0][0] + u[0][1] * vT[1][0] + u[0][2] * vT[2][0];
	h[1][0] = u[0][0] * vT[0][1] + u[0][1] * vT[1][1] + u[0][2] * vT[2][1];
	h[2][0] = u[0][0] * vT[0][2] + u[0][1] * vT[1][2] + u[0][2] * vT[2][2];
	h[0][1] = u[1][0] * vT[0][0] + u[1][1] * vT[1][0] + u[1][2] * vT[2][0];
	h[1][1] = u[1][0] * vT[0][1] + u[1][1] * vT[1][1] + u[1][2] * vT[2][1];
	h[2][1] = u[1][0] * vT[0][2] + u[1][1] * vT[1][2] + u[1][2] * vT[2][2];
	h[0][2] = u[2][0] * vT[0][0] + u[2][1] * vT[1][0] + u[2][2] * vT[2][0];
	h[1][2] = u[2][0] * vT[0][1] + u[2][1] * vT[1][1] + u[2][2] * vT[2][1];
	h[2][2] = u[2][0] * vT[0][2] + u[2][1] * vT[1][2] + u[2][2] * vT[2][2];
	// det( R )
	w[0] = detf( h[0], 3 );
	if( w[0] > 0.0f ){	// desired R = V * UT
		m[0][0] = u[0][0] * vT[0][0] + u[0][1] * vT[1][0] + u[0][2] * vT[2][0];
		m[1][0] = u[0][0] * vT[0][1] + u[0][1] * vT[1][1] + u[0][2] * vT[2][1];
		m[2][0] = u[0][0] * vT[0][2] + u[0][1] * vT[1][2] + u[0][2] * vT[2][2];
		m[0][1] = u[1][0] * vT[0][0] + u[1][1] * vT[1][0] + u[1][2] * vT[2][0];
		m[1][1] = u[1][0] * vT[0][1] + u[1][1] * vT[1][1] + u[1][2] * vT[2][1];
		m[2][1] = u[1][0] * vT[0][2] + u[1][1] * vT[1][2] + u[1][2] * vT[2][2];
		m[0][2] = u[2][0] * vT[0][0] + u[2][1] * vT[1][0] + u[2][2] * vT[2][0];
		m[1][2] = u[2][0] * vT[0][1] + u[2][1] * vT[1][1] + u[2][2] * vT[2][1];
		m[2][2] = u[2][0] * vT[0][2] + u[2][1] * vT[1][2] + u[2][2] * vT[2][2];
	}
	else{						// desired R = V' * UT, V'=[v1 v2 -v3]
		m[0][0] = u[0][0] * vT[0][0] + u[0][1] * vT[1][0] + u[0][2] * -vT[2][0];
		m[1][0] = u[0][0] * vT[0][1] + u[0][1] * vT[1][1] + u[0][2] * -vT[2][1];
		m[2][0] = u[0][0] * vT[0][2] + u[0][1] * vT[1][2] + u[0][2] * -vT[2][2];
		m[0][1] = u[1][0] * vT[0][0] + u[1][1] * vT[1][0] + u[1][2] * -vT[2][0];
		m[1][1] = u[1][0] * vT[0][1] + u[1][1] * vT[1][1] + u[1][2] * -vT[2][1];
		m[2][1] = u[1][0] * vT[0][2] + u[1][1] * vT[1][2] + u[1][2] * -vT[2][2];
		m[0][2] = u[2][0] * vT[0][0] + u[2][1] * vT[1][0] + u[2][2] * -vT[2][0];
		m[1][2] = u[2][0] * vT[0][1] + u[2][1] * vT[1][1] + u[2][2] * -vT[2][1];
		m[2][2] = u[2][0] * vT[0][2] + u[2][1] * vT[1][2] + u[2][2] * -vT[2][2];
	}
	// desire T = p2 - R * p1
	m[0][3] = ave[1][0] - m[0][0] * ave[0][0] - m[0][1] * ave[0][1] - m[0][2] * ave[0][2];
	m[1][3] = ave[1][1] - m[1][0] * ave[0][0] - m[1][1] * ave[0][1] - m[1][2] * ave[0][2];
	m[2][3] = ave[1][2] - m[2][0] * ave[0][0] - m[2][1] * ave[0][1] - m[2][2] * ave[0][2];
}

void RigidTransform_SVDMatfv( Point3fv_Pair pa, float m[3][4] ){
#define x ptr[0]
#define y ptr[1]
#define z ptr[2]
	const float ni = 1.0f / ( float )pa.n;
	float buf[2][3];
	Point3fv *p1, *p2, p1c = { buf[0] }, p2c = { buf[1] };
	float ave[2][3], h[3][3], u[3][3], vT[3][3], w[3], hh[3][3];
	int i;

	// calculate mean
	memset( ave[0], 0x00, 6 * sizeof( float ) );
	p1 = pa.pt1;
	p2 = pa.pt2;
	for( i = 0; i < pa.n; i++, p1++, p2++ ){
		ave[0][0] += p1->x;	ave[0][1] += p1->y;	ave[0][2] += p1->z;
		ave[1][0] += p2->x;	ave[1][1] += p2->y;	ave[1][2] += p2->z;	
	}
	ave[0][0] *= ni;	ave[0][1] *= ni;	ave[0][2] *= ni;
	ave[1][0] *= ni;	ave[1][1] *= ni;	ave[1][2] *= ni;	

	// calculate matrix H = ( mi - m ) * ( di - d )
	p1 = pa.pt1;
	p2 = pa.pt2;
	memset( h[0], 0x00, 9 * sizeof( float ) );
	for( i = 0; i < pa.n; i++, p1++, p2++ ){
		p1c.x = p1->x - ave[0][0];
		p1c.y = p1->y - ave[0][1];
		p1c.z = p1->z - ave[0][2];
		p2c.x = p2->x - ave[1][0];
		p2c.y = p2->y - ave[1][1];
		p2c.z = p2->z - ave[1][2];
		h[0][0] += p1c.x * p2c.x;	h[0][1] += p1c.x * p2c.y;	h[0][2] += p1c.x * p2c.z;
		h[1][0] += p1c.y * p2c.x;	h[1][1] += p1c.y * p2c.y;	h[1][2] += p1c.y * p2c.z;
		h[2][0] += p1c.z * p2c.x;	h[2][1] += p1c.z * p2c.y;	h[2][2] += p1c.z * p2c.z;
	}
	// H = U * A * VT
	svdf( h[0], u[0], w, vT[0], 3, 3, 0.00001f, hh[0] );
	// desired R = V * UT
	h[0][0] = u[0][0] * vT[0][0] + u[0][1] * vT[1][0] + u[0][2] * vT[2][0];
	h[1][0] = u[0][0] * vT[0][1] + u[0][1] * vT[1][1] + u[0][2] * vT[2][1];
	h[2][0] = u[0][0] * vT[0][2] + u[0][1] * vT[1][2] + u[0][2] * vT[2][2];
	h[0][1] = u[1][0] * vT[0][0] + u[1][1] * vT[1][0] + u[1][2] * vT[2][0];
	h[1][1] = u[1][0] * vT[0][1] + u[1][1] * vT[1][1] + u[1][2] * vT[2][1];
	h[2][1] = u[1][0] * vT[0][2] + u[1][1] * vT[1][2] + u[1][2] * vT[2][2];
	h[0][2] = u[2][0] * vT[0][0] + u[2][1] * vT[1][0] + u[2][2] * vT[2][0];
	h[1][2] = u[2][0] * vT[0][1] + u[2][1] * vT[1][1] + u[2][2] * vT[2][1];
	h[2][2] = u[2][0] * vT[0][2] + u[2][1] * vT[1][2] + u[2][2] * vT[2][2];
	// det( R )
	w[0] = detf( h[0], 3 );
	if( w[0] > 0.0f ){	// desired R = V * UT
		m[0][0] = u[0][0] * vT[0][0] + u[0][1] * vT[1][0] + u[0][2] * vT[2][0];
		m[1][0] = u[0][0] * vT[0][1] + u[0][1] * vT[1][1] + u[0][2] * vT[2][1];
		m[2][0] = u[0][0] * vT[0][2] + u[0][1] * vT[1][2] + u[0][2] * vT[2][2];
		m[0][1] = u[1][0] * vT[0][0] + u[1][1] * vT[1][0] + u[1][2] * vT[2][0];
		m[1][1] = u[1][0] * vT[0][1] + u[1][1] * vT[1][1] + u[1][2] * vT[2][1];
		m[2][1] = u[1][0] * vT[0][2] + u[1][1] * vT[1][2] + u[1][2] * vT[2][2];
		m[0][2] = u[2][0] * vT[0][0] + u[2][1] * vT[1][0] + u[2][2] * vT[2][0];
		m[1][2] = u[2][0] * vT[0][1] + u[2][1] * vT[1][1] + u[2][2] * vT[2][1];
		m[2][2] = u[2][0] * vT[0][2] + u[2][1] * vT[1][2] + u[2][2] * vT[2][2];
	}
	else{						// desired R = V' * UT, V'=[v1 v2 -v3]
		m[0][0] = u[0][0] * vT[0][0] + u[0][1] * vT[1][0] + u[0][2] * -vT[2][0];
		m[1][0] = u[0][0] * vT[0][1] + u[0][1] * vT[1][1] + u[0][2] * -vT[2][1];
		m[2][0] = u[0][0] * vT[0][2] + u[0][1] * vT[1][2] + u[0][2] * -vT[2][2];
		m[0][1] = u[1][0] * vT[0][0] + u[1][1] * vT[1][0] + u[1][2] * -vT[2][0];
		m[1][1] = u[1][0] * vT[0][1] + u[1][1] * vT[1][1] + u[1][2] * -vT[2][1];
		m[2][1] = u[1][0] * vT[0][2] + u[1][1] * vT[1][2] + u[1][2] * -vT[2][2];
		m[0][2] = u[2][0] * vT[0][0] + u[2][1] * vT[1][0] + u[2][2] * -vT[2][0];
		m[1][2] = u[2][0] * vT[0][1] + u[2][1] * vT[1][1] + u[2][2] * -vT[2][1];
		m[2][2] = u[2][0] * vT[0][2] + u[2][1] * vT[1][2] + u[2][2] * -vT[2][2];
	}
	// desire T = p2 - R * p1
	m[0][3] = ave[1][0] - m[0][0] * ave[0][0] - m[0][1] * ave[0][1] - m[0][2] * ave[0][2];
	m[1][3] = ave[1][1] - m[1][0] * ave[0][0] - m[1][1] * ave[0][1] - m[1][2] * ave[0][2];
	m[2][3] = ave[1][2] - m[2][0] * ave[0][0] - m[2][1] * ave[0][1] - m[2][2] * ave[0][2];
#undef x
#undef y
#undef z
}