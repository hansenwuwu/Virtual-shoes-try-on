
#include "3dconfig.h"
#include <math.h>

//
// macro
//

#define SQR(x)			( (x) * (x) )
#define SMALL			0.00000001
#define M_PI			3.1415926535897932384626433832795
#define M_PIf			3.1415926535897932384626433832795f
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

//
// thread variables
//

__declspec( thread ) double gxd, gyd, gzd, gcad, gcbd, gccd, gsad, gsbd, gscd;
__declspec( thread ) float gxf, gyf, gzf, gcaf, gcbf, gccf, gsaf, gsbf, gscf;

//
// Functions
//

//
// Retate and Translation
//

void CompRTMati( float M[3][4], RT_Parami theta ){
	// compute cos ans sin
	gsaf = ( float )sin( DegToRec( ( float )theta.a ) );
	gsbf = ( float )sin( DegToRec( ( float )theta.b ) );
	gscf = ( float )sin( DegToRec( ( float )theta.c ) );
	gcaf = ( float )cos( DegToRec( ( float )theta.a ) );
	gcbf = ( float )cos( DegToRec( ( float )theta.b ) );
	gccf = ( float )cos( DegToRec( ( float )theta.c ) );
	
	// update matrix M
	M[0][0] = gccf * gcbf;	M[0][1] = gccf * gsbf * gsaf - gscf * gcaf;	M[0][2] = gccf * gsbf * gcaf + gscf * gsaf;	M[0][3] = ( float )theta.x;
	M[1][0] = gscf * gcbf;	M[1][1] = gscf * gsbf * gsaf + gccf * gcaf;	M[1][2] = gscf * gsbf * gcaf - gccf * gsaf;	M[1][3] = ( float )theta.y;
	M[2][0] = -gsbf;	M[2][1] = gcbf * gsaf;	M[2][2] = gcbf * gcaf;	M[2][3] = ( float )theta.z;
}

void CompRTMatf( float M[3][4], RT_Paramf theta ){
	// compute cos ans sin
	gsaf = ( float )sin( theta.a );
	gsbf = ( float )sin( theta.b );
	gscf = ( float )sin( theta.c );
	gcaf = ( float )cos( theta.a );
	gcbf = ( float )cos( theta.b );
	gccf = ( float )cos( theta.c );
	
	// update matrix M
	M[0][0] = gccf * gcbf;	M[0][1] = gccf * gsbf * gsaf - gscf * gcaf;	M[0][2] = gccf * gsbf * gcaf + gscf * gsaf;	M[0][3] = theta.x;
	M[1][0] = gscf * gcbf;	M[1][1] = gscf * gsbf * gsaf + gccf * gcaf;	M[1][2] = gscf * gsbf * gcaf - gccf * gsaf;	M[1][3] = theta.y;
	M[2][0] = -gsbf;	M[2][1] = gcbf * gsaf;	M[2][2] = gcbf * gcaf;	M[2][3] = theta.z;
}

void CompRTMatd( double M[3][4], RT_Paramd theta ){
	// compute cos ans sin
	gsad = sin( theta.a );
	gsbd = sin( theta.b );
	gscd = sin( theta.c );
	gcad = cos( theta.a );
	gcbd = cos( theta.b );
	gccd = cos( theta.c );
	
	// update matrix M
	M[0][0] = gccd * gcbd;	M[0][1] = gccd * gsbd * gsad - gscd * gcad;	M[0][2] = gccd * gsbd * gcad + gscd * gsad;	M[0][3] = theta.x;
	M[1][0] = gscd * gcbd;	M[1][1] = gscd * gsbd * gsad + gccd * gcad;	M[1][2] = gscd * gsbd * gcad - gccd * gsad;	M[1][3] = theta.y;
	M[2][0] = -gsbd;	M[2][1] = gcbd * gsad;	M[2][2] = gcbd * gcad;	M[2][3] = theta.z;
}

//
// Rotate at Point
//

void CompRPMati( float M[3][4], RP_Parami theta ){
	// compute cos ans sin
	gsaf = ( float )sin( DegToRec( ( float )theta.a ) );
	gsbf = ( float )sin( DegToRec( ( float )theta.b ) );
	gscf = ( float )sin( DegToRec( ( float )theta.c ) );
	gcaf = ( float )cos( DegToRec( ( float )theta.a ) );
	gcbf = ( float )cos( DegToRec( ( float )theta.b ) );
	gccf = ( float )cos( DegToRec( ( float )theta.c ) );
	gxf = ( float )theta.x;
	gyf = ( float )theta.y;
	gzf = ( float )theta.z;
	
	// update matrix M
	M[0][0] = gccf * gcbf;	M[0][1] = gccf * gsbf * gsaf - gscf * gcaf;	M[0][2] = gccf * gsbf * gcaf + gscf * gsaf;
	M[1][0] = gscf * gcbf;	M[1][1] = gscf * gsbf * gsaf + gccf * gcaf;	M[1][2] = gscf * gsbf * gcaf - gccf * gsaf;
	M[2][0] = -gsbf;	M[2][1] = gcbf * gsaf;	M[2][2] = gcbf * gcaf;
	M[0][3] = gxf - M[0][0] * gxf - M[0][1] * gyf - M[0][2] * gzf;
	M[1][3] = gyf - M[1][0] * gxf - M[1][1] * gyf - M[1][2] * gzf;
	M[2][3] = gzf - M[2][0] * gxf - M[2][1] * gyf - M[2][2] * gzf;
}

void CompRPMatf( float M[3][4], RP_Paramf theta ){
	// compute cos ans sin
	gsaf = ( float )sin( theta.a );
	gsbf = ( float )sin( theta.b );
	gscf = ( float )sin( theta.c );
	gcaf = ( float )cos( theta.a );
	gcbf = ( float )cos( theta.b );
	gccf = ( float )cos( theta.c );

	// update matrix M
	M[0][0] = gccf * gcbf;	M[0][1] = gccf * gsbf * gsaf - gscf * gcaf;	M[0][2] = gccf * gsbf * gcaf + gscf * gsaf;
	M[1][0] = gscf * gcbf;	M[1][1] = gscf * gsbf * gsaf + gccf * gcaf;	M[1][2] = gscf * gsbf * gcaf - gccf * gsaf;
	M[2][0] = -gsbf;	M[2][1] = gcbf * gsaf;	M[2][2] = gcbf * gcaf;
	M[0][3] = theta.x - M[0][0] * theta.x - M[0][1] * theta.y - M[0][2] * theta.z;
	M[1][3] = theta.y - M[1][0] * theta.x - M[1][1] * theta.y - M[1][2] * theta.z;
	M[2][3] = theta.z - M[2][0] * theta.x - M[2][1] * theta.y - M[2][2] * theta.z;
}

void CompRPMatd( double M[3][4], RP_Paramd theta ){
	// compute cos ans sin
	gsad = sin( theta.a );
	gsbd = sin( theta.b );
	gscd = sin( theta.c );
	gcad = cos( theta.a );
	gcbd = cos( theta.b );
	gccd = cos( theta.c );

	// update matrix M
	M[0][0] = gccd * gcbd;	M[0][1] = gccd * gsbd * gsad - gscd * gcad;	M[0][2] = gccd * gsbd * gcad + gscd * gsad;	
	M[1][0] = gscd * gcbd;	M[1][1] = gscd * gsbd * gsad + gccd * gcad;	M[1][2] = gscd * gsbd * gcad - gccd * gsad; 
	M[2][0] = -gsbd;	M[2][1] = gcbd * gsad;	M[2][2] = gcbd * gcad;
	M[0][3] = theta.x - M[0][0] * theta.x - M[0][1] * theta.y - M[0][2] * theta.z;
	M[1][3] = theta.y - M[1][0] * theta.x - M[1][1] * theta.y - M[1][2] * theta.z;
	M[2][3] = theta.z - M[2][0] * theta.x - M[2][1] * theta.y - M[2][2] * theta.z;
}

// 
// Axis Rotation
//

void CompRAMati( float M[3][4], RA_Parami theta ){
	gcaf = ( float )cos( DegToRec( ( float )theta.a ) );
	gsaf = ( float )sin( DegToRec( ( float )theta.a ) );
	gcbf = 1.0f - gcaf;	// 1 - cos
	// normalize vector
	gzf = 1.0f / ( float )sqrt( ( float )( SQR( theta.v.x ) + SQR( theta.v.y ) + SQR( theta.v.z ) ) );
	gxf = theta.v.x * gzf;
	gyf = theta.v.y * gzf;
	gzf = theta.v.z * gzf;
	// rotation
	M[0][0] = gcaf + SQR( gxf ) * gcbf;	M[0][1] = gxf * gyf * gcbf - gzf * gsaf;	M[0][2] = gxf * gyf * gcbf + gyf * gsaf;
	M[1][0] = gxf * gyf * gcbf + gzf * gsaf;	M[1][1] = gcaf + SQR( gyf ) * gcbf;	M[1][2] = gyf * gzf * gcbf - gxf * gsaf;
	M[2][0] = gxf * gzf * gcbf - gyf * gsaf;	M[2][1] = gyf * gzf * gcbf + gxf * gsaf;	M[2][2] = gcaf + SQR( gzf ) * gcbf;
	// position
	M[0][3] = theta.p.x - M[0][0] * theta.p.x - M[0][1] * theta.p.y - M[0][2] * theta.p.z;
	M[1][3] = theta.p.y - M[1][0] * theta.p.x - M[1][1] * theta.p.y - M[1][2] * theta.p.z;
	M[2][3] = theta.p.z - M[2][0] * theta.p.x - M[2][1] * theta.p.y - M[2][2] * theta.p.z;
}

void CompRAMatf( float M[3][4], RA_Paramf theta ){
	gcaf = ( float )cos( theta.a );
	gsaf = ( float )sin( theta.a );
	gcbf = 1.0f - gcaf;	// 1 - cos
	// normalize vector
	gzf = 1.0f / ( float )sqrt( SQR( theta.v.x ) + SQR( theta.v.y ) + SQR( theta.v.z ) );
	gxf = theta.v.x * gzf;
	gyf = theta.v.y * gzf;
	gzf = theta.v.z * gzf;
	// rotation
	M[0][0] = gcaf + SQR( gxf ) * gcbf;	M[0][1] = gxf * gyf * gcbf - gzf * gsaf;	M[0][2] = gxf * gyf * gcbf + gyf * gsaf;
	M[1][0] = gxf * gyf * gcbf + gzf * gsaf;	M[1][1] = gcaf + SQR( gyf ) * gcbf;	M[1][2] = gyf * gzf * gcbf - gxf * gsaf;
	M[2][0] = gxf * gzf * gcbf - gyf * gsaf;	M[2][1] = gyf * gzf * gcbf + gxf * gsaf;	M[2][2] = gcaf + SQR( gzf ) * gcbf;
	// position
	M[0][3] = theta.p.x - M[0][0] * theta.p.x - M[0][1] * theta.p.y - M[0][2] * theta.p.z;
	M[1][3] = theta.p.y - M[1][0] * theta.p.x - M[1][1] * theta.p.y - M[1][2] * theta.p.z;
	M[2][3] = theta.p.z - M[2][0] * theta.p.x - M[2][1] * theta.p.y - M[2][2] * theta.p.z;
}

void CompRAMatd( double M[3][4], RA_Paramd theta ){
	gcad = cos( theta.a );
	gsad = sin( theta.a );
	gcbd = 1.0 - gcad;	// 1 - cos
	// normalize vector
	gzd = 1.0 / sqrt( SQR( theta.v.x ) + SQR( theta.v.y ) + SQR( theta.v.z ) );
	gxd = theta.v.x * gzd;
	gyd = theta.v.y * gzd;
	gzd = theta.v.z * gzd;
	// rotation
	M[0][0] = gcad + SQR( gxd ) * gcbd;	M[0][1] = gxd * gyd * gcbd - gzd * gsad;	M[0][2] = gxd * gyd * gcbd + gyd * gsad;
	M[1][0] = gxd * gyd * gcbd + gzd * gsad;	M[1][1] = gcad + SQR( gyd ) * gcbd;	M[1][2] = gyd * gzd * gcbd - gxd * gsad;
	M[2][0] = gxd * gzd * gcbd - gyd * gsad;	M[2][1] = gyd * gzd * gcbd + gxd * gsad;	M[2][2] = gcad + SQR( gzd ) * gcbd;
	// position
	M[0][3] = theta.p.x - M[0][0] * theta.p.x - M[0][1] * theta.p.y - M[0][2] * theta.p.z;
	M[1][3] = theta.p.y - M[1][0] * theta.p.x - M[1][1] * theta.p.y - M[1][2] * theta.p.z;
	M[2][3] = theta.p.z - M[2][0] * theta.p.x - M[2][1] * theta.p.y - M[2][2] * theta.p.z;
}

//
// Rotate and Translate Parameter From Matrix
//

void CompRTFromMati( RT_Parami *theta, float m[3][4] ){
	// combute y-rotation = asin( -m[2][0] )
	if( m[2][0] >= 1.0f ) theta->b = 180;
	else if( m[2][0] <= -1.0f ) theta->b = -180;
	else theta->b = ( int )RecToDeg( asin( -m[2][0] ) );
	
	if( m[2][0] >= 1.0f ){
		// asin( m[0][1] )
		if( m[0][1] >= 1.0f ) gxf = 0.5f * M_PIf;
		else if( m[0][1] <= -1.0f ) gxf = -0.5f * M_PIf;
		else gxf = ( float )asin( m[0][1] );
		// acos( m[0][2] )
		if( m[0][2] >= 1.0f ) gyf = 0.0f;
		else if( m[0][2] <= -1.0f ) gyf = M_PIf;
		else gyf = ( float )acos( m[0][2] );
		theta->a = ( int )RecToDeg( 0.5f * ( gxf + gyf ) );
		theta->c = ( int )RecToDeg( 0.5f * ( gyf - gxf ) );
	}
	else if( m[2][0] <= -1.0f ){
		// asin( -m[0][1] )
		if( m[0][1] >= 1.0f ) gxf = -0.5f * M_PIf;
		else if( m[0][1] <= -1.0f ) gxf = 0.5f * M_PIf;
		else gxf = ( float )asin( -m[0][1] );
		// acos( -m[0][2] )
		if( m[0][2] >= 1.0f ) gyf = M_PIf;
		else if( m[0][2] <= -1.0f ) gyf = 0.0f;
		else gyf = ( float )acos( -m[0][2] );
		theta->a = ( int )RecToDeg( 0.5f * ( gxf + gyf ) );
		theta->c = ( int )RecToDeg( 0.5f * ( gxf - gyf ) );
	}
	else{
		// x-rotation = acos( m[2][2] / cos(b) )
		gxf = m[2][2] / ( float )cos( ( float )theta->b );
		if( gxf >= 1.0f ) theta->a = 0;
		else if( gxf <= -1.0f ) theta->a = 180;
		else theta->a = ( int )RecToDeg( ( float )acos( gxf ) );
		// y-rotation = acos( m[0][0] / cos(b) )
		gxf = m[0][0] / ( float )cos( ( float )theta->b );
		if( gxf >= 1.0f ) theta->c = 0;
		else if( gxf <= -1.0f ) theta->c = 180;
		else theta->c = ( int )RecToDeg( ( float )acos( gxf ) );
	}
	theta->x = ( int )m[0][3];	theta->y = ( int )m[1][3];	theta->z = ( int )m[2][3];
}

void CompRTFromMatf( RT_Paramf *theta, float m[3][4] ){
	// combute y-rotation = asin( -m[2][0] )
	if( m[2][0] >= 1.0f ) theta->b = 0.5f * M_PIf;
	else if( m[2][0] <= -1.0f ) theta->b = -0.5f * M_PIf;
	else theta->b = ( float )asin( -m[2][0] );
	
	if( m[2][0] >= 1.0f ){
		// asin( m[0][1] )
		if( m[0][1] >= 1.0f ) gxf = 0.5f * M_PIf;
		else if( m[0][1] <= -1.0f ) gxf = -0.5f * M_PIf;
		else gxf = ( float )asin( m[0][1] );
		// acos( m[0][2] )
		if( m[0][2] >= 1.0f ) gyf = 0.0f;
		else if( m[0][2] <= -1.0f ) gyf = M_PIf;
		else gyf = ( float )acos( m[0][2] );
		theta->a = 0.5f * ( gxf + gyf );
		theta->c = 0.5f * ( gyf - gxf );
	}
	else if( m[2][0] <= -1.0f ){
		// asin( -m[0][1] )
		if( m[0][1] >= 1.0f ) gxf = -0.5f * M_PIf;
		else if( m[0][1] <= -1.0f ) gxf = 0.5f * M_PIf;
		else gxf = ( float )asin( -m[0][1] );
		// acos( -m[0][2] )
		if( m[0][2] >= 1.0f ) gyf = M_PIf;
		else if( m[0][2] <= -1.0f ) gyf = 0.0f;
		else gyf = ( float )acos( -m[0][2] );
		theta->a = 0.5f * ( gxf + gyf );
		theta->c = 0.5f * ( gxf - gyf );
	}
	else{
		gcbf = ( float )cos( theta->b );
		// x-rotation = acos( m[2][2] / cos(b) )
		gxf = m[2][2] / gcbf;
		if( gxf >= 1.0f ) theta->a = 0.0f;
		else if( gxf <= -1.0f ) theta->a = M_PIf;
		else theta->a = ( float )acos( gxf );
		if( sin( theta->a ) * gcbf * m[2][1] < 0.0f ) theta->a = -theta->a;
		// y-rotation = acos( m[0][0] / cos(b) )
		gxf = m[0][0] / gcbf;
		if( gxf >= 1.0f ) theta->c = 0.0f;
		else if( gxf <= -1.0f ) theta->c = M_PIf;
		else theta->c = ( float )acos( gxf );
		if( sin( theta->c ) * gcbf * m[1][0] < 0.0f ) theta->c = -theta->c;
	}
	theta->x = m[0][3];	theta->y = m[1][3];	theta->z = m[2][3];
}

void CompRTFromMatd( RT_Paramd *theta, double m[3][4] ){
	// combute y-rotation = asin( -m[2][0] )
	if( m[2][0] >= 1.0 ) theta->b = 0.5 * M_PI;
	else if( m[2][0] <= -1.0 ) theta->b = -0.5 * M_PI;
	else theta->b = asin( -m[2][0] );
	
	if( m[2][0] >= 1.0 ){
		// asin( m[0][1] )
		if( m[0][1] >= 1.0 ) gxd = 0.5 * M_PI;
		else if( m[0][1] <= -1.0 ) gxd = -0.5 * M_PI;
		else gxd = asin( m[0][1] );
		// acos( m[0][2] )
		if( m[0][2] >= 1.0 ) gyd = 0.0;
		else if( m[0][2] <= -1.0 ) gyd = M_PI;
		else gyd = acos( m[0][2] );
		theta->a = 0.5 * ( gxd + gyd );
		theta->c = 0.5 * ( gyd - gxd );
	}
	else if( m[2][0] <= -1.0 ){
		// asin( -m[0][1] )
		if( m[0][1] >= 1.0 ) gxd = -0.5 * M_PI;
		else if( m[0][1] <= -1.0 ) gxd = 0.5 * M_PI;
		else gxd = asin( -m[0][1] );
		// acos( -m[0][2] )
		if( m[0][2] >= 1.0 ) gyd = M_PI;
		else if( m[0][2] <= -1.0 ) gyd = 0.0;
		else gyd = acos( -m[0][2] );
		theta->a = 0.5 * ( gxd + gyd );
		theta->c = 0.5 * ( gxd - gyd );
	}
	else{
		gcbd = cos( theta->b );
		// x-rotation = acos( m[2][2] / cos(b) )
		gxd = m[2][2] / gcbd;
		if( gxd >= 1.0 ) theta->a = 0.0;
		else if( gxd <= -1.0 ) theta->a = M_PI;
		else theta->a = acos( gxd );
		if( sin( theta->a ) * gcbd * m[2][1] < 0.0 ) theta->a = -theta->a;
		// y-rotation = acos( m[0][0] / cos(b) )
		gxd = m[0][0] / gcbd;
		if( gxd >= 1.0 ) theta->c = 0.0;
		else if( gxd <= -1.0 ) theta->c = M_PI;
		else theta->c = acos( gxd );
		if( sin( theta->c ) * gcbd * m[1][0] < 0.0 ) theta->c = -theta->c;
	}
	theta->x = m[0][3];	theta->y = m[1][3];	theta->z = m[2][3];
}



//
// Transformation
//

void MultiPointf( Point3f *p, float m[3][4] ){
	gxf = m[0][0] * p->x + m[0][1] * p->y + m[0][2] * p->z + m[0][3];
	gyf = m[1][0] * p->x + m[1][1] * p->y + m[1][2] * p->z + m[1][3];
	gzf = m[2][0] * p->x + m[2][1] * p->y + m[2][2] * p->z + m[2][3];
	p->x = gxf;
	p->y = gyf;
	p->z = gzf;
}

void MultiPointd( Point3d *p, double m[3][4] ){
	gxd = m[0][0] * p->x + m[0][1] * p->y + m[0][2] * p->z + m[0][3];
	gyd = m[1][0] * p->x + m[1][1] * p->y + m[1][2] * p->z + m[1][3];
	gzd = m[2][0] * p->x + m[2][1] * p->y + m[2][2] * p->z + m[2][3];
	p->x = gxd;
	p->y = gyd;
	p->z = gzd;
}

void MultiPointfv( Point3fv *p, float m[3][4] ){
	gxf = m[0][0] * p->ptr[0] + m[0][1] * p->ptr[1] + m[0][2] * p->ptr[2] + m[0][3];
	gyf = m[1][0] * p->ptr[0] + m[1][1] * p->ptr[1] + m[1][2] * p->ptr[2] + m[1][3];
	gzf = m[2][0] * p->ptr[0] + m[2][1] * p->ptr[1] + m[2][2] * p->ptr[2] + m[2][3];
	p->ptr[0] = gxf;
	p->ptr[1] = gyf;
	p->ptr[2] = gzf;
}

void MultiPointdv( Point3dv *p, double m[3][4] ){
	gxd = m[0][0] * p->ptr[0] + m[0][1] * p->ptr[1] + m[0][2] * p->ptr[2] + m[0][3];
	gyd = m[1][0] * p->ptr[0] + m[1][1] * p->ptr[1] + m[1][2] * p->ptr[2] + m[1][3];
	gzd = m[2][0] * p->ptr[0] + m[2][1] * p->ptr[1] + m[2][2] * p->ptr[2] + m[2][3];
	p->ptr[0] = gxd;
	p->ptr[1] = gyd;
	p->ptr[2] = gzd;
}

void MultiMatf( float c[3][4], const float a[3][4], const float b[3][4] ){
	c[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0];
	c[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1];
	c[0][2] = a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2];
	c[0][3] = a[0][0] * b[0][3] + a[0][1] * b[1][3] + a[0][2] * b[2][3] + a[0][3];
	c[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0];
	c[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1];
	c[1][2] = a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2];
	c[1][3] = a[1][0] * b[0][3] + a[1][1] * b[1][3] + a[1][2] * b[2][3] + a[1][3];
	c[2][0] = a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0];
	c[2][1] = a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1];
	c[2][2] = a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2];
	c[2][3] = a[2][0] * b[0][3] + a[2][1] * b[1][3] + a[2][2] * b[2][3] + a[2][3];
}

void MultiMatd( double c[3][4], const double a[3][4], const double b[3][4] ){
	c[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0];
	c[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1];
	c[0][2] = a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2];
	c[0][3] = a[0][0] * b[0][3] + a[0][1] * b[1][3] + a[0][2] * b[2][3] + a[0][3];
	c[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0];
	c[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1];
	c[1][2] = a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2];
	c[1][3] = a[1][0] * b[0][3] + a[1][1] * b[1][3] + a[1][2] * b[2][3] + a[1][3];
	c[2][0] = a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0];
	c[2][1] = a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1];
	c[2][2] = a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2];
	c[2][3] = a[2][0] * b[0][3] + a[2][1] * b[1][3] + a[2][2] * b[2][3] + a[2][3];
}

void TransformPointf( Point3f *p, RT_Paramf theta ){
	// compute cos ans sin
	gsaf = ( float )sin( theta.a );
	gsbf = ( float )sin( theta.b );
	gscf = ( float )sin( theta.c );
	gcaf = ( float )cos( theta.a );
	gcbf = ( float )cos( theta.b );
	gccf = ( float )cos( theta.c );
	
	// update matrix M
	gxf = ( gccf * gcbf ) * p->x + ( gccf * gsbf * gsaf - gscf * gcaf ) * p->y + ( gccf * gsbf * gcaf + gscf * gsaf ) * p->z + theta.x;
	gyf = ( gscf * gcbf ) * p->x + ( gscf * gsbf * gsaf + gccf * gcaf ) * p->y + ( gscf * gsbf * gcaf - gccf * gsaf ) * p->z + theta.y;
	gzf = ( -gsbf ) * p->x + ( gcbf * gsaf ) * p->y + ( gcbf * gcaf ) * p->z + theta.z;
	p->x = gxf;
	p->y = gyf;
	p->z = gzf;
}

void TransformPointd( Point3d *p, RT_Paramd theta ){
	// compute cos ans sin
	gsad = sin( theta.a );
	gsbd = sin( theta.b );
	gscd = sin( theta.c );
	gcad = cos( theta.a );
	gcbd = cos( theta.b );
	gccd = cos( theta.c );
	
	// update matrix M
	gxd = ( gccd * gcbd ) * p->x + ( gccd * gsbd * gsad - gscd * gcad ) * p->y + ( gccd * gsbd * gcad + gscd * gsad ) * p->z + theta.x;
	gyd = ( gscd * gcbd ) * p->x + ( gscd * gsbd * gsad + gccd * gcad ) * p->y + ( gscd * gsbd * gcad - gccd * gsad ) * p->z + theta.y;
	gzd = ( -gsbd ) * p->x + ( gcbd * gsad ) * p->y + ( gcbd * gcad ) * p->z + theta.z;
	p->x = gxd;
	p->y = gyd;
	p->z = gzd;
}

void TransformPointfv( Point3fv *p, RT_Paramf theta ){
	// compute cos ans sin
	gsaf = ( float )sin( theta.a );
	gsbf = ( float )sin( theta.b );
	gscf = ( float )sin( theta.c );
	gcaf = ( float )cos( theta.a );
	gcbf = ( float )cos( theta.b );
	gccf = ( float )cos( theta.c );
	
	// update matrix M
	gxf = ( gccf * gcbf ) * p->ptr[0] + ( gccf * gsbf * gsaf - gscf * gcaf ) * p->ptr[1] + ( gccf * gsbf * gcaf + gscf * gsaf ) * p->ptr[2] + theta.x;
	gyf = ( gscf * gcbf ) * p->ptr[0] + ( gscf * gsbf * gsaf + gccf * gcaf ) * p->ptr[1] + ( gscf * gsbf * gcaf - gccf * gsaf ) * p->ptr[2] + theta.y;
	gzf = ( -gsbf ) * p->ptr[0] + ( gcbf * gsaf ) * p->ptr[1] + ( gcbf * gcaf ) * p->ptr[2] + theta.z;
	p->ptr[0] = gxf;
	p->ptr[1] = gyf;
	p->ptr[2] = gzf;
}

void TransformPointdv( Point3dv *p, RT_Paramd theta ){
	// compute cos ans sin
	gsad = sin( theta.a );
	gsbd = sin( theta.b );
	gscd = sin( theta.c );
	gcad = cos( theta.a );
	gcbd = cos( theta.b );
	gccd = cos( theta.c );
	
	// update matrix M
	gxd = ( gccd * gcbd ) * p->ptr[0] + ( gccd * gsbd * gsad - gscd * gcad ) * p->ptr[1] + ( gccd * gsbd * gcad + gscd * gsad ) * p->ptr[2] + theta.x;
	gyd = ( gscd * gcbd ) * p->ptr[0] + ( gscd * gsbd * gsad + gccd * gcad ) * p->ptr[1] + ( gscd * gsbd * gcad - gccd * gsad ) * p->ptr[2] + theta.y;
	gzd = ( -gsbd ) * p->ptr[0] + ( gcbd * gsad ) * p->ptr[1] + ( gcbd * gcad ) * p->ptr[2] + theta.z;
	p->ptr[0] = gxd;
	p->ptr[1] = gyd;
	p->ptr[2] = gzd;
}

//
// Inverse
//

void InvRTMatf( float out[3][4], const float in[3][4] ){
	out[0][0] = in[0][0];	out[0][1] = in[1][0];	out[0][2] = in[2][0];
	out[1][0] = in[0][1];	out[1][1] = in[1][1];	out[1][2] = in[2][1];
	out[2][0] = in[0][2];	out[2][1] = in[1][2];	out[2][2] = in[2][2];
	out[0][3] = -( in[0][0] * in[0][3] + in[1][0] * in[1][3] + in[2][0] * in[2][3] );
	out[1][3] = -( in[0][1] * in[0][3] + in[1][1] * in[1][3] + in[2][1] * in[2][3] );
	out[2][3] = -( in[0][2] * in[0][3] + in[1][2] * in[1][3] + in[2][2] * in[2][3] );
}

void InvRTMatd( double out[3][4], const double in[3][4] ){
	out[0][0] = in[0][0];	out[0][1] = in[1][0];	out[0][2] = in[2][0];
	out[1][0] = in[0][1];	out[1][1] = in[1][1];	out[1][2] = in[2][1];
	out[2][0] = in[0][2];	out[2][1] = in[1][2];	out[2][2] = in[2][2];
	out[0][3] = -( in[0][0] * in[0][3] + in[1][0] * in[1][3] + in[2][0] * in[2][3] );
	out[1][3] = -( in[0][1] * in[0][3] + in[1][1] * in[1][3] + in[2][1] * in[2][3] );
	out[2][3] = -( in[0][2] * in[0][3] + in[1][2] * in[1][3] + in[2][2] * in[2][3] );
}

//
// Basic Linear Algebra
//

// Dot: a * b

int Doti( struct _point_3d_int a, struct _point_3d_int b ){
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

float Dotf( struct _point_3d_float a, struct _point_3d_float b ){
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

double Dotd( struct _point_3d_double a, struct _point_3d_double b ){
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

int Dotiv( struct _point_3d_int_ptr a, struct _point_3d_int_ptr b ){
	return a.ptr[0] * b.ptr[0] + a.ptr[1] * b.ptr[1] + a.ptr[2] * b.ptr[2];
}

float Dotfptr( struct _point_3d_float_ptr a, struct _point_3d_float_ptr b ){
	return a.ptr[0] * b.ptr[0] + a.ptr[1] * b.ptr[1] + a.ptr[2] * b.ptr[2];
}

double Dotdptr( struct _point_3d_double_ptr a, struct _point_3d_double_ptr b ){
	return a.ptr[0] * b.ptr[0] + a.ptr[1] * b.ptr[1] + a.ptr[2] * b.ptr[2];
}

// Cross: c = a x b

void Crossi( struct _point_3d_int a, struct _point_3d_int b, struct _point_3d_int *c ){
	(*c).x = a.y * b.z - a.z * b.y;
	(*c).y = a.z * b.x - a.x * b.z;
	(*c).z = a.x * b.y - a.y * b.x;
}

void Crossf( struct _point_3d_float a, struct _point_3d_float b, struct _point_3d_float *c ){
	(*c).x = a.y * b.z - a.z * b.y;
	(*c).y = a.z * b.x - a.x * b.z;
	(*c).z = a.x * b.y - a.y * b.x;
}

void Crossd( struct _point_3d_double a, struct _point_3d_double b, struct _point_3d_double *c ){
	(*c).x = a.y * b.z - a.z * b.y;
	(*c).y = a.z * b.x - a.x * b.z;
	(*c).z = a.x * b.y - a.y * b.x;
}

void Crossiptr( struct _point_3d_int_ptr a, struct _point_3d_int_ptr b, struct _point_3d_int_ptr c ){
	c.ptr[0] = a.ptr[1] * b.ptr[2] - a.ptr[2] * b.ptr[1];
	c.ptr[1] = a.ptr[2] * b.ptr[0] - a.ptr[0] * b.ptr[2];
	c.ptr[2] = a.ptr[0] * b.ptr[1] - a.ptr[1] * b.ptr[0];
}

void Crossfptr( struct _point_3d_float_ptr a, struct _point_3d_float_ptr b, struct _point_3d_float_ptr c ){
	c.ptr[0] = a.ptr[1] * b.ptr[2] - a.ptr[2] * b.ptr[1];
	c.ptr[1] = a.ptr[2] * b.ptr[0] - a.ptr[0] * b.ptr[2];
	c.ptr[2] = a.ptr[0] * b.ptr[1] - a.ptr[1] * b.ptr[0];
}

void Crossdptr( struct _point_3d_double_ptr a, struct _point_3d_double_ptr b, struct _point_3d_double_ptr c ){
	c.ptr[0] = a.ptr[1] * b.ptr[2] - a.ptr[2] * b.ptr[1];
	c.ptr[1] = a.ptr[2] * b.ptr[0] - a.ptr[0] * b.ptr[2];
	c.ptr[2] = a.ptr[0] * b.ptr[1] - a.ptr[1] * b.ptr[0];
}


// out[16] = A[16] * B[16]
void Mult16Matf(float out[16],float A[16],float B[16])
{
	out[ 0] =A[ 0] *B[ 0] +A[ 1]*B[ 4] +A[ 2] *B[ 8] +A[ 3]*B[12];
	out[ 1] =A[ 0] *B[ 1] +A[ 1]*B[ 5] +A[ 2] *B[ 9] +A[ 3]*B[13];
	out[ 2] =A[ 0] *B[ 2] +A[ 1]*B[ 6] +A[ 2] *B[10] +A[ 3]*B[14];
	out[ 3] =A[ 0] *B[ 3] +A[ 1]*B[ 7] +A[ 2] *B[11] +A[ 3]*B[15];
	out[ 4] =A[ 4] *B[ 0] +A[ 5]*B[ 4] +A[ 6] *B[ 8] +A[ 7]*B[12];
	out[ 5] =A[ 4] *B[ 1] +A[ 5]*B[ 5] +A[ 6] *B[ 9] +A[ 7]*B[13];
	out[ 6] =A[ 4] *B[ 2] +A[ 5]*B[ 6] +A[ 6] *B[10] +A[ 7]*B[14];
	out[ 7] =A[ 4] *B[ 3] +A[ 5]*B[ 7] +A[ 6] *B[11] +A[ 7]*B[15];
	out[ 8] =A[ 8] *B[ 0] +A[ 9]*B[ 4] +A[10] *B[ 8] +A[11]*B[12];
	out[ 9] =A[ 8] *B[ 1] +A[ 9]*B[ 5] +A[10] *B[ 9] +A[11]*B[13];
	out[10] =A[ 8] *B[ 2] +A[ 9]*B[ 6] +A[10] *B[10] +A[11]*B[14];
	out[11] =A[ 8] *B[ 3] +A[ 9]*B[ 7] +A[10] *B[11] +A[11]*B[15];
	out[12] =A[12] *B[ 0] +A[13]*B[ 4] +A[14] *B[ 8] +A[15]*B[12];
	out[13] =A[12] *B[ 1] +A[13]*B[ 5] +A[14] *B[ 9] +A[15]*B[13];
	out[14] =A[12] *B[ 2] +A[13]*B[ 6] +A[14] *B[10] +A[15]*B[14];
	out[15] =A[12] *B[ 3] +A[13]*B[ 7] +A[14] *B[11] +A[15]*B[15];
}


void Trans16Matf(float m[16])
{
	float tempqq[16];
	tempqq[ 0]=m[ 0];	tempqq[ 1]=m[ 4];	tempqq[ 2]=m[ 8];	tempqq[ 3]=m[12];
	tempqq[ 4]=m[ 1];	tempqq[ 5]=m[ 5];	tempqq[ 6]=m[ 9];	tempqq[ 7]=m[13];
	tempqq[ 8]=m[ 2];	tempqq[ 9]=m[ 6];	tempqq[10]=m[10];	tempqq[11]=m[14];
	tempqq[12]=m[ 3];	tempqq[13]=m[ 7];	tempqq[14]=m[11];	tempqq[15]=m[15];
	memcpy(m,tempqq,16*sizeof(float));


}