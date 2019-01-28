
#include "depthimage.h"
#include "track.h"
#include <string.h>
#include <math.h>
#include <float.h>
#include "numerical\linear.h"

//
// Macro
//

#define SQR( x )	( (x) * (x) )

//
// Class Functions
//

DepthImage::DepthImage( void ): ImageR< int >(), CAMERA_PARAM(){
	_color = 1;
}

int DepthImage::Initialize( int w, int h, int color ){
	return this->ImageR< int >::Initialize( w, h, 1 );	// must be gray-scale
}

void DepthImage::Clear( void ){
	this->ImageR< int >::Clear();
	this->CAMERA_PARAM::Clear();
}

int DepthImage::ReadFile( const char *path ){
	if( this->ImageR< int >::ReadFile( path ) == 0 ) return 0;
	return SetColor( 1 );			// change to gray-scale
}

int DepthImage::ReadBMP( const char *path ){
	if( this->ImageR< int >::ReadBMP( path ) == 0 ) return 0;
	return SetColor( 1 );			// change to gray-scale
}

int DepthImage::ReadJPEG( const char *path ){
	if( this->ImageR< int >::ReadJPEG( path ) == 0 ) return 0;
	return SetColor( 1 );			// change to gray-scale
}

int DepthImage::ReadPPM( const char *path ){
	if( this->ImageR< int >::ReadPPM( path ) == 0 ) return 0;
	return SetColor( 1 );			// change to gray-scale
}

int DepthImage::SetImage( int *bits, int w, int h, int color ){
	this->ImageR< int >::SetImage( bits, w, h, color );
	return SetColor( 1 );			// change to gray-scale
}

// 3D informations

inline Point3f DepthImage::Comp3DPosf( int u, int v ){
	Point3f buf;
	buf.z = ( float )_bits[ v * _stride + u ];
	buf.x = ( float )( ( ( float )u - _u0 ) * _beta - ( ( float )v - _v0 ) * _gamma );
	buf.x *= ( float )( buf.z / ( _alpha * _beta ) );
	buf.y = ( float )( ( ( float )v - _v0 ) * buf.z / _beta );
	buf.z *= -1.0f;

	/*!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
	//buf.x = -buf.x;
	return buf;
}

inline Point3d DepthImage::Comp3DPosd( int u, int v ){
	Point3d buf;
	buf.z = ( double )_bits[ v * _stride + u ];
	buf.x = ( ( ( double )u - _u0 ) * _beta - ( ( double )v - _v0 ) * _gamma );
	buf.x *= buf.z / ( _alpha * _beta );
	buf.y = ( ( double )v - _v0 ) * buf.z / _beta;
	buf.z *= -1.0;
	return buf;
}

inline Point3f DepthImage::CompNormalf( int x, int y ){
	Point3f p[3][3], n, du, dv;

	// out of range
	if( x < 1 || y < 1 || x >= this->ImageR< int >::_w - 1 || y >= this->ImageR< int >::_h - 1 ) {
		memset( &( p[0][0] ), 0, sizeof( Point3f ) );
		return p[0][0];
	}
	
	// get 3D coordinate around ( x, y )
	p[0][0] = Comp3DPosf( x - 1, y - 1 );
	p[1][0] = Comp3DPosf( x, y - 1 );
	p[2][0] = Comp3DPosf( x + 1, y - 1 );
	p[0][1] = Comp3DPosf( x - 1, y );
	p[1][1] = Comp3DPosf( x, y );
	p[2][1] = Comp3DPosf( x + 1, y );
	p[0][2] = Comp3DPosf( x - 1, y + 1 );
	p[1][2] = Comp3DPosf( x, y + 1 );
	p[2][2] = Comp3DPosf( x + 1, y + 1 );

	// Sobel operation
	du.x = 0.125f * ( p[2][0].x + p[2][2].x + p[2][1].x + p[2][1].x - p[0][0].x - p[0][2].x - p[0][1].x - p[0][1].x );
	du.y = 0.125f * ( p[2][0].y + p[2][2].y + p[2][1].y + p[2][1].y - p[0][0].y - p[0][2].y - p[0][1].y - p[0][1].y );
	du.z = 0.125f * ( p[2][0].z + p[2][2].z + p[2][1].z + p[2][1].z - p[0][0].z - p[0][2].z - p[0][1].z - p[0][1].z );
	dv.x = 0.125f * ( p[0][2].x + p[2][2].x + p[1][2].x + p[1][2].x - p[0][0].x - p[2][0].x - p[1][0].x - p[1][0].x );
	dv.y = 0.125f * ( p[0][2].y + p[2][2].y + p[1][2].y + p[1][2].y - p[0][0].y - p[2][0].y - p[1][0].y - p[1][0].y );
	dv.z = 0.125f * ( p[0][2].z + p[2][2].z + p[1][2].z + p[1][2].z - p[0][0].z - p[2][0].z - p[1][0].z - p[1][0].z );

	Crossf( du, dv, &n );
	return n;
}

inline Point3d DepthImage::CompNormald( int x, int y ){
	Point3d p[3][3], n, du, dv;

	// out of range
	if( x < 1 || y < 1 || x >= this->ImageR< int >::_w - 1 || y >= this->ImageR< int >::_h - 1 ) {
		memset( &( p[0][0] ), 0, sizeof( Point3d ) );
		return p[0][0];
	}
	
	// get 3D coordinate around ( x, y )
	p[0][0] = Comp3DPosd( x - 1, y - 1 );
	p[1][0] = Comp3DPosd( x, y - 1 );
	p[2][0] = Comp3DPosd( x + 1, y - 1 );
	p[0][1] = Comp3DPosd( x - 1, y );
	p[1][1] = Comp3DPosd( x, y );
	p[2][1] = Comp3DPosd( x + 1, y );
	p[0][2] = Comp3DPosd( x - 1, y + 1 );
	p[1][2] = Comp3DPosd( x, y + 1 );
	p[2][2] = Comp3DPosd( x + 1, y + 1 );

	// Sobel operation
	du.x = 0.125 * ( p[2][0].x + p[2][2].x + p[2][1].x + p[2][1].x - p[0][0].x - p[0][2].x - p[0][1].x - p[0][1].x );
	du.y = 0.125 * ( p[2][0].y + p[2][2].y + p[2][1].y + p[2][1].y - p[0][0].y - p[0][2].y - p[0][1].y - p[0][1].y );
	du.z = 0.125 * ( p[2][0].z + p[2][2].z + p[2][1].z + p[2][1].z - p[0][0].z - p[0][2].z - p[0][1].z - p[0][1].z );
	dv.x = 0.125 * ( p[0][2].x + p[2][2].x + p[1][2].x + p[1][2].x - p[0][0].x - p[2][0].x - p[1][0].x - p[1][0].x );
	dv.y = 0.125 * ( p[0][2].y + p[2][2].y + p[1][2].y + p[1][2].y - p[0][0].y - p[2][0].y - p[1][0].y - p[1][0].y );
	dv.z = 0.125 * ( p[0][2].z + p[2][2].z + p[1][2].z + p[1][2].z - p[0][0].z - p[2][0].z - p[1][0].z - p[1][0].z );

	Crossd( du, dv, &n );
	return n;
}

int DepthImage::Curvature( double *k1, double *k2, int x, int y, int r ){
	Point3d p[3][3], n, du, dv, duu, duv, dvv;
	double L, M, N, b, d;

	// out of range
	if( x < r || y < r || x >= this->ImageR< int >::_w - r || y >= this->ImageR< int >::_h - r ) return 0;
	
	// get 3D coordinate around ( x, y )
	p[0][0] = Comp3DPosd( x - r, y - r );
	p[1][0] = Comp3DPosd( x, y - r );
	p[2][0] = Comp3DPosd( x + r, y - r );
	p[0][1] = Comp3DPosd( x - r, y );
	p[1][1] = Comp3DPosd( x, y );
	p[2][1] = Comp3DPosd( x + r, y );
	p[0][2] = Comp3DPosd( x - r, y + r );
	p[1][2] = Comp3DPosd( x, y + r );
	p[2][2] = Comp3DPosd( x + r, y + r );

	// compute du, dv, duu, duv, dvv from 9 points
	d = 0.5 / ( double )r;
	du.x = d * ( p[2][1].x - p[0][1].x );
	du.y = d * ( p[2][1].y - p[0][1].y );
	du.z = d * ( p[2][1].z - p[0][1].z );
	dv.x = d * ( p[1][2].x - p[1][0].x );
	dv.y = d * ( p[1][2].y - p[1][0].y );
	dv.z = d * ( p[1][2].z - p[1][0].z );
	d = SQR( 1.0 / ( double )r );
	duu.x = d * ( p[2][1].x + p[0][1].x - p[1][1].x - p[1][1].x );
	duu.y = d * ( p[2][1].y + p[0][1].y - p[1][1].y - p[1][1].y );
	duu.z = d * ( p[2][1].z + p[0][1].z - p[1][1].z - p[1][1].z );
	dvv.x = d * ( p[1][2].x + p[1][0].x - p[1][1].x - p[1][1].x );
	dvv.y = d * ( p[1][2].y + p[1][0].y - p[1][1].y - p[1][1].y );
	dvv.z = d * ( p[1][2].z + p[1][0].z - p[1][1].z - p[1][1].z );
	d = SQR( 0.5 / ( double )r );
	duv.x = d * ( p[2][2].x + p[0][0].x - p[2][0].x - p[0][2].x );
	duv.y = d * ( p[2][2].y + p[0][0].y - p[2][0].y - p[0][2].y );
	duv.z = d * ( p[2][2].z + p[0][0].z - p[2][0].z - p[0][2].z );

	// compute n = u x v
	Crossd( du, dv, &n );
	
	// second fundamental form 
	L = Dotd( duu, n );
	M = Dotd( duv, n );
	N = Dotd( dvv, n );

	d = L - M;
	d = 0.5 * sqrt( SQR( d ) + 4.0 * SQR( M ) );	// d = sqrt( b2 - 4ac )
	b = -0.5 * ( L + M );

	*k1 = b + d;
	*k2 = b - d;
	return 1;
}

double DepthImage::GaussCurvature( int x, int y ){
	Point3d p[3][3], n, du, dv, duu, duv, dvv;
	double L, M, N;

	// out of range
	if( x <= 0 || y <= 0 || x > this->ImageR< int >::_w - 2 || y > this->ImageR< int >::_h - 2 ) return 0.0;
	
	// get 3D coordinate around ( x, y )
	p[0][0] = Comp3DPosd( x - 1, y - 1 );
	p[1][0] = Comp3DPosd( x, y - 1 );
	p[2][0] = Comp3DPosd( x + 1, y - 1 );
	p[0][1] = Comp3DPosd( x - 1, y );
	p[1][1] = Comp3DPosd( x, y );
	p[2][1] = Comp3DPosd( x + 1, y );
	p[0][2] = Comp3DPosd( x - 1, y + 1 );
	p[1][2] = Comp3DPosd( x, y + 1 );
	p[2][2] = Comp3DPosd( x + 1, y + 1 );

	// compute du, dv, duu, duv, dvv from 9 pointspoints
	du.x = 0.5 * ( p[2][1].x - p[0][1].x );
	du.y = 0.5 * ( p[2][1].y - p[0][1].y );
	du.z = 0.5 * ( p[2][1].z - p[0][1].z );
	dv.x = 0.5 * ( p[1][2].x - p[1][0].x );
	dv.y = 0.5 * ( p[1][2].y - p[1][0].y );
	dv.z = 0.5 * ( p[1][2].z - p[1][0].z );
	duu.x = p[2][1].x + p[0][1].x - p[1][1].x - p[1][1].x;
	duu.y = p[2][1].y + p[0][1].y - p[1][1].y - p[1][1].y;
	duu.z = p[2][1].z + p[0][1].z - p[1][1].z - p[1][1].z;
	dvv.x = p[1][2].x + p[1][0].x - p[1][1].x - p[1][1].x;
	dvv.y = p[1][2].y + p[1][0].y - p[1][1].y - p[1][1].y;
	dvv.z = p[1][2].z + p[1][0].z - p[1][1].z - p[1][1].z;
	duv.x = 0.25 * ( p[2][2].x + p[0][0].x - p[2][0].x - p[0][2].x );
	duv.y = 0.25 * ( p[2][2].y + p[0][0].y - p[2][0].y - p[0][2].y );
	duv.z = 0.25 * ( p[2][2].z + p[0][0].z - p[2][0].z - p[0][2].z );

	// compute n = u x v
	Crossd( du, dv, &n );
	
	// second fundamental form 
	L = Dotd( duu, n );
	M = Dotd( duv, n );
	N = Dotd( dvv, n );

	return L * N - SQR( M );
}

double DepthImage::MeanCurvature( int x, int y ){
	Point3d p[3][3], n, du, dv, duu, dvv;
	double L, N;

	// out of range
	if( x <= 0 || y <= 0 || x > this->ImageR< int >::_w - 2 || y > this->ImageR< int >::_h - 2 ) return 0.0;
	
	// get 3D coordinate around ( x, y )
	p[0][0] = Comp3DPosd( x - 1, y - 1 );
	p[1][0] = Comp3DPosd( x, y - 1 );
	p[2][0] = Comp3DPosd( x + 1, y - 1 );
	p[0][1] = Comp3DPosd( x - 1, y );
	p[1][1] = Comp3DPosd( x, y );
	p[2][1] = Comp3DPosd( x + 1, y );
	p[0][2] = Comp3DPosd( x - 1, y + 1 );
	p[1][2] = Comp3DPosd( x, y + 1 );
	p[2][2] = Comp3DPosd( x + 1, y + 1 );

	// compute du, dv, duu, duv, dvv from 9 pointspoints
	du.x = 0.5 * ( p[2][1].x - p[0][1].x );
	du.y = 0.5 * ( p[2][1].y - p[0][1].y );
	du.z = 0.5 * ( p[2][1].z - p[0][1].z );
	dv.x = 0.5 * ( p[1][2].x - p[1][0].x );
	dv.y = 0.5 * ( p[1][2].y - p[1][0].y );
	dv.z = 0.5 * ( p[1][2].z - p[1][0].z );
	duu.x = p[2][1].x + p[0][1].x - p[1][1].x - p[1][1].x;
	duu.y = p[2][1].y + p[0][1].y - p[1][1].y - p[1][1].y;
	duu.z = p[2][1].z + p[0][1].z - p[1][1].z - p[1][1].z;
	dvv.x = p[1][2].x + p[1][0].x - p[1][1].x - p[1][1].x;
	dvv.y = p[1][2].y + p[1][0].y - p[1][1].y - p[1][1].y;
	dvv.z = p[1][2].z + p[1][0].z - p[1][1].z - p[1][1].z;

	// compute n = u x v
	Crossd( du, dv, &n );
	
	// second fundamental form 
	L = Dotd( duu, n );
	N = Dotd( dvv, n );

	return L + N;
}

int DepthImage::CompRotateTranslate( Point2i_Array pa, double Rt[16] ){
	Point2i *p2d;
	Point3d_Array pa3d;
	const int w = this->ImageR< int >::_w, h = this->ImageR< int >::_h;
	double m[4][4];
	int i, j;

	// initialize
	i = pa.n * sizeof( Point3d );
	if( ( int )_m_size < i ){
		_m_size = ( unsigned char )i;
		_mem = realloc( _mem, _m_size );
		if( _mem == NULL ) return 0;
	}
	pa3d.pt = ( Point3d * )_mem;

	p2d = pa.pt;
	for( i = j = 0; i < pa.n; i++, p2d++ ){
		if( _bits[ ( *p2d ).y * _stride + ( *p2d ).x ] == 0 ) continue;
		pa3d.pt[ j++ ] = Comp3DPosd( ( *p2d ).x, ( *p2d ).y );
	}
	pa3d.n = j;
	CompPCACoordd( pa3d, m );
	memcpy( Rt, m, 16 * sizeof( double ) );
	return 1;
}

//void 