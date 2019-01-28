//
// measure.cpp: about measureimg foot anthropometry data
//

#include "tracking.h"
#include "config.h"
#include <image\imgconf.h>
#include <image\imgproc.h>
#include <lib3D\depthimage.h>
#include <lib3D\3dconfig.h>
#include <float.h>
#include <gl\gl.h>

using namespace FootTracker;

//
// Global variables
//

extern DepthImage gdimg[2];
extern ImageR< int > gtexture[2];
extern ImageR< unsigned int > *gbg;

//
// Global Memory
//

extern Memory gmem;

// 
// inline functions
//

inline static int DepthImageSegment( DepthImage &img, ImagePoint p, Region region ){
	static ImageBuffer imgbuf;
	Color1i zero = { 0 };
	int i, j;

	// setting image buffer
	imgbuf.img = img.GetImage();
	imgbuf.w = region.w = img.Width();
	imgbuf.h = region.h = img.Height();
	imgbuf.stride = img.Stride();
	imgbuf.color = 1;

	// select feet region
	if( imgbuf.img == NULL ) return 0;
	Selecti( imgbuf, region, p, 100 );

	// inverse selection
	for( i = 0, j = imgbuf.w * imgbuf.h; i < j; i++ ) region.img[i] = !region.img[i];
	Filling1i( imgbuf, region, zero );
	return 1;
}

float Tracker::CompFootLength( int color, int threshold ){
	const float r = 150.0f;
	const int 
		w = gdimg->Width(),
		h = gdimg->Height();
	Point3f_Array pts;
	Point3f p3d, pmax, pmin;
	ImagePoint p;
	Region buf;
	double param[9];
	float min = FLT_MAX, max = -FLT_MAX;
	int i;

	// intrinsic camera matrix
	gdimg->GetIntrinsic( param );

	if( ( pts.n = _CountColorPoints( color, threshold ) ) == 0 ) 		// counting points
		return 0.0f;
	Allocate( pts.n * sizeof( Point3f ) );				// allocate
	Assign( 1, &pts.pt, pts.n * sizeof( Point3f ) );		// assign to points array
	_CompColorPoints( color, threshold, pts.pt );		// get data

	// compute mean of color points at pts.pt[0]
	glBegin( GL_POINTS );
	glColor3f( 0.0, 1.0, 0.0 );
	for( i = 1; i < pts.n; i++ ){
		pts.pt[0].x += pts.pt[i].x;
		pts.pt[0].y += pts.pt[i].y;
		pts.pt[0].z += pts.pt[i].z;
		glVertex3f( pts.pt[i].x, pts.pt[i].y, pts.pt[i].z );
	}
	glEnd();
	p3d.x = pts.pt[0].x / ( float )pts.n;	
	p3d.y = pts.pt[0].y / ( float )pts.n;	
	p3d.z = pts.pt[0].z / ( float )pts.n;
	// project to 2D image
	p.x = ( int )( ( param[0] * pts.pt[0].x + param[1] * pts.pt[0].y - param[2] * pts.pt[0].z ) / -pts.pt[0].z );
	p.y = ( int )( ( param[4] * pts.pt[0].y - param[5] * pts.pt[0].z ) / -pts.pt[0].z );

	// segmentation
	Allocate( w * h );
	Assign( 1, &buf.img, w * h );
	DepthImageSegment( gdimg[0], p, buf );

	// Select points in radius
	if( ( pts.n = _CountRadiusPoints( p3d.x, p3d.y, p3d.z, r, 1 ) ) == 0 ) // counting points
		return 0.0f;
	Allocate( pts.n * sizeof( Point3f ) );					// allocate
	Assign( 1, &pts.pt, pts.n * sizeof( Point3f ) );		// assign to points array
	_CompRadiusPoints( p3d.x, p3d.y, p3d.z, r, pts.pt, pts.n, 1 );	// compute 3D points

	// find minimum and maximum x values
	for( i = 0; i < pts.n; i++ ){
		if( pts.pt[i].x < min ) {
			min = pts.pt[i].x;
			pmin = pts.pt[i];
		}
		if( pts.pt[i].x > max ) {
			max = pts.pt[i].x;
			pmax = pts.pt[i];
		}
	}
	glBegin( GL_LINES );
	glColor3f( 0.0f, 1.0f, 0.0f );
	glVertex3fv( ( float * )&pmin );
	glVertex3fv( ( float * )&pmax );
	glEnd();
	return max - min;
}

float Tracker::CompFootWidth( int color, int threshold ){
	const float r = 90.0f;
	const int 
		w = gdimg->Width(),
		h = gdimg->Height();
	Point3f_Array pts;
	Point3f p3d, pmax, pmin;
	ImagePoint p;
	Region buf;
	double param[9];
	float min = FLT_MAX, max = -FLT_MAX;
	int i;

	// intrinsic camera matrix
	gdimg->GetIntrinsic( param );

	if( ( pts.n = _CountColorPoints( color, threshold ) ) == 0 ) 		// counting points
		return 0.0f;
	Allocate( pts.n * sizeof( Point3f ) );				// allocate
	Assign( 1, &pts.pt, pts.n * sizeof( Point3f ) );		// assign to points array
	_CompColorPoints( color, threshold, pts.pt );		// get data

	// compute mean of color points at pts.pt[0]
	glBegin( GL_POINTS );
	glColor3f( 0.0, 1.0, 0.0 );
	for( i = 1; i < pts.n; i++ ){
		pts.pt[0].x += pts.pt[i].x;
		pts.pt[0].y += pts.pt[i].y;
		pts.pt[0].z += pts.pt[i].z;
		glVertex3f( pts.pt[i].x, pts.pt[i].y, pts.pt[i].z );
	}
	glEnd();
	p3d.x = pts.pt[0].x / ( float )pts.n;	
	p3d.y = pts.pt[0].y / ( float )pts.n;	
	p3d.z = pts.pt[0].z / ( float )pts.n;
	// project to 2D image
	p.x = ( int )( ( param[0] * pts.pt[0].x + param[1] * pts.pt[0].y - param[2] * pts.pt[0].z ) / -pts.pt[0].z );
	p.y = ( int )( ( param[4] * pts.pt[0].y - param[5] * pts.pt[0].z ) / -pts.pt[0].z );

	// segmentation
	Allocate( w * h );
	Assign( 1, &buf.img, w * h );
	DepthImageSegment( gdimg[0], p, buf );

	// Select points in radius
	if( ( pts.n = _CountRadiusPoints( p3d.x, p3d.y, p3d.z, r, 1 ) ) == 0 ) // counting points
		return 0.0f;
	Allocate( pts.n * sizeof( Point3f ) );					// allocate
	Assign( 1, &pts.pt, pts.n * sizeof( Point3f ) );		// assign to points array
	_CompRadiusPoints( p3d.x, p3d.y, p3d.z, r, pts.pt, pts.n, 1 );	// compute 3D points

	// find minimum and maximum x values
	for( i = 0; i < pts.n; i++ ){
		if( pts.pt[i].x < min ) {
			min = pts.pt[i].x;
			pmin = pts.pt[i];
		}
		if( pts.pt[i].x > max ) {
			max = pts.pt[i].x;
			pmax = pts.pt[i];
		}
	}
	glBegin( GL_LINES );
	glColor3f( 0.0f, 1.0f, 0.0f );
	glVertex3fv( ( float * )&pmin );
	glVertex3fv( ( float * )&pmax );
	glEnd();
	return max - min;
}
/**/