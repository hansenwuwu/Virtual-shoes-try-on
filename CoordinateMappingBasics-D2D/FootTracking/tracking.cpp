// tracking.cpp : 定義 DLL 應用程式的匯出函式。
//
#define _FOOT_TRACKING_CPP_

#include <iostream>
#include "../stdafx.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <fstream>
#include <windows.h>
//#include <GL/glut.h>
#include "../include/GL/glut.h"
//#include <GL/gl.h>
#include "../numerical/linear.h"
#include "../lib3D/depthimage.h"
#include "../lib3D/3dconfig.h"
#include "../lib3D/track.h"
#include "../lib3D/calibration.h"
#include "../image/vision.h"
#include "../image/imgproc.h"
#include "../image/imgproc2.h"
#include "../FootTracking/tracking.h"
#include "../FootTracking/config.h"
//#include <cv.h>
//#include <highgui.h>

using namespace std;
using namespace FootTracker;


//using namespace KinectScene;

//
// Macro
//

#define SQR(x)	( (x) * (x) )
#define M_PIf			3.1415926535897932384626433832795f
#define M_HalfPIf		1.5707963267948966192313216916398f
//
// Global variables
//

DepthImage gdimg[2];
ImageR< int > gtexture[2];
ImageR< unsigned int > *gbg = NULL;

//
// Global Memory
//

Memory gmem;

//
// CUDA
//
//
//
/*#pragma comment(lib,"CudaDll.lib") 
extern "C" float CudaICPf( Point3f_Array patch, Point3f_Array base, RT_Paramf *theta, const ICP_Paramf param, const Point3f_Array buf, float icpe[100] ); 
*/

//
// inline function
//

void getRotateAngle(double *x, double *y, double *angle){
	double dot, len1, len2;
	dot = x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
	len1 = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
	len2 = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
	*angle = acos(dot / (len1 * len2)) * 180 / M_PIf;
}

void initPosition(RT_Paramf *theta, int LR){
	if (LR == FT_LEFT){
		theta[0].x = 230;
		theta[0].y = -137;
		theta[0].z = -1100;
		theta[0].a = -0.345780;
		theta[0].b = -0.791593;
		theta[0].c = -1.168923;
	}
	else{
		theta[0].x = -12;
		theta[0].y = -137;
		theta[0].z = -1100;
		theta[0].a = 0.568555;
		theta[0].b = -0.743519;
		theta[0].c = -2.108131;
	}
}

inline int ThisPos( int color, int threshold, float &x, float &y, float &z ){
	// need to initialize
	if( gdimg->Empty() || gtexture->Empty() ) return 0;

	const int 
		w = gdimg->Width(), 
		h = gdimg->Height();
	ImageBuffer imgbuf, dbuf;
	Point3f buf;
	bool (*pIsColor)( int* );
	int *dbit, *gbit, i, j, n;

	// select color
	switch( color ){
		case FT_GREEN:
			pIsColor = IsGreen;
			break;
		case FT_SKIN:
			pIsColor = IsSkin;
			break;
		default:
			return 0;
	}

	// setting image buffer
	imgbuf.img = gtexture[0].GetImage();
	imgbuf.w = dbuf.w = w;
	imgbuf.h = dbuf.h = h;
	imgbuf.stride = gtexture[0].Stride();
	imgbuf.color = 3;
	dbuf.img = gdimg[0].GetImage();
	dbuf.stride = gdimg[0].Stride();
	dbuf.color = 1;

	// mean of green points
	x = y = z = 0.0f;
	n = 0;
	for( i = 0; i < h; i++ ){
		gbit = ( int* )imgbuf.img + i * imgbuf.stride;
		dbit = ( int* )dbuf.img + i * dbuf.stride;
		for( j = 0; j < w; j++, gbit += 3, dbit++ ){
			if( dbit[0] != 0 && ( pIsColor )( gbit ) ){	// check color
				buf = gdimg[0].Comp3DPosf( j, i );
				x += buf.x;
				y += buf.y;
				z += buf.z;
				n++;
			}
		}
	}
	x /= ( float )n;
	y /= ( float )n;
	z /= ( float )n;
	return n;
}

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
	if ( imgbuf.img == NULL ) return 0;
	//time_t tt,tt2,tt3;
	//tt = clock();
	Selecti( imgbuf, region, p, 500 ); //
	//printf("tt: %d\n", clock() - tt);
	
	// inverse selection
	//tt2 = clock();
	for( i = 0, j = imgbuf.w * imgbuf.h; i < j; i++ ) region.img[i] = !region.img[i];
	//printf("tt2: %d\n", clock() - tt2);
	//tt3 = clock();
	Filling1i( imgbuf, region, zero );
	//printf("tt3: %d\n", clock() - tt3);
	//-----測試切割結果-----
	int *temp = (int *)imgbuf.img;
	int count = 0;
	int col, row;
	for (size_t i = 0; i < 1920 * 1080; i++){
		row = 1080 - 1 - i / 1920;
		col = i % 1920;
		//color img
		if (temp[i] != 0){
			count++;
		}
	}
	printf("count = %d\n", count);
//	int *temp = (int *)imgbuf.img;
	color3uc *img1 = (color3uc *)malloc(1920 * 1080 * sizeof(color3uc));
//	int col, row;
	for (size_t i = 0; i < 1920 * 1080; i++){
		row = 1080 - 1 - i / 1920;
		col = i % 1920;
		//color img
		if (temp[i] == 0){
			img1[row * 1920 + col].r = 255;
			img1[row * 1920 + col].g = 255;
			img1[row * 1920 + col].b = 255;
		}
		else{
			img1[row * 1920 + col].r = 0;
			img1[row * 1920 + col].g = 0;
			img1[row * 1920 + col].b = 0;
		}		
	}
	glDisable(GL_DEPTH_TEST);
	glPixelZoom((float)1, (float)1);
	//glRasterPos2f( -1.f, 1.f );
	glPixelStorei(GL_PACK_ALIGNMENT, 4);
	glDrawPixels(1920, 1080, GL_RGB, GL_UNSIGNED_BYTE, img1);
	glEnable(GL_DEPTH_TEST);
	
	return 1;
}

//
// Function
//

// Setting

int   FootTracker::InitTracker( point3f* _depth, int w, int h ){
	double *a, *pa1, *pa2, *b, *pb1, *pb2, temp;
	point3f *scp = NULL, *pscp;
	int i, j, k, n;

	// construct classes
	if( gdimg->Empty() ){
		gdimg[0].Initialize( w, h, 1 );
		gdimg[1].Initialize( w, h, 1 );
	}
	if( gtexture->Empty() ){
		gtexture[0].Initialize( w, h, 3 );
		gtexture[1].Initialize( w, h, 3 );
	}
	if( gbg == NULL ){
		// 3-dimension: mean, variance, n
		gbg = new ImageR< unsigned int >;
		gbg->Initialize( w, h, 3 );
		// initialize to zero
		memset( gbg->GetImage(), 0x00, h * gbg->Stride() * sizeof( int ) );
	}

	// memory pool
	gmem.Allocate( 12 * w * h * sizeof( double ) + 
		6 * 6 * sizeof( double ) + 
		6 * ( sizeof( int ) + sizeof( double ) ) +
		w * h * sizeof( point3f ) );
	if( gmem.Size() == 0 ) return 0;
	a = ( double* )gmem.Offset( 6 * ( sizeof( double ) + sizeof( int ) ) );	// w * h * 12
	b = a + w * h * 12;														// 6 * 6
	scp = ( point3f * )( b + 36 );

//	ni.GetDepth( scp );
	for( i = 0; i < h; i++ )
	{
		for( j = 0; j < w; j++ )
		{			
			(*scp) = point3f( (float)(_depth)->x, (float)(_depth)->y , (float)(_depth)->z );
			scp++;
			_depth++;
		}
	}
	scp -= (w*h);

	// construct matrix A
	n = 0;								// count the meaningful data
	pa1 = a;	
	pa2 = a + 6;	
	pscp = scp;
	for( i = 0; i < h; i++ ){
		for( j = 0; j < w; j++, pscp++ ){
			if( fabs( pscp[0].z ) < FLT_MIN ) continue;
			// row 1
			pa1[0] = -pscp[0].x;//!!!!!!!!!!!!!!!!1
			pa1[1] = pscp[0].y;
			pa1[2] = pscp[0].z;
			pa1[3] = pa1[4] = 0.0;
			pa1[5] = -pscp[0].z * ( float )j;
			// row 2
			pa2[0] = pa2[1] = pa2[2] = 0.0;
			pa2[3] = pscp[0].y;
			pa2[4] = pscp[0].z;
			pa2[5] = -pscp[0].z * ( float )i;
			pa1 += 12;	pa2 += 12;	n += 2;
		}
	}

	// B = AT * A
	for( i = 0; i < 6; i++ ){
		pb1 = pb2 = b + i * 6 + i;
		for( j = i; j < 6; j++, pb1++, pb2 += 6 ){
			pa1 = a + i;
			pa2 = a + j;
			temp = 0.0;
			for( k = 0; k < n; k++, pa1 += 6, pa2 += 6 ) temp += pa1[0] * pa2[0];
			pb1[0] = pb2[0] = temp;
		}
	}

	// Find Minimal eigenvector
	InvIter( b, a, 0.0, 6, 30, ( double* )gmem.Offset(0) );
	// Set projection
	for( i = 0; i < 5; i++ ) a[i] /= a[5];
	gdimg[0].SetIntrinsic( a[0], a[3], a[1], a[2], a[4] );
	gdimg[1].SetIntrinsic( a[0], a[3], a[1], a[2], a[4] );
	gdimg[0].SetSize( h, w );
	gdimg[1].SetSize( h, w );

	// free memory
	return 1;
}


int   FootTracker::ExitTracker( void ){
	if( gdimg->Empty() == 0 ) {
		gdimg[0].Clear();
		gdimg[1].Clear();
	}
	if( gtexture->Empty() == 0 ) {
		gtexture[0].Clear();
		gtexture[1].Clear();
	}
	if( gbg ){
		delete gbg;
		gbg = NULL;
	}
	gmem.Clear();
	return 1;
}

int   FootTracker::SetTrackerDepth( float *depth_img ){
	if( gdimg->Empty() ) return 0;
	
	float *depth_bit;
	const int 
		w = gdimg->Width(), 
		h = gdimg->Height(),
		stride = gdimg->Stride();
	unsigned *gbg_bit;
	int *gbit, i, j;

	// save privious image
	memcpy( gdimg[1].GetImage(), gdimg[0].GetImage(), gdimg[0].Stride() * gdimg[0].Height() * sizeof( float ) );

	// Setting Depth ImageR
	if( gbg->Color() != 1 ){			// still training
		gbit = gdimg->GetImage();
		for( i = 0; i < h; i++ ){
			depth_bit = depth_img + ( h - i - 1 ) * w;
			//depth_bit = depth_img ;
			for( j = 0; j < w; j++, depth_bit++, gbit++ ) {
				// filter the null point and point that is too far
				gbit[0] = ( depth_bit[0] < 500 || depth_bit[0] > 1500 ) ? 0 : ( int )depth_bit[0]; //50 3000
			}
		}		
	}
	else{
		gbg_bit = gbg->GetImage();
		gbit = gdimg->GetImage();
		for( i = 0; i < h; i++ ){
			//depth_bit = depth_img + ( h - i - 1 ) * w;
			depth_bit = depth_img;
			for( j = 0; j < w; j++, depth_bit++, gbit++, gbg_bit++ ) {
				// filter the null point and point that is too far
				gbit[0] = ( depth_bit[0] < 500 || depth_bit[0] > 1500 || depth_bit[0] > gbg_bit[0] ) ?  //50 3000
					0 : ( int )depth_bit[0]; 
			}
		}

	}
	return 1;
}

int   FootTracker::SetTrackerTexture( color3uc *img ){
	if( gtexture->Empty() ) return 0;
	color3uc *imgbit;
	const int 
		w = gtexture->Width(), 
		h = gtexture->Height(),
		stride = gtexture->Stride();
	int *gbit, *gbits, i, j;

	// save privious image
	memcpy( gtexture[1].GetImage(), gtexture[0].GetImage(), stride * h * sizeof( int ) );

	// Setting Texture ImageR
	gbits = ( int* )gtexture->GetImage();
	for( i = 0; i < h; i++ ){
		//imgbit = img + ( h - i - 1 ) * w;
		imgbit = img;
		gbit = gbits + i * stride;
		for( j = 0; j < w; j++, imgbit++, gbit += 3 ){
			// from RGB to NTSC color space
			gbit[0] = ( 299 * ( int )imgbit[0].r + 587 * ( int )imgbit[0].g + 114 * ( int )imgbit[0].b ) / 1000;
			gbit[1] = ( 596 * ( int )imgbit[0].r - 274 * ( int )imgbit[0].g - 322 * ( int )imgbit[0].b ) / 1000;
			gbit[2] = ( 211 * ( int )imgbit[0].r - 523 * ( int )imgbit[0].g + 312 * ( int )imgbit[0].b ) / 1000;
		}
	}
	return 1;
}

int   FootTracker::SetTrackerBG( XnDepthPixel *depth_img ){
	if( gbg == NULL || gbg->Empty() || gbg->Color() != 3 ) return 0;

	XnDepthPixel *depth_bit;
	const int 
		w = gbg->Width(), 
		h = gbg->Height(),
		stride = gbg->Stride();
	unsigned int *gbit;
	int i, j;
	
	// Setting Depth ImageR
	for( i = 0; i < h; i++ ){
		gbit = gbg->GetImage() + i * stride;
		depth_bit = depth_img + ( h - i - 1 ) * w;
		for( j = 0; j < w; j++, depth_bit++, gbit += 3 ) {
			if( depth_bit[0] < 50 ) continue;
			gbit[0] += ( unsigned int )depth_bit[0];
			gbit[1] += ( unsigned int )SQR( depth_bit[0] );
			gbit[2] ++;
		}
	}
	return 1;
}

int   FootTracker::TrackerBGTraining( int p ){
	const int 
		w = gbg->Width(), 
		h = gbg->Height(), 
		c = gbg->Color();
	unsigned int *bit, *bgbits;
	float ni, u, v;
	int stride, i, j;
	switch( p ){
		case FT_BEGIN:
			gbg->SetColor( 3 );
			memset( gbg->GetImage(), 0x00, h * gbg->Stride() * sizeof( int ) );
			break;
		case FT_END:
			// FT_END without FT_BEGIN first
			if( gbg->Color() != 3 ) return 0;
			
			bgbits = gbg->GetImage();
			stride = gbg->Stride();
			for( i = 0; i < h; i++ ){
				bit = bgbits + i * stride;
				for( j = 0; j < w; j++, bit += c ){
					// case of never detect depth image
					if( bit[2] == 0 ){
						memset( bit, 0xff, 3 * sizeof( int ) );	// set to maximum value
						continue;
					}
					// V[x] = E[x^2] - [E(x)]2
					ni = 1.0f / ( float )bit[2];
					u = ( float )bit[0] * ni;
					v = ( float )bit[1] * ni - SQR( u );
					// set to u - 3v
					bit[0] = bit[1] = bit[2] = ( unsigned int )( u - 0.5 * u );
				}
			}
			gbg->SetColor( 1 );
			break;
	}
	return p;
}

//
// Tracker
//

//
// Protected function
//

int FootTracker::Tracker::_CountColorPoints( int &color , int &threshold ){
	const int 
		w = gdimg->Width(), 
		h = gdimg->Height();
	Point3f_Array p;
	ImageBuffer imgbuf, dbuf;
	bool (*pIsColor)( int* );
	int *dbit, *gbit, i, j;

	// select color
	switch( color ){
		case FT_GREEN:
			pIsColor = IsGreen;
			break;
		case FT_SKIN:
			pIsColor = IsSkin;
			break;
		default:
			return 0;
	}

	// setting image buffer
	imgbuf.img = gtexture->GetImage();
	imgbuf.w = dbuf.w = w;
	imgbuf.h = dbuf.h = h;
	imgbuf.stride = gtexture->Stride();
	imgbuf.color = 3;
	dbuf.img = gdimg->GetImage();
	dbuf.stride = gdimg->Stride();
	dbuf.color = 1;

	// counting green points
	p.n = 0;
	p.pt = NULL;
	for( i = 0; i < h; i++ ){
		gbit = ( int* )imgbuf.img + i * imgbuf.stride;
		dbit = ( int* )dbuf.img + i * dbuf.stride;
		for( j = 0; j < w; j++, gbit += 3, dbit++ ){
			if( dbit[0] != 0 && ( pIsColor )( gbit ) )	// check color
				p.n++;
		}
	}
	// if not enough point
	return ( p.n < threshold ) ? 0 : p.n;
}
	
int FootTracker::Tracker::_CompColorPoints( int &color, int &threshold, void *ptr ){
	const int 
		w = gdimg->Width(), 
		h = gdimg->Height();
	Point3f_Array p = { ( Point3f * )ptr, 0 };
	ImageBuffer imgbuf, dbuf;
	bool (*pIsColor)( int* );
	int *dbit, *gbit, i, j;

	// select color
	switch( color ){
		case FT_GREEN:
			pIsColor = IsGreen;
			break;
		case FT_SKIN:
			pIsColor = IsSkin;
			break;
		default:
			return 0;
	}

	// setting image buffer
	imgbuf.img = gtexture->GetImage();
	imgbuf.w = dbuf.w = w;
	imgbuf.h = dbuf.h = h;
	imgbuf.stride = gtexture->Stride();
	imgbuf.color = 3;
	dbuf.img = gdimg->GetImage();
	dbuf.stride = gdimg->Stride();
	dbuf.color = 1;

	// record points
	for( i = 0; i < h; i++ ){
		gbit = ( int* )imgbuf.img + i * imgbuf.stride;
		dbit = ( int* )dbuf.img + i * dbuf.stride;
		for( j = 0; j < w; j++, gbit += 3, dbit++ ){
			if( dbit[0] != 0 && ( pIsColor )( gbit ) ){	// check color
				p.pt[ p.n++ ] = gdimg->Comp3DPosf( j, i );
			}
		}
	}
	return p.n;
}

int FootTracker::Tracker::_CountRadiusPoints( float &x, float &y, float &z, float r, int jmp ){
	const int 
		w = gdimg->Width(), 
		h = gdimg->Height();
	const float r2 = SQR( r );
	Point3f temp, pt;
	int i, j, n;
	
	// counting the point number of foot's depth image
	n = 0;
	
	for( i = 0; i < h; i += jmp ){
		for( j = 0; j < w; j += jmp ){
			// compute distance to center
			pt = gdimg->Comp3DPosf( j, i );
			temp.x = pt.x - x;
			temp.y = pt.y - y;
			temp.z = pt.z - z;
			if (SQR(temp.x) + SQR(temp.y) + SQR(temp.z) > r2) continue;	// counting
			n++;
		}
	}
	return n;

}

int FootTracker::Tracker::_CompRadiusPoints( float &x, float &y, float &z, float r, void *ptr, int &n, int jmp ){
	const int 
		w = gdimg->Width(), 
		h = gdimg->Height();
	const float r2 = SQR( r );
	Point3f_Array depth_array = { ( Point3f * )ptr, n };
	Point3f temp, pt;
	int i, j, k;

	glBegin(GL_POINTS);
	glColor3f((float)1.0, (float)1.0, (float)1.0);
	// record the point of foot's depth image
	for( i = k = 0; i < h; i += jmp ){
		for( j = 0; j < w; j += jmp ){
			// compute distance to center
			pt = gdimg->Comp3DPosf( j, i );
			temp.x = pt.x - x;
			temp.y = pt.y - y;
			temp.z = pt.z - z;
			if( SQR( temp.x ) + SQR( temp.y ) + SQR( temp.z ) > r2 ) continue;
			depth_array.pt[ k++ ] = pt; //如果該點與中心點的距離小於r2，則把該點加入depth_array(足部分割後的深度點群)
			if( k == depth_array.n ) break;
		}
		if( j < w ) break;
	}
	glEnd();	
	return k;
}


void FootTracker::Tracker::_DepthImgSegment( float &x, float &y, float &z ){
	double param[9];
	Region buf;
	ImagePoint p;

	// intrinsic matrix
	gdimg->GetIntrinsic( param );
	// memory pool
	buf.img = ( char* )this->Allocate( gdimg->Width() * gdimg->Height() );
	// image point
	p.x = ( int )( ( param[0] * x + param[1] * y - param[2] * z ) / -z );
	p.y = ( int )( ( param[4] * y - param[5] * z ) / -z );
	// segmentation
	DepthImageSegment( gdimg[0], p, buf );	
}
//
// Public functions
//

  FootTracker::Tracker::Tracker( void ): Memory(){
	memset( _prev, 0x00, sizeof( _prev ) );
}

int   FootTracker::Tracker::TrackColor( int color, int threshold, float m[16] ){
	// need to initialize
	if( gdimg->Empty() || gtexture->Empty() ) return 0;
	Point3f_Array p = { NULL, 0 };
	float temp[4][4];
	int i;

	// count
	if( ( p.n = _CountColorPoints( color, threshold ) ) == 0 ) return 0;

	// memory allocate
	p.pt = ( Point3f * )Allocate( p.n * sizeof( Point3f ) );
	if( p.pt == NULL ) return 0;
	// record
	_CompColorPoints( color, threshold, p.pt );
	// draw white points
	
	glBegin( GL_POINTS );
	glColor3f( 1.0, 0.0, 0.0 );
	
	for( i = 0; i < p.n; i++ ) glVertex3f( p.pt[i].x, p.pt[i].y, p.pt[i].z );
	glEnd();
	//
	glBegin( GL_LINES );
		glColor3f( 1.0, 0.0, 0.0 );
		glVertex3f( 0.0, 0.0, 0.0 );	glVertex3f( 100.0, 0.0, 0.0 );
		glColor3f( 0.0, 1.0, 0.0 );
		glVertex3f( 0.0, 0.0, 0.0 );	glVertex3f( 0.0, 100.0, 0.0 );
		glColor3f( 0.0, 0.0, 1.0 );
		glVertex3f( 0.0, 0.0, 0.0 );	glVertex3f( 0.0, 0.0, 100.0 );
	glEnd();
	//
	CompPCACoordf( p, temp );

	if( temp[2][2] < 0 ){	// Zz must > 0
		temp[0][1] *= -1.0;	temp[1][1] *= -1.0;	temp[2][1] *= -1.0;
		temp[0][2] *= -1.0;	temp[1][2] *= -1.0;	temp[2][2] *= -1.0;
	}
	if( temp[1][0] > 0 ){		// Xy must < 0
		temp[0][0] *= -1.0;	temp[1][0] *= -1.0;	temp[2][0] *= -1.0;
		temp[0][1] *= -1.0;	temp[1][1] *= -1.0;	temp[2][1] *= -1.0;
	}
	memcpy( m, temp[0], 16 * sizeof( float ) );

	// save transformation
	memcpy( _prev, m, 12 * sizeof( float ) );

	return p.n;
}
int	  FootTracker::Tracker::TrackInitSolution(float m[16], point3f center)
{
	m[0] = -0.239761; m[1] = 0.968747; m[2] = -0.063598;
	m[4] = -0.672416; m[5] = -0.118455; m[6] = 0.730633;
	m[8] = 0.700265; m[9] = 0.217942; m[10] = 0.679802;
	m[3] = center.x;
	m[7] = center.y;
	m[11] = center.z;
	memcpy(_prev, m, 12 * sizeof(float));
	return 0;
}

int   FootTracker::Tracker::TrackSelect( Point3f_Array spoints,float m[16] ){
	// need to initialize
	if( gdimg->Empty() || gtexture->Empty() ) return 0;

	Point3f_Array p = { NULL, 0 };
	float temp[4][4];

	CompPCACoordf( spoints, temp );

	if( temp[2][2] < 0 ){	// Zz must > 0
		temp[0][1] *= -1.0;	temp[1][1] *= -1.0;	temp[2][1] *= -1.0;
		temp[0][2] *= -1.0;	temp[1][2] *= -1.0;	temp[2][2] *= -1.0;
	}
	if( temp[1][0] > 0 ){		// Xy must < 0
		temp[0][0] *= -1.0;	temp[1][0] *= -1.0;	temp[2][0] *= -1.0;
		temp[0][1] *= -1.0;	temp[1][1] *= -1.0;	temp[2][1] *= -1.0;
	}
	memcpy( m, temp[0], 16 * sizeof( float ) );

	// save transformation
	memcpy( _prev, m, 12 * sizeof( float ) );
	return p.n;
}

float FootTracker::Tracker::TrackICP( MeshBuffer mesh, float m[16], float icpe[100], int LR ){
	if( gdimg->Empty() ) return 0;
	//time_t t,t2,t_total;
	Point3f_Array mesh_array, depth_array, buf_array;
	Point3f center;
	float *ptr, m1[3][4], e1;
	bool *v_list;
	RT_Paramf theta[1];
	ICP_Paramf param;
	int i;
	//t_total = clock();
	// initialize
	mesh_array.n = mesh.n;
	depth_array.n = 0;
	center.x = m[3] + 50.0f * m[0];
	center.y = m[7] + 50.0f * m[4];
	center.z = m[11] + 50.0f * m[8];
	memset( theta, 0x00, 1 * sizeof( RT_Paramf ) );

	// depth image segmentation 用顏色辨識取代
	//t = clock();
//	if( m[11] < -50 ) 
//		_DepthImgSegment( center.x, center.y, center.z ); // 深度切割 傳入中心點	
	// count points in radius
	depth_array.n = _CountRadiusPoints(center.x, center.y, center.z, 170.0f);
	if (depth_array.n == 0 || depth_array.n > 10000) {
		/*_e = FLT_MAX;
		return 0;*/
		// initial position
		initPosition(&theta[0], LR);
	}else{
		// memory allocate from memory pool
		mesh_array.pt = (Point3f*)Allocate(2 * (mesh_array.n + depth_array.n) * sizeof(Point3f));
		if (mesh_array.pt == NULL) {									// mesh_array.n
			_e = FLT_MAX;
			return 0;
		}
		depth_array.pt = mesh_array.pt + mesh_array.n;					// depth_array.n
		buf_array.pt = depth_array.pt + depth_array.n;					// mesh_array.n
		// record points
		_CompRadiusPoints(center.x, center.y, center.z, 170.0f, depth_array.pt, depth_array.n);
		//printf("_DepthImgSegment: %d\n",  clock() - t );
		// 用顏色辨識取代到此

		// trim the vertex that can't see 
		v_list = (bool*)buf_array.pt;
		memset(v_list, 0x00, mesh.n * sizeof(bool));
		// record the points can be see
		for (i = 0; i < mesh.nf; i++){
			ptr = mesh.vn + 3 * mesh.f[i].normal[0];	// normal of first point
			if (m[8] * ptr[0] + m[9] * ptr[1] + m[10] * ptr[2] > -0.f) v_list[mesh.f[i].vertex[0]] = true;
			ptr = mesh.vn + 3 * mesh.f[i].normal[1];	// normal of second point
			if (m[8] * ptr[0] + m[9] * ptr[1] + m[10] * ptr[2] > -0.f) v_list[mesh.f[i].vertex[1]] = true;
			ptr = mesh.vn + 3 * mesh.f[i].normal[2];	// normal of third point
			if (m[8] * ptr[0] + m[9] * ptr[1] + m[10] * ptr[2] > -0.f) v_list[mesh.f[i].vertex[2]] = true;
		}

		// save the points
		ptr = mesh.v;
		for (i = mesh_array.n = 0; i < mesh.n; i++, ptr += 3){
			if (v_list[i] == false) continue;
			memcpy(mesh_array.pt + mesh_array.n, ptr, 3 * sizeof(float));
			mesh_array.n++;
		}
		// compute transform matrix
		//(5,10)
		param.iter = 50;
		param.trans.iter = 20;
		param.trans.err = 0.0001f;
		param.trans.lr1 = 0.05f;
		param.trans.lr2 = 0.4f;
		memcpy(m1[0], m, 12 * sizeof(float));
		CompRTFromMatf(theta, m1);
		//t = clock();
		//e1 = ICPf(mesh_array, depth_array, theta, param, buf_array, icpe);
		e1 = ICP_kdtree(mesh_array, depth_array, theta, param, buf_array, icpe);
	}
	
	//printf("ICP: %d\n", clock() - t);
	//if (theta[0].a>M_HalfPIf || theta[0].a<-M_HalfPIf)//abs(theta[0].a) < 1.5 && abs(theta[0].b) < 1.5 && abs(theta[0].c) < 1.5)
	//{
	//	//(0, 0, -pi/2) 鞋尖向下
	//	//(-pi,-pi,0) 鞋尖向前
	//	theta[0].a = -M_HalfPIf;
	//	theta[0].b = -M_HalfPIf; //繞上方旋轉,左右轉
	//	theta[0].c = 0.0;
	//}
	//if (1)
	//{
	//	theta[0].a = -M_HalfPIf;
	//	theta[0].b = -M_HalfPIf;
	//	theta[0].c = 0;
	//}
	/*if (theta[0].a>M_HalfPIf || theta[0].a<-M_HalfPIf)
	{
		theta[0].a = -M_HalfPIf;		
	}*/
	/*if (theta[0].b>0.0 || theta[0].b<-M_PIf)
	{
		theta[0].b = -M_HalfPIf;
	}*/
	/*if (theta[0].a>0.0 || theta[0].a<-M_PIf)
	{
		if (theta[0].c<M_HalfPIf && theta[0].c>-M_HalfPIf)
			theta[0].a = -M_HalfPIf;
	}
	if (theta[0].c>M_HalfPIf || theta[0].c<-M_HalfPIf)
	{
		if(theta[0].a<0.0 && theta[0].a>-M_PIf)
			theta[0].c = 0;
	}*/
	CompRTMatf( m1, theta[0] );
	//printf("(%f, %f, %f) (%f, %f, %f)\n", theta[0].x, theta[0].y, theta[0].z, theta[0].a, theta[0].b, theta[0].c);

	// check the invalid situations
	double nv[3], nv1[3], rotateAngle, rotateAngle1, rotateAngle2, ny[3], nz[3];
	nv[0] = m1[0][2];
	nv[1] = m1[1][2];
	nv[2] = m1[2][2];
	nv1[0] = m1[0][0];
	nv1[1] = m1[1][0];
	nv1[2] = m1[2][0];
	ny[0] = 0;
	ny[1] = 1;
	ny[2] = 0;
	nz[0] = 0;
	nz[1] = 0;
	nz[2] = 1;

	getRotateAngle(nv, ny, &rotateAngle);
	getRotateAngle(nv1, nz, &rotateAngle1);
	getRotateAngle(nv1, ny, &rotateAngle2);
	//printf("nv = (%lf, %lf, %lf), rotateAngle = %lf\n", nv[0], nv[1], nv[2], rotateAngle);
	//printf("rotateAngle = %lf, rotateAngle = %lf\n", rotateAngle, rotateAngle1);
	if (rotateAngle > 90){
		// 鞋子(0, 1, 0)和y軸夾角 避免鞋底朝上
		//printf("==========================================================\n");
		initPosition(&theta[0], LR);
		CompRTMatf(m1, theta[0]);
	}

	if (rotateAngle1 > 120){
		// 鞋尖(0, 0, 1)和z軸夾角 避免鞋尖朝後
		//printf("==========================================================\n");
		initPosition(&theta[0], LR);
		CompRTMatf(m1, theta[0]);
	}

	if (rotateAngle2 < 70){
		// 鞋尖(0, 0, 1)和y軸夾角 避免鞋尖朝上
		//printf("==========================================================\n");
		initPosition(&theta[0], LR);
		CompRTMatf(m1, theta[0]);
	}

	/*if (LR == FT_LEFT){
		printf("Left Position = (%f, %f, %f) (%f, %f, %f)\n", theta[0].x, theta[0].y, theta[0].z, theta[0].a, theta[0].b, theta[0].c);
	}
	else{
		printf("Right Position = (%f, %f, %f) (%f, %f, %f)\n", theta[0].x, theta[0].y, theta[0].z, theta[0].a, theta[0].b, theta[0].c);
	}*/

	memcpy( m, m1[0], 12 * sizeof( float ) );
	// save transformation
	memcpy( _prev, m, 12 * sizeof( float ) );
	//printf("Total: %d\n", clock() - t_total);
	_e = FLT_MAX;
	return e1;
}
/*
void   FootTracker::Tracker::GetDepthArray( float m[16], Point3f_Array *DArray ){
	if( gdimg->Empty() ) return ;

	Point3f center;
	
	// initialize
	
	DArray->n = 0;
	center.x = m[3] + 50.0f * m[0];
	center.y = m[7] + 50.0f * m[4];
	center.z = m[11] + 50.0f * m[8];
	
	// depth image segmentation
	_DepthImgSegment( center.x, center.y, center.z );
	// count points in radius
	DArray->n = _CountRadiusPoints( center.x, center.y, center.z, 180.0f );
	if (DArray->n == 0)	return;
		
	// record points
	_CompRadiusPoints( center.x, center.y, center.z, 180.0f, DArray->pt, DArray->n );

	return ;
}

float   FootTracker::Tracker::TrackMergeICP( MeshBuffer mesh, float m[16], Point3f_Array *DArray1, Point3f_Array *DArray2 ){
	if( gdimg->Empty() ) return 0;

	Point3f_Array mesh_array, depth_array, buf_array;
	Point3f center;
	float *ptr, m1[3][4], e1, icpe[100];
	bool *v_list;
	RT_Paramf theta[1];
	ICP_Paramf param;
	int i;

	// initialize
	mesh_array.n = mesh.n;
	depth_array.n = 0;
	center.x = m[3] + 50.0f * m[0];
	center.y = m[7] + 50.0f * m[4];
	center.z = m[11] + 50.0f * m[8];
	memset( theta, 0x00, 1 * sizeof( RT_Paramf ) );


	depth_array.n = DArray1->n + DArray2->n;
	// memory allocate from memory pool
	mesh_array.pt = ( Point3f* )Allocate( 2 * ( mesh_array.n + depth_array.n ) * sizeof( Point3f ) );
	if( mesh_array.pt == NULL ) {									// mesh_array.n
		_e = FLT_MAX;
		return 0;
	}
	depth_array.pt = mesh_array.pt + mesh_array.n;					// depth_array.n
	buf_array.pt = depth_array.pt + depth_array.n;					// mesh_array.n
	// record points
	_CompRadiusPoints( center.x, center.y, center.z, 150.0f, depth_array.pt, depth_array.n );
	memcpy( depth_array.pt, DArray1->pt, DArray1->n * sizeof( Point3f ) );
	memcpy( depth_array.pt + DArray1->n, DArray2->pt, DArray2->n * sizeof( Point3f ) );
	// trim the vertex that can't see 
	v_list = ( bool* )buf_array.pt;
	memset( v_list, 0x00, mesh.n * sizeof( bool ) );
	// record the points can be see
	for( i = 0; i < mesh.nf; i++ ){
		ptr = mesh.vn + 3 * mesh.f[i].normal[0];	// normal of first point
		if( m[8] * ptr[0] + m[9] * ptr[1] + m[10] * ptr[2] > -0.f ) v_list[ mesh.f[i].vertex[0] ] = true;
		ptr = mesh.vn + 3 * mesh.f[i].normal[1];	// normal of second point
		if( m[8] * ptr[0] + m[9] * ptr[1] + m[10] * ptr[2] > -0.f ) v_list[ mesh.f[i].vertex[1] ] = true;
		ptr = mesh.vn + 3 * mesh.f[i].normal[2];	// normal of third point
		if( m[8] * ptr[0] + m[9] * ptr[1] + m[10] * ptr[2] > -0.f ) v_list[ mesh.f[i].vertex[2] ] = true;
	}
	// save the points
	ptr = mesh.v;
	for( i = mesh_array.n = 0; i < mesh.n; i++, ptr += 3 ){
		//if( v_list[i] == false ) continue;
		memcpy( mesh_array.pt + mesh_array.n, ptr, 3 * sizeof( float ) );
		mesh_array.n++;
	}
	
	//	memcpy( mesh_array.pt, mesh.v, mesh.n * 3 * sizeof( float ) );
	// compute transform matrix
	//(5,10)
	param.iter = 5;
	param.trans.iter = 10;
	param.trans.err = .5f;
	param.trans.lr1 = 0.05f;
	param.trans.lr2 = 0.4f;
	memcpy( m1[0], m, 12 * sizeof( float ) );
	CompRTFromMatf( theta, m1 );
	e1=ICPf( mesh_array, depth_array, theta, param, buf_array, icpe );
	CompRTMatf( m1, theta[0] );
	
	memcpy( m, m1[0], 12 * sizeof( float ) );
	// save transformation
	memcpy( _prev, m, 12 * sizeof( float ) );


	_e = FLT_MAX;
	return e1;
}*/

//int   FootTracker::Tracker::TrackICP( MeshBuffer mesh, float m[16], float s[3] ){
//	if( gdimg->Empty() ) return 0;
//
//	Point3f_Array mesh_array, depth_array, buf_array;
//	Point3f center;
//	float *ptr, m1[3][4], e1;
//	bool *v_list;
//	RT_Paramf theta[1];
//	Vector3f scale = { 1.0f, 1.0f, 1.0f };
//	ICP_Paramf param;
//	int i;
//
//	// initialize
//	mesh_array.n = mesh.n;
//	depth_array.n = 0;
//	center.x = m[3] + 50.0f * m[0];
//	center.y = m[7] + 50.0f * m[4];
//	center.z = m[11] + 50.0f * m[8];
//	memset( theta, 0x00, 1 * sizeof( RT_Paramf ) );
//
//	// count points in radius
//	if( ( depth_array.n = _CountRadiusPoints( center.x, center.y, center.z, 150.0f ) ) == 0 ) return 0;
//	// memory allocate from memory pool
//	mesh_array.pt = ( Point3f* )Allocate( 2 * ( mesh_array.n + depth_array.n ) * sizeof( Point3f ) );
//	if( mesh_array.pt == NULL ) return 0;							// mesh_array.n
//	depth_array.pt = mesh_array.pt + mesh_array.n;					// depth_array.n
//	buf_array.pt = depth_array.pt + depth_array.n;					// mesh_array.n
//	// record points
//	_CompRadiusPoints( center.x, center.y, center.z, 150.0f, depth_array.pt, depth_array.n );
//
//	// trim the vertex that can't see
//	v_list = ( bool* )buf_array.pt;
//	memset( v_list, 0x00, mesh.n * sizeof( bool ) );
//	// record the points can be see
//	for( i = 0; i < mesh.nf; i++ ){
//		ptr = mesh.vn + 3 * mesh.f[i].normal[0];	// normal of first point
//		if( m[8] * ptr[0] + m[9] * ptr[1] + m[10] * ptr[2] > 0.f ) v_list[ mesh.f[i].vertex[0] ] = true;
//		ptr = mesh.vn + 3 * mesh.f[i].normal[1];	// normal of second point
//		if( m[8] * ptr[0] + m[9] * ptr[1] + m[10] * ptr[2] > 0.f ) v_list[ mesh.f[i].vertex[1] ] = true;
//		ptr = mesh.vn + 3 * mesh.f[i].normal[2];	// normal of third point
//		if( m[8] * ptr[0] + m[9] * ptr[1] + m[10] * ptr[2] > 0.f ) v_list[ mesh.f[i].vertex[2] ] = true;
//	}
//	// save the points
//	ptr = mesh.v;
//	for( i = mesh_array.n = 0; i < mesh.n; i++, ptr += 3 ){
//		if( v_list[i] == false ) continue;
//		memcpy( mesh_array.pt + mesh_array.n, ptr, 3 * sizeof( float ) );
//		mesh_array.n++;
//	}
//	param.iter = 2;
//	param.trans.iter = 500;
//	param.trans.err = 1.0f;
//	param.trans.lr1 = 0.01f;
//	param.trans.lr2 = 0.6f;
//	memcpy( m1[0], m, 12 * sizeof( float ) );
//	CompRTFromMatf( theta, m1 );
////	e1 = CompRTSFromPointsf( mesh_array, depth_array, theta, &scale, param, buf_array );
//	CompRTMatf( m1, theta[0] );
//	memcpy( m, m1[0], 12 * sizeof( float ) );
//	memcpy( s, &scale, 3 * sizeof( float ) );
//	// save transformation
//	memcpy( _prev, m, 12 * sizeof( float ) );
//	return depth_array.n;
//}

int   FootTracker::Tracker::TrackColorICP( MeshBuffer mesh, int color, float m[16], float icpe[100] ){
	if( gdimg->Empty() ) return 0;

	Point3f_Array mesh_array, depth_array, buf_array;
	Point3f center;
	float *ptr, m0[16], m1[3][4], e1;
	bool *v_list;
	RT_Paramf theta[1];
	ICP_Paramf param;
	int i;

	// initialize
	memcpy( m1[0], m, 12 * sizeof( float ) );
	CompRTFromMatf( theta, m1 );
	memcpy( m, _prev, 12 * sizeof( float ) );						// prev PCA
	if( ThisPos( color, 50, m0[3], m0[7], m0[11] ) == 0 ) {			// this PCA
		_e = FLT_MAX;
		return 0;
	}
	// movement flow
	e1 = exp( -( SQR( m0[3] - m[3] ) + SQR( m0[7] - m[7] ) + SQR( m0[11] - m[11] ) ) / 1600.0f );
	e1 = 1.0f - e1;
	theta[0].x = m[3] + e1 * ( m0[3] - m[3] );	
	theta[0].y = m[7] + e1 * ( m0[7] - m[7] );	
	theta[0].z = m[11] + e1 * ( m0[11] - m[11] );
	mesh_array.n = mesh.n;
	depth_array.n = 0;
	center.x = m0[3] + 20.0f * m1[0][0];
	center.y = m0[7] + 20.0f * m1[1][0];
	center.z = m0[11] + 20.0f * m1[2][0];

	// depth image segmentation
//	if( m[11] < -50 ) 
		_DepthImgSegment( center.x, center.y, center.z );
	// count points in radius
	if( ( depth_array.n = _CountRadiusPoints( center.x, center.y, center.z, 150.0f ) ) == 0 ) {
		_e = FLT_MAX;
		return 0;
	}
	// memory allocate from memory pool
	mesh_array.pt = ( Point3f* )Allocate( 2 * ( mesh_array.n + depth_array.n ) * sizeof( Point3f ) );
	if( mesh_array.pt == NULL ) {									// mesh_array.n
		_e = FLT_MAX;
		return 0;
	}
	depth_array.pt = mesh_array.pt + mesh_array.n;					// depth_array.n
	buf_array.pt = depth_array.pt + depth_array.n;					// mesh_array.n
	// record points
	_CompRadiusPoints( center.x, center.y, center.z, 150.0f, depth_array.pt, depth_array.n );
/**/
	// trim the vertex that can't see
	v_list = ( bool* )buf_array.pt;
	memset( v_list, 0x00, mesh.n * sizeof( bool ) );
	// record the points can be see
	for( i = 0; i < mesh.nf; i++ ){
		ptr = mesh.vn + 3 * mesh.f[i].normal[0];	// normal of first point
		if( m1[2][0] * ptr[0] + m1[2][1] * ptr[1] + m1[2][2] * ptr[2] > -0.2f ) v_list[ mesh.f[i].vertex[0] ] = true;
		ptr = mesh.vn + 3 * mesh.f[i].normal[1];	// normal of second point
		if( m1[2][0] * ptr[0] + m1[2][1] * ptr[1] + m1[2][2] * ptr[2] > -0.2f ) v_list[ mesh.f[i].vertex[1] ] = true;
		ptr = mesh.vn + 3 * mesh.f[i].normal[2];	// normal of third point
		if( m1[2][0] * ptr[0] + m1[2][1] * ptr[1] + m1[2][2] * ptr[2] > -0.2f ) v_list[ mesh.f[i].vertex[2] ] = true;
	}
	// save the points
	ptr = mesh.v;
	for( i = mesh_array.n = 0; i < mesh.n; i++, ptr += 3 ){
		if( v_list[i] == false ) continue;
		memcpy( mesh_array.pt + mesh_array.n, ptr, 3 * sizeof( float ) );
		mesh_array.n++;
	};
/**/
//	memcpy( mesh_array.pt, mesh.v, mesh.n * 3 * sizeof( float ) );
	mesh_array.n = mesh.n;
	// compute transform matrix
	param.iter = 3;
	param.trans.iter = 10;
	param.trans.err = 1.0f;
	param.trans.lr1 = 0.05f;
	param.trans.lr2 = 0.4f;
	e1 = ICPf( mesh_array, depth_array, theta, param, buf_array, icpe );
	CompRTMatf( m1, theta[0] );
	memcpy( m, m1[0], 12 * sizeof( float ) );
	// save transformation
	memcpy( _prev, m, 12 * sizeof( float ) );

	_e = FLT_MAX;
	return depth_array.n;
}

#define WIN 45

float   FootTracker::Tracker::Residual( void ){
	return _e;
}



void FootTracker::Tracker::InitialCam( float mcam[16] )
{
	mcam[ 0]=1; mcam[ 1]=0; mcam[ 2]=0; mcam[ 3]=0;
	mcam[ 4]=0; mcam[ 5]=1; mcam[ 6]=0; mcam[ 7]=0;
	mcam[ 8]=0; mcam[ 9]=0; mcam[10]=1; mcam[11]=0;
	mcam[12]=0; mcam[13]=0; mcam[14]=0; mcam[15]=1;
}




//
// OpenGL Setting
//

void   FootTracker::Draw( void ){
	if( gdimg->Empty() ) return;
	printf("///////////////");
	const float m[16] = {
		1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		0.0f, 0.0f, -1.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	};
	const float diffuse[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
	const float position[4] = { 0.0f, 100.0f, 1.0f, 1.0f };
	int *pdepth[4], *dbits, *ptex[4], *tbits;
	int i, j;
	Point3f buf;

	dbits = gdimg->GetImage();
	tbits = gtexture->GetImage();
	
	glPushMatrix();
	glMultMatrixf( m );
	glEnable( GL_LIGHTING );
	glLightfv( GL_LIGHT0, GL_DIFFUSE, diffuse ); 
	glLightfv( GL_LIGHT0, GL_AMBIENT, diffuse ); 
	glLightfv( GL_LIGHT0, GL_POSITION, position ); 
	glEnable( GL_LIGHT0 );
	glDisable( GL_LIGHTING );

	glBegin( GL_TRIANGLES );
	for(i = 1; i < 1080; i++) {
		pdepth[0] = dbits + ( i - 1 ) * 1920;
		pdepth[1] = dbits + (i - 1) * 1920 + 1;
		pdepth[2] = dbits + i * 1920;
		pdepth[3] = dbits + i * 1920 + 1;
		ptex[0] = tbits + ( i - 1 ) * 1920;
		ptex[1] = tbits + ( i - 1 ) * 1920 + 3;
		ptex[2] = tbits + i * 1920;
		ptex[3] = tbits + i * 1920 + 3;
		for (j = 1; j < 1920; j++, pdepth[0]++, pdepth[1]++, pdepth[2]++, pdepth[3]++, ptex[0] += 3, ptex[1] += 3, ptex[2] += 3, ptex[3] += 3){
			if( *pdepth[1] < 50 || *pdepth[2] < 50 ) continue;
			if( *pdepth[0] > 50 ){
				glColor3f( ( float )ptex[0][0] / 255.0f, ( float )ptex[0][1] / 255.0f, ( float )ptex[0][2] / 255.0f );
				buf = gdimg->Comp3DPosf( j - 1, i - 1 );
				glVertex3f( buf.x, buf.y, buf.z );
				glColor3f( ( float )ptex[1][0] / 255.0f, ( float )ptex[1][1] / 255.0f, ( float )ptex[1][2] / 255.0f );
				buf = gdimg->Comp3DPosf( j, i - 1 );
				glVertex3f( buf.x, buf.y, buf.z );
				glColor3f( ( float )ptex[2][0] / 255.0f, ( float )ptex[2][1] / 255.0f, ( float )ptex[2][2] / 255.0f );
				buf = gdimg->Comp3DPosf( j - 1, i );
				glVertex3f( buf.x, buf.y, buf.z );
			}
			if( *pdepth[3] > 50 ){
				glColor3f( ( float )ptex[1][0] / 255.0f, ( float )ptex[1][1] / 255.0f, ( float )ptex[1][2] / 255.0f );
				buf = gdimg->Comp3DPosf( j, i - 1 );
				glVertex3f( buf.x, buf.y, buf.z );
				glColor3f( ( float )ptex[3][0] / 255.0f, ( float )ptex[3][1] / 255.0f, ( float )ptex[3][2] / 255.0f );
				buf = gdimg->Comp3DPosf( j, i );
				glVertex3f( buf.x, buf.y, buf.z );
				glColor3f( ( float )ptex[2][0] / 255.0f, ( float )ptex[2][1] / 255.0f, ( float )ptex[2][2] / 255.0f );
				buf = gdimg->Comp3DPosf( j - 1, i );
				glVertex3f( buf.x, buf.y, buf.z );
			}
		}
	}
	glEnd();
	glDisable( GL_LIGHTING );
	glPopMatrix();
}

void   FootTracker::CompGLMatf( float mf[16], float m[16] ){
	mf[0] = m[0];	mf[1] = m[4];	mf[2] = m[8];	mf[3] = 0.0f;
	mf[4] = m[1];	mf[5] = m[5];	mf[6] = m[9];	mf[7] = 0.0f;
	mf[8] = m[2];	mf[9] = m[6];	mf[10] = m[10];	mf[11] = 0.0f;
	mf[12] = m[3];	mf[13] = m[7];	mf[14] = m[11];	mf[15] = 1.0f;
}

void   FootTracker::CompGLMatd( double mf[16], double m[16] ){
	mf[0] = m[0];	mf[1] = m[4];	mf[2] = m[8];	mf[3] = 0.0;
	mf[4] = m[1];	mf[5] = m[5];	mf[6] = m[9];	mf[7] = 0.0;
	mf[8] = m[2];	mf[9] = m[6];	mf[10] = m[10];	mf[11] = 0.0;
	mf[12] = m[3];	mf[13] = m[7];	mf[14] = m[11];	mf[15] = 1.0;
}

void   FootTracker::DrawPit( int r ){
	double k1, k2, gauss, mean;
	const int h = gdimg->Height() - r, w = gdimg->Width() - r;
	Point3f buf;
	int i, j;

	glBegin( GL_POINTS );
	for( i = r; i < h; i++ ){
		for( j = r; j < h; j++ ){
			gdimg->Curvature( &k1, &k2, j, i, r );
			mean = 0.5 * ( k1 + k2 );
			gauss = k1 * k2;
			if( mean < 0.0 && gauss < 0.0 ) {
				buf = gdimg->Comp3DPosf( j, i );
				glVertex3f( buf.x, buf.y, buf.z );
			}
		}
	}
	glEnd();
}

int FootTracker::Tracker::test(float m[16]){

	Point3f_Array mesh_array, depth_array, buf_array;
	Point3f center;
	float *ptr, m1[3][4], e1;
	bool *v_list;
	RT_Paramf theta[1];
	ICP_Paramf param;
	int i;

	// initialize
	depth_array.n = 0;
	center.x = m[3] + 50.0f * m[0];
	center.y = m[7] + 50.0f * m[4];
	center.z = m[11] + 50.0f * m[8];
	memset(theta, 0x00, 1 * sizeof(RT_Paramf));

	// depth image segmentation
	//	if( m[11] < -50 ) 
	_DepthImgSegment(center.x, center.y, center.z); // 深度切割 傳入中心點	

	return 0;
}

Point3f FootTracker::DepthOfAnkle(float x, float y, float z, float r, float tx, float ty, float tz, float color, float calibx, float caliby, float calibz){

	const int
		w = gdimg->Width(),
		h = gdimg->Height();

	const float r2 = SQR(r);

	Point3f temp, pt;
	int i, j;
	Point3f tempsum = { 0 };

	for (i = 0; i < h; i += 3){
		for (j = 0; j < w; j += 3){

			// compute distance to center
			pt = gdimg->Comp3DPosf(j, i);
			temp.x = pt.x - x;
			temp.y = pt.y - y;
			temp.z = pt.z - z;

			if (SQR(temp.x) + SQR(temp.y) + SQR(temp.z) > r2) continue;	// counting

			pt.x += (calibx + 10);
			pt.y += (caliby);
			pt.z += (calibz);
			
			// draw foot angle point
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glPointSize(2.75);
			glBegin(GL_POINTS);
			glColor4f(0.0f, 0.0f, 0.0f, 0.5f);
			glVertex3f(pt.x, pt.y, pt.z);
			glEnd();
			glPointSize(1.0);

		}
	}

	Point3f line_vector = { tempsum.x, tempsum.y, tempsum.z};
	return line_vector;

}