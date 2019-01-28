#ifndef _IMGPROC_
#define _IMGPROC_

#include <windows.h>

// Warping parameters //

#define IMG_NEAREST		0x01
#define IMG_BILINEAR	0x02
#define IMG_HERMITE		0x04

// Image Space //

#define IMG_GREY		0x00000000
#define IMG_RED			0x00000001
#define IMG_BLUE		0x00000010
#define IMG_GREEN		0x00000100
#define IMG_RGB			0x00000111
#define IMG_X			0x00001000
#define IMG_Y			0x00010000

// Image Selection //

#define IMG_SELECTED	1
#define IMG_UNSELECTED	0
#define IMG_BOUNDARY	-1

//#ifdef __cplusplus

//               //
// Image process //
//               //

void AveImage( float *img, float *ave, int w, int h );
void Ave2Image( float *img, float *ave2, int w, int h );
void VarImage( float *ave, float *ave2, float *var, int w, int h );

void DrawLine( BYTE *img, int w, int h, int stride, int x1, int y1, int x2, int y2, int r, int g, int b );
void DrawLine( BITMAP *img, int x1, int y1, int x2, int y2, int r, int g, int b );

void Dither( BYTE *img, int w, int h, int stride, double* threshold, int m_w, int m_h );	// need to free it's output
void Dither( BITMAP *img, double* threshold, int m_w, int m_h );	// need to free it's output

void Contrast( BYTE *img, int w, int h, int stride, BYTE threshold, double ratio );	// need to free it's output
void Contrast( BITMAP *img, BYTE threshold, double ratio );	// need to free it's output

void Convolution( BYTE *img, int w, int h, int stride, int size, double *kernel );
void Convolution( BYTE *img, int w, int h, int stride, int size, float *kernel );
void Convolution( BITMAP *img, int size, double *kernel );
void Convolution( BITMAP *img, int size, float *kernel );

void GaussFilter( BYTE *img, int w, int h, int stride, int win );
void GaussFilter( BITMAP *img, int win );

void GradientX( BYTE *img, int w, int h, int stride, int *g );
void GradientY( BYTE *img, int w, int h, int stride, int *g );
void HoG( BYTE *img, int w, int h, int stride, int win, int *hog );

void BilateralFilterGrey( BYTE *, int w, int h, int win, int segment );			// only for 8-bit grey image
void BilateralFilterGrey( BYTE *, BYTE *, int w, int h, int win, int segment );	// only for 8-bit grey image
void BilateralFilterGrey( float *img, float *src, float *ref, int w, int h, int win, int segment );
void BilateralFilterGrey( double *img, double *src, double *ref, int w, int h, int win, int segment );
void BilateralFilter( BYTE *img, int w, int h, int stride, int win, int segment );
void BilateralFilter( BYTE *src, BYTE *tar, int w, int h, int stride, int win, int segment );
void BilateralFilter( BITMAP *img, int win, int segment );
int BilateralFilter( BITMAP *src, const BITMAP &tar, int win, int segment );

void Histogram( BITMAP *img );
void Histogram( BYTE *img, int w, int h, int stride );
void Histogram( BITMAP *source, const BITMAP &target, int space );
void Histogram( BYTE *source, int w1, int h1, int stride1, BYTE *target, int w2, int h2, int stride2, int space );

void Warping( BYTE* img, int w, int h, int stride, double H[], char param = IMG_NEAREST );
void Warping( BITMAP *img, double H[9], char param = IMG_NEAREST );

void Select( BYTE *img, int w, int h, int stride, int x, int y );
void Select( BITMAP *img, int x, int y );
void Select( BYTE *img, char *region, int w, int h, int stride, int p[], int n );
void Select( BITMAP *img, char *region, int p[], int n );
void Select( BYTE *img, int w, int h, int stride, int p[], int n );
void Select( BITMAP *img, int p[], int n );

void ColorTransfer( BYTE *img1, int w1, int h1, int stride1, BYTE *img2, int w2, int h2, int stride2 );
void ColorTransfer( BITMAP *source, const BITMAP &target );

void Scale( BYTE *img1, int w1, int h1, int stride1, BYTE *img2, int w2, int h2, int stride2 );
void Scale8( float *img1, int w1, int h1, float *img2, int w2, int h2 );
void Scale8( double *img1, int w1, int h1, double *img2, int w2, int h2 );

void Colorization( BYTE *img, int w, int h, int stride, int scribble[], int color[], int n );
void Colorization( BYTE *img1, int w1, int h1, int stride1, BYTE *img2, int w2, int h2, int stride2, 
				 int scrib[], int col[], int n );
void Colorization( BITMAP *img, int scribble[], int color[], int n );
void Colorization( BITMAP *img, const BITMAP &scrib );
void Colorization( BITMAP *img, const BITMAP &sml, int scribble[], int color[], int n );
void Colorization( BITMAP *img, const BITMAP &sml, const BITMAP &scrible ); 

// Seam Carving
void SeamCarving( BYTE *img, int w1, int h1, int stride, int w2, int h2 );	// resizing
void SeamEnlarge( BYTE *img1, int w1, int h1, int stride1, 
				 BYTE *img2, int w2, int h2, int stride2 );	// enlarging
void SeamResize( BYTE *img1, int w1, int h1, int stride, BYTE *img2, int w2, int h2, int stride2 );
void SeamObjRemove( BYTE *img1, BYTE *img2, int w, int h, int stride );

// Poisson Image Editing
int SeamlessClone( BYTE *img1, int w1, int h1, int stride1, BYTE *img2, int w2, int h2, int stride2,
				   int p[], int n, int offsetx, int offsety );
int MixSeamlessClone( BYTE *img1, int w1, int h1, int stride1, BYTE *img2, int w2, int h2, int stride2,
				   int p[], int n, int offsetx, int offsety );

// Matte
int Matte( BYTE *img, BYTE *alpha, int w, int h, int stride, int f[], int nf, int b[], int nb );
int Matte( BYTE *img, BYTE *scrib, BYTE *alpha, int w, int h, int stride );

//
// imgproc2.cpp
//

#include "imgconf.h"

//
// Image Painting
//

int DrawLine1b( ImageBuffer buf, ImagePoint, ImagePoint, Color1b );
int DrawLine3b( ImageBuffer buf, ImagePoint, ImagePoint, Color3b );
int DrawLine4b( ImageBuffer buf, ImagePoint, ImagePoint, Color4b );
int DrawLine1i( ImageBuffer buf, ImagePoint, ImagePoint, Color1i );
int DrawLine3i( ImageBuffer buf, ImagePoint, ImagePoint, Color3i );
int DrawLine4i( ImageBuffer buf, ImagePoint, ImagePoint, Color4i );
int DrawLine1f( ImageBuffer buf, ImagePoint, ImagePoint, Color1f );
int DrawLine3f( ImageBuffer buf, ImagePoint, ImagePoint, Color3f );
int DrawLine4f( ImageBuffer buf, ImagePoint, ImagePoint, Color4f );
int DrawLine1d( ImageBuffer buf, ImagePoint, ImagePoint, Color1d );
int DrawLine3d( ImageBuffer buf, ImagePoint, ImagePoint, Color3d );
int DrawLine4d( ImageBuffer buf, ImagePoint, ImagePoint, Color4d );

int DrawTriangle1b( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color1b, Color1b, Color1b );
int DrawTriangle3b( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color3b, Color3b, Color3b );
int DrawTriangle4b( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color4b, Color4b, Color4b );
int DrawTriangle1i( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color1i, Color1i, Color1i );
int DrawTriangle3i( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color3i, Color3i, Color3i );
int DrawTriangle4i( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color4i, Color4i, Color4i );
int DrawTriangle1f( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color1f, Color1f, Color1f );
int DrawTriangle3f( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color3f, Color3f, Color3f );
int DrawTriangle4f( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color4f, Color4f, Color4f );
int DrawTriangle1d( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color1d, Color1d, Color1d );
int DrawTriangle3d( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color3d, Color3d, Color3d );
int DrawTriangle4d( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color4d, Color4d, Color4d );

int Filling1b( ImageBuffer, Region, Color1b );
int Filling3b( ImageBuffer, Region, Color3b );
int Filling4b( ImageBuffer, Region, Color4b );
int Filling1i( ImageBuffer, Region, Color1i );
int Filling3i( ImageBuffer, Region, Color3i );
int Filling4i( ImageBuffer, Region, Color4i );
int Filling1f( ImageBuffer, Region, Color1f );
int Filling3f( ImageBuffer, Region, Color3f );
int Filling4f( ImageBuffer, Region, Color4f );
int Filling1d( ImageBuffer, Region, Color1d );
int Filling3d( ImageBuffer, Region, Color3d );
int Filling4d( ImageBuffer, Region, Color4d );

//
// Editing
//

int Selectb( ImageBuffer, Selected, ImagePoint, float threshold );
int Selecti( ImageBuffer, Selected, ImagePoint, float threshold );
int Selectf( ImageBuffer, Selected, ImagePoint, float threshold );
int Selectd( ImageBuffer, Selected, ImagePoint, float threshold );

//
// Image Processing
//

int Convolutionb( ImageBuffer buf, int win, float *kernel );
int Convolutioni( ImageBuffer buf, int win, float *kernel );
int Convolutionf( ImageBuffer buf, int win, float *kernel );
int Convolutiond( ImageBuffer buf, int win, float *kernel );

int GaussFilterb( ImageBuffer buf, int win );
int GaussFilteri( ImageBuffer buf, int win );
int GaussFilterf( ImageBuffer buf, int win );
int GaussFilterd( ImageBuffer buf, int win );

int Scaleb( ImageBuffer src, ImageBuffer tar );
int Scalei( ImageBuffer src, ImageBuffer tar );
int Scalef( ImageBuffer src, ImageBuffer tar );
int Scaled( ImageBuffer src, ImageBuffer tar );

void Ditherb( ImageBuffer, unsigned char * );
void Ditheri( ImageBuffer, int *, int max, int min );
void Ditherf( ImageBuffer, float *, float max, float min );
void Ditherd( ImageBuffer, double *, double max, double min );

//
// Color Converting
//

int RGBtoNTSCi( ImageBuffer );
int RGBtoNTSCf( ImageBuffer );
int RGBtoNTSCd( ImageBuffer );

int NTSCtoRGBi( ImageBuffer );
int NTSCtoRGBf( ImageBuffer );
int NTSCtoRGBd( ImageBuffer );

#endif