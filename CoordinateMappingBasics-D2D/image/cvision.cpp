
#include "bitmap.h"
#include "cvision.h"
#include <float.h>
#include <math.h>
#include <list>
#include "algo.h"
using namespace std;

#define SQUARE(x)	(x)*(x)
#define SQR(x)		SQUARE(x)
#define ERR_BOUND	0.0000000001
#define M_PI		3.14159265358979323846

//
// inline functions
//

inline void SetGrayPixel( BITMAP &bmp, int x, int y, BYTE value ){
	BYTE *bits = (BYTE*) bmp.bmBits;
	bits[ y * bmp.bmWidthBytes + x * 3 + 0 ] = 
	bits[ y * bmp.bmWidthBytes + x * 3 + 1 ] = 
	bits[ y * bmp.bmWidthBytes + x * 3 + 2 ] = value; 
}

inline int GrayPixel( const BITMAP &bmp, int x, int y ){
	BYTE *bits = (BYTE*) bmp.bmBits;
	return (
		(int)bits[ y * bmp.bmWidthBytes + x * 3 + 0 ] + 
		(int)bits[ y * bmp.bmWidthBytes + x * 3 + 1 ] + 
		(int)bits[ y * bmp.bmWidthBytes + x * 3 + 2 ]
		) / 3;
}

inline void SetPixel( BITMAP &bmp, int x, int y, BGR c, BYTE value ){
	BYTE *bits = (BYTE*) bmp.bmBits;
	bits[ y * bmp.bmWidthBytes + x * 3 + (int)c ] = value;
}

inline int Pixel( const BITMAP &bmp, int x, int y, BGR c ){
	BYTE *bits = (BYTE*) bmp.bmBits;
	return (int)bits[ y * bmp.bmWidthBytes + x * 3 + (int)c ];

}

inline double Pixel( const BITMAP &bmp, int x, int y, XYZ c ){
	double b, g, r;
	b = (double)Pixel( bmp, x, y, B ) / 255.0;
	g = (double)Pixel( bmp, x, y, G ) / 255.0;
	r = (double)Pixel( bmp, x, y, R ) / 255.0;
	switch( c ){
	case X: return ( 0.49 * r + 0.31 * g + 0.2 * b ) / 0.17697;
	case Y: return ( 0.17697 * r + 0.81240 * g + 0.01063 * b ) / 0.17697;
	case Z: return ( 0.00 * r + 0.01 * g + 0.88 * b ) / 0.17697;
	}
	return 0.0;
}

inline double Pixel( const BITMAP &bmp, int x, int y, LUV c ){
#define u (4.0*Xc/(Xc+15.0*Yc+3.0*Zc))
#define v (9.0*Yc/(Xc+15.0*Yc+3.0*Zc))
#define un 0.2009
#define vn 0.4610
	double Xc, Yc, Zc, l;
	Xc = Pixel( bmp, x, y, X );
	Yc = Pixel( bmp, x, y, Y );
	Zc = Pixel( bmp, x, y, Z );
	l = ( Yc > 0.00885645 ) ? ( 116.0 * pow( Yc, 0.333 ) - 16.0 ) : 903.2962963 * Yc;
	switch( c ){
	case L: return l;
	case U: return 13.0 * l * ( u - un );
	case V: return 13.0 * l * ( v - vn );
	}
	return 0.0;
#undef u
#undef v
#undef vn
#undef un
}

inline double SobelX( BITMAP& img, int x, int y ){
	return (
		1.0 * GrayPixel( img, x+1, y+1 ) + 2.0 * GrayPixel( img, x+1, y ) + 1.0 * GrayPixel( img, x+1, y-1 ) -
		1.0 * GrayPixel( img, x-1, y+1 ) - 2.0 * GrayPixel( img, x-1, y ) - 1.0 * GrayPixel( img, x-1, y-1 ) 
		) / 8.0;
}
inline double SobelY( BITMAP& img, int x, int y ){
	return (
		1.0 * GrayPixel( img, x+1, y+1 ) + 2.0 * GrayPixel( img, x, y+1 ) + 1.0 * GrayPixel( img, x-1, y+1 ) -
		1.0 * GrayPixel( img, x+1, y-1 ) - 2.0 * GrayPixel( img, x, y-1 ) - 1.0 * GrayPixel( img, x-1, y-1 ) 
		) / 8.0;
}
inline double funcR( double A[2][2], double k){
	return ( A[0][0] * A[1][1] - A[1][0] * A[0][1] ) - k * SQUARE( ( A[0][0] + A[1][1] ) );
}

inline bool NonMax( double **grad, double **dx, double **dy, int i, int j ){
	double x = dx[i][j];
	double y = dy[i][j];			// vector (x, y) on the gradient
	double temp = max( x, y );
	x /= temp; y /= temp;
	if( x == 1.0 ){
		if( y > 0.0 ){
			if( ( 1.0 - y ) * grad[i][j + 1] + y * grad[i + 1][j + 1] < grad[i][j] &&		// local maximum?
				( 1.0 - y ) * grad[i][j - 1] + y * grad[i - 1][j - 1] < grad[i][j] ) 
				return true;
		}
		else{
			if( ( 1.0 + y ) * grad[i][j + 1] + -(y) * grad[i - 1][j + 1] < grad[i][j] &&
				( 1.0 + y ) * grad[i][j - 1] + -(y) * grad[i + 1][j - 1] < grad[i][j] ) 
				return true;
		}
	}else{
		if( x > 0.0 ){
			if( ( 1.0 - x ) * grad[i + 1][j] + x * grad[i + 1][j + 1] < grad[i][j] &&
				( 1.0 - x ) * grad[i - 1][j] + x * grad[i - 1][j - 1] < grad[i][j] ) 
				return true;
		}
		else{
			if( ( 1.0 + x ) * grad[i + 1][j] + -(x) * grad[i + 1][j - 1] < grad[i][j] &&
				( 1.0 + x ) * grad[i - 1][j] + -(x) * grad[i - 1][j + 1] < grad[i][j] ) 
				return true;
		}
	}
	return false;
}

void EdgeDetect( BITMAP &img, double threshold, bool onlyedge ){	
	int i, j;
	int h = img.bmHeight;
	int w = img.bmWidth;
	double **dx, **dy, **grad;

	// initialize the dx, dy, and gradient magnitude matrix
	dx = new double* [h];
	dy = new double* [h];
	grad = new double* [h];
	for( i = 0; i < h; i++ ){
		dx[i] = new double [w];
		dy[i] = new double [w];
		grad[i] = new double [w];
		for( j = 0; j < w; j++ ) dx[i][j] = dy[i][j] = grad[i][j] = 0.0;
	}

	for( i = 1; i < h - 1; i++ ) {
		for( j = 1; j < w - 1; j++ ) {
			// 1B: Apply Sobel Operator		
			dx[i][j] = SobelX( img, j, i );
			dy[i][j] = SobelY( img, j, i );

			// 1C: Compute the gradient magnitude
			grad[i][j] = sqrt( SQUARE( dx[i][j] ) + SQUARE( dy[i][j] ) );
		}
	}

	// 1D: Non-maximum suppression
	for( i = 0; i < h; i++ ){
		for( j = 0; j < w; j++ ){
			if( i == 0 || j == 0 || i == h - 1 || j == w - 1 ){
				SetGrayPixel( img, j, i, 0xff );
				continue;
			}
			SetGrayPixel( img, j, i, 0xFF );
			if( grad[i][j] > threshold ){
				if( !onlyedge || NonMax( grad, dx, dy, i, j ) ) {
					// 1E: Extract straight line segments
					SetGrayPixel( img, j, i, 0x00 );
				}
			}
		}
	}

	// delete matrices dx and dy
	for( i = 0; i < h; i++ ){
		delete dx[i];
		delete dy[i];
		delete grad[i];
	}
	delete dx;
	delete dy;
	delete grad;
}

void CornerDetect( BITMAP &bmp, POINT **c, int *n_corner, int nw, double threadshold ){
	CornerDetect( bmp, c, n_corner, nw, nw, nw, bmp.bmWidth - nw, bmp.bmHeight - nw,threadshold );
}

void CornerDetect( BITMAP &bmp, POINT **c, int *n_corner, int nw, int left, int bottom, int right, int up, double threadshold ){
	// Ix[i][j]: the deviation by dx at I[i][j]
	// Iy[i][j]: the deviation by dy at I[i][j]
	double **Ix, **Iy, **R;
	double H[2][2], temp;
	int i, j, ii, jj, h, w, n;
	int height, width;
	list<POINT> cor;								// Linked list of corner coordinate 
	list< POINT >::iterator	it = cor.begin();		// Iterator of list
	BYTE* pixel;									// array of pixel data

	h = bmp.bmHeight;							// Height of image
	w = bmp.bmWidth;							// Width of image
	pixel = (BYTE*) bmp.bmBits;
	cor.push_back( POINT() );
	*n_corner = 0;
	
	// Check the input rectangle and threashold
	threadshold = abs( threadshold );
	if( right < left ) swap( right, left );
	if( bottom > up ) swap( right, left );
	if( left < nw || left >= w - nw ) left = nw;
	if( right < nw || right >= w - nw ) right = w - 1 - nw;
	if( bottom < nw || bottom >= h - nw ) bottom = nw;
	if( up < nw || up >= h - nw ) up = h - 1 - nw;

	height = up - bottom + 2 * nw;
	width = right - left + 2 * nw;

	// initialize matrices
	Ix = new double* [height];
	Iy = new double* [height];
	R = new double* [height];
	for( i = 0; i < height; i++ ) {
		Ix[i] = new double [width]; 
		Iy[i] = new double [width];
		R[i] = new double [width];
		for( j = 0; j < width; j++ ) R[i][j] = 0.0;
	}
	
	// Get all of the deviation values 
	for( i = 0; i < height; i++ ) {
		for( j = 0; j < width; j++ ) {
			if( i + bottom - nw > 0 && i + bottom - nw < h-1 && j + left - nw > 0 && j + left < w-1 ) {
				if( GrayPixel( bmp, j + left - nw, i + bottom - nw ) > 0x40 ) Ix[i][j] = Iy[i][j] = 0.0;
				else{
					Ix[i][j] = SobelX( bmp, j + left - nw, i + bottom - nw );
					Iy[i][j] = SobelY( bmp, j + left - nw, i + bottom - nw );
				}
			}
			else{
				if( j + left - nw == 0 ) Ix[i][j] = 
					0.5 * GrayPixel( bmp, j + left - nw+1, i + bottom - nw ) - 
					0.5 * GrayPixel( bmp, j + left - nw, i + bottom - nw );
				else if( j + left - nw == w-1 ) Ix[i][j] = 
					0.5 * GrayPixel( bmp, j + left - nw, i + bottom - nw ) - 
					0.5 * GrayPixel( bmp, j + left - nw-1, i + bottom - nw );
				else Ix[i][j] = 
					0.5 * GrayPixel( bmp, j + left - nw+1, i + bottom - nw ) - 
					0.5 * GrayPixel( bmp, j + left - nw-1, i + bottom - nw );
				if( i + bottom - nw == 0 ) Iy[i][j] = 
					0.5 * GrayPixel( bmp, j + left - nw, i + bottom - nw+1 ) -
					0.5 *  GrayPixel( bmp, j + left - nw, i + bottom - nw );
				else if( i + bottom - nw == h-1 ) Iy[i][j] = 
					0.5 * GrayPixel( bmp, j + left - nw, i + bottom - nw ) - 
					0.5 * GrayPixel( bmp, j + left - nw, i + bottom - nw-1 );
				else Iy[i][j] = 
					0.5 * GrayPixel( bmp, j + left - nw, i + bottom - nw+1 ) - 
					0.5 * GrayPixel( bmp, j + left - nw, i + bottom - nw-1 );
			}
		}
	}

#define A 1.0
#define B ( -H[1][1] - H[0][0] )
#define C ( H[0][0] * H[1][1] - H[1][0] * H[0][1] )
	// Corner Detector
	for( i = 0 ; i < height ; i++ ) {
		for( j = 0 ; j < width; j++ ) {
			// Calculate H
			H[0][0] = H[0][1] = H[1][0] = H[1][1] = 0.0;		// Initialize H = 0
			n = 0;												// count the number of iteration
			w = ( int )( nw * 0.5f );
			for( ii = -w; ii < w; ii++ ) {
				for( jj = -w ; jj < w; jj++ ) {
					if( i + ii >= 0 && i + ii < height && j + jj >= 0 && j + jj < width ){
						H[0][0] += Ix[i + ii][j + jj] * Ix[i + ii][j + jj];
						H[1][0] += Ix[i + ii][j + jj] * Iy[i + ii][j + jj];
						H[1][1] += Iy[i + ii][j + jj] * Iy[i + ii][j + jj];
						n++;
					}
				}
			}
			// A = H/n
			H[0][0] /= (double)n;
			H[1][0] /= (double)n;
			H[1][1] /= (double)n;
			H[0][1] = H[1][0];			// for the off-diagonal elements of H are equal
			temp = fabs( -B - sqrt( B * B - 4.0 * A * C ) ) / (2.0 * A );	
			if( temp > threadshold && temp > R[i][j] ) 
				R[i][j] = temp;		// test by threshold and choose the max one along RGB colors
		}
	}
#undef A
#undef B
#undef C

	// delete matrices Ix and Iy
	for( i = 0; i < height; i++ ){
		delete Ix[i];
		delete Iy[i];
	}
	delete Ix;
	delete Iy;

	// Find local maximum
	for( i = nw; i < height - nw; i++ ) {
		for( j = nw; j < width - nw; j++ ) {
			for( ii = -1; ii <= 1; ii++ ){					// compare the I(i, j) with 8 of its neighbors
				for( jj = -1; jj <= 1; jj++ ){
					if( ii == 0 && jj == 0 ) continue;
					else if( R[i][j] > R[i + ii][j + jj] ) continue;
					else break;
				}
				if( jj < 2 ) break;							// it isn't a local maximum
			}
			if( ii == 2 ){									// it is a local maximum
				it = cor.begin();
				POINT buf={ j + left - nw, i + bottom - nw };
				if( *n_corner > 0 && R[i][j] > R[it->y + nw - bottom][it->x + nw - left] ) cor.push_front( buf );
				else{
					cor.pop_back();
					cor.push_back( buf );
					cor.push_back( POINT() );
				}
				(*n_corner)++;
			}
		}
	}

	// copy the linked list to the input array referance
	it = cor.begin();
	*c = new POINT [ *n_corner ];
	for( i = 0; i < *n_corner; i++, it++ ) (*c)[i] = (*it);

	// Sorting
	for( i = 0; i < *n_corner; i++ ){
		for( j = i + 1; j < *n_corner; j++ ){
			if( R[(*c)[i].y + nw - bottom][(*c)[i].x + nw - left] < R[(*c)[j].y + nw - bottom][(*c)[j].x + nw - left] )
				swap( (*c)[i], (*c)[j] );
		}
	}

	// delete matrix R
	for( i = 0; i < height; i++ ) delete R[i];
	delete R;
	cor.clear();
}

inline float GetMinor( POINT &p1, POINT &p2, POINT &p, float x, float y, float a, float min ){
	float b, d2, f2, cos2, sin2, temp1, temp2, e, e1, e2;
	// checking
	temp1 = ( float )( p2.x - p1.x );	temp2 = ( float )( p2.y - p1.y );	e = SQR( temp1 ) + SQR( temp2 );
	temp1 = ( float )( p.x - p1.x );	temp2 = ( float )( p.y - p1.y );	e1 = SQR( temp1 ) + SQR( temp2 );
	temp1 = ( float )( p2.x - p.x );	temp2 = ( float )( p2.y - p.y );	e2 = SQR( temp1 ) + SQR( temp2 );
	if( e + e1 <= e2 || e + e2 <= e1 ) return 0.0f;
	// d2
	temp1 = ( float )p.x - x; temp2 = ( float )p.y - y;
	d2 = SQR( temp1 ) + SQR( temp2 );
	if( d2 < SQR( min ) ) return 0.0f;			// distance is too small
	// f2
	f2 = e2;
	// cos and sin
	temp1 = a * a;
	cos2 = ( temp1 + d2 - f2 ) / ( 2.0f * a * sqrt( d2 ) );
	cos2 = SQR( cos2 );
	sin2 = 1.0f - cos2;
	// b
	temp2 = temp1 - d2 * cos2;
	if( temp2 < 0.0001f ) return 0.0f;		// divided by 0
	b = ( temp1 * d2 * sin2 ) / temp2;
	return sqrt(b);
}

void EllipseDetect( BITMAP &img, int left, int bottom, int right, int top, 
				   float *x0, float *y0, float *a0, float *b0, float *angle, int threshold, float min ){
	POINT *edge;
	int i, j, k, l, n, na, *accumulator, max = threshold;
	float x, y, o, a, b, temp; 					// center, orientation, major and minor axices
	BYTE *bits = ( BYTE* )img.bmBits;
	int it = 0;

	// count edge pixel number
	n = 0;
	for( i = left; i <= right; i++ ){
		for( j = bottom; j <= top; j++ ){
			if( bits[j * img.bmWidthBytes + i * 3] == 0x00 ) n++;
		}
	}

	// memory allocate
	edge = ( POINT* )malloc( n * sizeof( POINT ) );
	na = max( right - left, top - bottom ) / 2;
	accumulator = ( int* )malloc( na * sizeof( int ) );

	// set edge pixel array
	k = 0;
	for( i = 0; i < img.bmWidth; i++ ){
		for( j = 0; j < img.bmHeight; j++ ){
			if( bits[j * img.bmWidthBytes + i * 3] == 0x00 ) {
				edge[ k ].x = i;
				edge[ k++ ].y = j;
			}
		}
	}

	// hough transform
	while( it++ < 1000 ){
//	for( i = 0; i < n; i++ ){				// pick a pair of edge pixel
		i = rand() % n;						// RANSEC
		j = rand() % n;
		while( i == j ) { j = rand() % n; }
//		for( j = i + 1; j < n; j++ ){
			// initialize accumulator array
			for( k = 0; k < na; k++ ) accumulator[k] = 0;

			// calculate center, orientation, and major axis
			x = ( float )( edge[i].x + edge[j].x ) / 2.0f;
			y = ( float )( edge[i].y + edge[j].y ) / 2.0f;
			a = sqrt( ( float )SQR( edge[i].x - edge[j].x ) + ( float )SQR( edge[i].y - edge[j].y ) ) / 2.0f;
			if( 2.0f * a < min ) continue;						// length to small
			temp = ( float )( edge[i].x - edge[j].x );
			o = ( fabs( temp ) < 0.0001f ) ? 0.5f * M_PI : atan( ( float )( edge[i].y - edge[j].y ) / temp );
	
			// estimate minor axis
			for( k = 0; k < n; k++ ){					// for every
				if( k == i || k == j ) continue;		// third edge pixel
				b = GetMinor( edge[i], edge[j], edge[k], x, y, a, min );
				if( b <= 1.0f || (int)b >= na ) continue;
				// voting
				accumulator[ (int)b ]++;
			}

			// find maxium element in accumulator array
			for( k = l = 0; k < na; k++ ){
				if( accumulator[k] > max ){		// detect a ellipse
					max = accumulator[k];
					*x0 = x;
					*y0 = y;
					*a0 = a;
					*b0 = (float)k;
					*angle = o;
					l = 1;
				}
			}
			/*
			if( l != 0 ) continue;		// eliminate pixel only if ellipse is detected
			// eliminate all pixel belong this ellipse
			for( k = 0; k < n; k++ ){
				if( k == i || k == j ) continue;
				// eliminate 
				if( (int)GetMinor( edge[i], edge[j], edge[k], x, y, a, min ) == (int)*b0 ){
					for( l = k; l < n - 1; l++ ) {
						edge[l] = edge[l + 1];
						if( l + 1 == i ) i--;
						if( l + 1 == j ) j--;
					}
					k--;
					n--;
				}
			}
			/**/
//		}
	}
	// free memory
	free( edge );
	free( accumulator );
}

//
// K-Mean
//

// macor

#define InitRGB	\
	for( i = 0; i < h; i++ ){\
		for( j = 0; j < w; j++ ){\
			data[( i * w + j ) * 5 + 0] = ( double )j / ( double )( w - 1 );	\
			data[( i * w + j ) * 5 + 1] = ( double )i / ( double )( h - 1 );	\
			data[( i * w + j ) * 5 + 2] = ( double )bits[i * stride + j * 3 + 0] / 255.0;	\
			data[( i * w + j ) * 5 + 3] = ( double )bits[i * stride + j * 3 + 1] / 255.0;	\
			data[( i * w + j ) * 5 + 4] = ( double )bits[i * stride + j * 3 + 2] / 255.0;	\
		}\
	}\
	for( i = 0; i < c; i++ ){\
		center[i * 5 + 0] = ( double )cen[i * 5 + 0] / ( double )( w - 1 );\
		center[i * 5 + 1] = ( double )cen[i * 5 + 1] / ( double )( h - 1 );\
		center[i * 5 + 2] = ( double )cen[i * 5 + 2] / 255.0;\
		center[i * 5 + 3] = ( double )cen[i * 5 + 3] / 255.0;\
		center[i * 5 + 4] = ( double )cen[i * 5 + 4] / 255.0;\
	}
#define InitGray \
	for( i = 0; i < h; i++ ){	\
		for( j = 0; j < w; j++ ){	\
			data[( i * w + j ) * 3 + 0] = ( double )j / ( double )( w - 1 );	\
			data[( i * w + j ) * 3 + 1] = ( double )i / ( double )( h - 1 );	\
			data[( i * w + j ) * 3 + 2] = ( double )bits[i * stride + j * 3] / 255.0;	\
		}\
	}\
	for( i = 0; i < c; i++ ){\
		center[i * 3 + 0] = ( double )cen[i * 3 + 0] / ( w - 1 );\
		center[i * 3 + 1] = ( double )cen[i * 3 + 1] / ( h - 1 );\
		center[i * 3 + 2] = ( double )cen[i * 3 + 2] / 255.0;\
	}

// functions

void KmeanSegment( BITMAP *img, int cen[], int c, int color ){
	int *cluster, *count, i, j, max,
		w = img->bmWidth, h = img->bmHeight, stride = img->bmWidthBytes;
	double domain[10] = { 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0 },
		*data, *center;
	void *memory;
	BYTE *bits = ( BYTE* )img->bmBits;

	// memory allocate
	memory = malloc( ( h * w * 5 + c * 5 ) * sizeof( double ) + ( c + h * w ) * sizeof( int ) );	
	data = ( double* )memory;				// h * w * 5
	center = data + h * w * 5;				// c * 5
	cluster = ( int* )( center + c * 5 );	// h * w
	count = cluster + h * w;

	switch( color ){
	case RGB_COLOR:
		InitRGB;
		// k-mean algorithm
		Kmean( data, center, domain, cluster, h * w, 5, c, false );
		break;

	case GRAY_COLOR:
		InitGray;
		// k-mean algorithm
		Kmean( data, center, domain, cluster, h * w, 3, c, false );
		break;
	}
	// segment
	for( i = 0; i < c; i++ ) count[i] = 0;
	for( i = 0; i < h * w; i++ ) count[ cluster[i] ]++;
	max = 0;
	for( i = 1; i < c; i++ ) if( count[i] > count[max] ) max = i;
	for( i = 0; i < h; i++ ){
		for( j = 0; j < w; j++ ){
			if( cluster[i * w + j] != max )
				bits[i * stride + j * 3 + 0] = bits[i * stride + j * 3 + 1] = bits[i * stride + j * 3 + 2] = 0x00;
		}
	}
	free( memory );
}

void KmeanSegment( BITMAP *img, int c, int color ){
	int *center, i, w = img->bmWidth, h = img->bmHeight;
	center = ( int* )malloc( 5 * c * sizeof( int ) );
	// randomize cluster center
	switch( color ){
	case RGB_COLOR:
		for( i = 0; i < c; i++ ){
			center[i * 5 + 0] = rand() % w;
			center[i * 5 + 1] = rand() % h;
			center[i * 5 + 2] = rand() % 256;
			center[i * 5 + 3] = rand() % 256;
			center[i * 5 + 4] = rand() % 256;
		}
		break;
	case GRAY_COLOR:
		for( i = 0; i < c; i++ ){
			center[i * 3 + 0] = rand() % w;
			center[i * 3 + 1] = rand() % h;
			center[i * 3 + 2] = rand() % 256;
		}
		break;
	}
	KmeanSegment( img, center, c, color );	// kmean
	free( center );
}

void KmeanSegment( BITMAP *img, int cen[], int c, int x, int y, int color ){
	int *cluster, *count, i, j, max,
		w = img->bmWidth, h = img->bmHeight, stride = img->bmWidthBytes;
	double domain[10] = { 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0 },
		*data, *center;
	void *memory;
	BYTE *bits = ( BYTE* )img->bmBits;

	// memory allocate
	memory = malloc( ( h * w * 5 + c * 5 ) * sizeof( double ) + ( c + h * w ) * sizeof( int ) );	
	data = ( double* )memory;				// h * w * 5
	center = data + h * w * 5;				// c * 5
	cluster = ( int* )( center + c * 5 );	// h * w
	count = cluster + h * w;

	switch( color ){
	case RGB_COLOR:
		InitRGB;
		// k-mean algorithm
		Kmean( data, center, domain, cluster, h * w, 5, c, true );
		break;

	case GRAY_COLOR:
		InitGray;
		// k-mean algorithm
		Kmean( data, center, domain, cluster, h * w, 3, c, false );
		break;
	}
	// segment
	max = cluster[y * w + x];
	for( i = 0; i < h; i++ ){
		for( j = 0; j < w; j++ ){
			if( cluster[i * w + j] != max )
				bits[i * stride + j * 3 + 0] = bits[i * stride + j * 3 + 1] = bits[i * stride + j * 3 + 2] = 0x00;
		}
	}
	free( memory );
}

void KmeanSegment( BITMAP *img, int c, int x, int y, int color ){
	int *center, i, w = img->bmWidth, h = img->bmHeight;
	center = ( int* )malloc( 5 * c * sizeof( int ) );
	// randomize cluster center
	switch( color ){
	case RGB_COLOR:
		for( i = 0; i < c; i++ ){
			center[i * 5 + 0] = rand() % w;
			center[i * 5 + 1] = rand() % h;
			center[i * 5 + 2] = rand() % 256;
			center[i * 5 + 3] = rand() % 256;
			center[i * 5 + 4] = rand() % 256;
		}
		break;
	case GRAY_COLOR:
		for( i = 0; i < c; i++ ){
			center[i * 3 + 0] = rand() % w;
			center[i * 3 + 1] = rand() % h;
			center[i * 3 + 2] = rand() % 256;
		}
		break;
	}
	KmeanSegment( img, center, c, x, y, color );	// kmean
	free( center );
}

void KmeanClustering( BITMAP *img, int cen[], int c, int color ){
	int *cluster, i, j,
		w = img->bmWidth, h = img->bmHeight, stride = img->bmWidthBytes;
	double domain[10] = { 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0 },
		*data, *center;
	void *memory;
	BYTE *bits = ( BYTE* )img->bmBits;

	// memory allocate
	memory = malloc( ( h * w * 5 + c * 5 ) * sizeof( double ) + h * w * sizeof( int ) );	
	data = ( double* )memory;				// h * w * 5
	center = data + h * w * 5;				// c * 5
	cluster = ( int* )( center + c * 5 );	// h * w

	switch( color ){
	case RGB_COLOR:
		InitRGB;
		// k-mean algorithm
		Kmean( data, center, domain, cluster, h * w, 5, c, false );
		// clustering
		for( i = 0; i < h; i++ ){
			for( j = 0; j < w; j++ ){
				bits[i * stride + j * 3 + 0] = ( BYTE )( center[ cluster[i * w + j] * 5 + 2 ] * 255.0 );
				bits[i * stride + j * 3 + 1] = ( BYTE )( center[ cluster[i * w + j] * 5 + 3 ] * 255.0 );
				bits[i * stride + j * 3 + 2] = ( BYTE )( center[ cluster[i * w + j] * 5 + 4 ] * 255.0 );
			}
		}
		break;
	case GRAY_COLOR:
		InitGray;
		// k-mean algorithm
		Kmean( data, center, domain, cluster, h * w, 3, c, false );
		// clustering
		for( i = 0; i < h; i++ ){
			for( j = 0; j < w; j++ ){
				bits[i * stride + j * 3 + 0] = 
				bits[i * stride + j * 3 + 1] = 
				bits[i * stride + j * 3 + 2] = ( BYTE )( center[ cluster[i * w + j] * 3 + 2 ] * 255.0 );
			}
		}
		break;
	}
	free( memory );
}

void KmeanClustering( BITMAP *img, int c, int color ){
	int *center, i, w = img->bmWidth, h = img->bmHeight;
	center = ( int* )malloc( 5 * c * sizeof( int ) );
	// randomize cluster center
	switch( color ){
	case RGB_COLOR:
		for( i = 0; i < c; i++ ){
			center[i * 5 + 0] = rand() % w;
			center[i * 5 + 1] = rand() % h;
			center[i * 5 + 2] = rand() % 256;
			center[i * 5 + 3] = rand() % 256;
			center[i * 5 + 4] = rand() % 256;
		}
		break;
	case GRAY_COLOR:
		for( i = 0; i < c; i++ ){
			center[i * 3 + 0] = rand() % w;
			center[i * 3 + 1] = rand() % h;
			center[i * 3 + 2] = rand() % 256;
		}
		break;
	}
	KmeanClustering( img, center, c, color );	// kmean
	free( center );
}

//
// Mean Shift
//

inline double Epanechnikov_g( double x[], int n ){
	double temp = 0.0;
	for( int i = 0; i < n; i++ ) temp += x[i] * x[i];
	return ( temp <= 1.0 ) ? 1.0 : 0.0;
}

inline double Normal_g( double x[], int n ){
	double temp = 0.0;
	for( int i = 0; i < n; i++ ) temp += x[i] * x[i];
	return 0.5 * exp( -0.5 * temp );
}

void MeanShift_Color( BITMAP &img, int ***output, int *n, double bw, int param ){
	int i, j, k, x, y;
	int h = img.bmHeight;
	int w = img.bmWidth;
	double ***c;				// cluster center vector c1, c2, .... , cn
	double m[3];				// translate vector m1, m2, ... , mn
	double buf[3], xi[3], sum, temp;
	double (*kernal)( double[], int );

	switch( param ){
	case EPANECHNIKOV_KERNAL: 
		kernal = &Epanechnikov_g; 
		break;
	case NORMAL_KERNAL: 
		kernal = &Normal_g; 
		break;
	}
	
	// initialize output matrix
	*output = new int* [h];
	for( i = 0; i < h; i++ ) {
		(*output)[i] = new int [w];
		for( j = 0; j < w; j++ ) (*output)[i][j] = -1;		// default value
	}
	// initialize cluster vector
	c = new double** [h];
	for( i = 0; i < h; i++ ) {
		c[i] = new double* [w];
		for( j = 0; j < w; j++ ) c[i][j] = new double [3];
	}

	for( x = 0; x < w; x++ ){
		for( y = 0; y < h; y++ ){
			// start mean-shift for I(x, y)
			c[y][x][0] = 
				( param & LUV_COLOR ) ? ( Pixel( img, x, y, L ) / 100.0 ) : (double)Pixel( img, x, y, B ) / 255.0;
			c[y][x][1] = 
				( param & LUV_COLOR ) ? ( Pixel( img, x, y, U ) / 200.0 + 0.5 ) : (double)Pixel( img, x, y, G ) / 255.0;
			c[y][x][2] =
				( param & LUV_COLOR ) ? ( Pixel( img, x, y, V ) / 200.0 + 0.5 ) : (double)Pixel( img, x, y, R ) / 255.0;
			while( 1 ){
				if( !( param & UNIFORM_KERNAL ) ){
					sum = 0.0;
					m[0] = m[1] = m[2] = 0.0;
					for( i = 0; i < h; i++ ){
						for( j = 0; j < w; j++ ){
							// generate data point x
							xi[0] = 
								( param & LUV_COLOR ) ? ( Pixel( img, j, i, L ) / 100.0 ) : (double)Pixel( img, j, i, B ) / 255.0;
							xi[1] = 
								( param & LUV_COLOR ) ? ( Pixel( img, j, i, U ) / 200.0 + 0.5 ) : (double)Pixel( img, j, i, G ) / 255.0;
							xi[2] =
								( param & LUV_COLOR ) ? ( Pixel( img, j, i, V ) / 200.0 + 0.5 ) : (double)Pixel( img, j, i, R ) / 255.0;
							// generate input of kernal
							for( k = 0; k < 3; k++ ) buf[k] = ( xi[k] - c[y][x][k] ) / bw;
							// calculate the kernal
							temp = (*kernal)( buf, 3 );
							//  xi * g 
							for( k = 0; k < 3; k++ ) m[k] += xi[k] * temp;
							sum += temp;
						}
					}
					for( k = 0; k < 3; k++ ) m[k] = m[k] / sum - c[y][x][k];
					temp = 0.0;
					for( k = 0; k < 3; k++ ) temp = m[k] * m[k];
				}
				else{
					temp = 0.0;
					for( k = 0; k < 3; k++ ) temp = c[y][x][k] * c[y][x][k];
				}
				if( temp < ERR_BOUND ) break;
				else for( k = 0; k < 3; k++ ) c[y][x][k] += m[k];
			}
		}
	}

	// clustering
	*n = 0;
	for( sum = h * w; sum > 0.0; (*n)++ ){
		x = 0; y = 0;
		// select a cluster
		for( i = 0; i < h; i++ ){
			for( j = 0; j < w; j++ ){
				if( (*output)[i][j] == -1 ){		// if dosn't clustering
					x = j; y = i;					// mark this data point
					(*output)[i][j] = *n;
					sum--;
					break;
				}
			}
			if( j != w ) break;
		}
		// clustering
		for( i = 0; i < h; i++ ){
			for( j = 0; j < w; j++ ){
				if( (*output)[i][j] == -1 ){		// if this data point dosn't clustered yet
					if( fabs( c[i][j][0] - c[y][x][0] ) < ERR_BOUND && 
						fabs( c[i][j][1] - c[y][x][1] ) < ERR_BOUND &&  
						fabs( c[i][j][2] - c[y][x][2] ) < ERR_BOUND ) {
						(*output)[i][j] = *n;
						sum--;
					}
				}
			}
		}
	}
	for( i = 0; i < h; i++ ) {
		for( j = 0; j < w; j++ ) delete c[i][j];
		delete c[i];
	}
	delete c;
	return;
}

void MeanShift_5D( BITMAP &img, int ***output, int *n, double bw, int param ){
	int i, j, k, x, y;
	int h = img.bmHeight;
	int w = img.bmWidth;
	double ***c;				// cluster center vector c1, c2, .... , cn
	double m[5];				// translate vector m1, m2, ... , mn
	double buf[5], xi[5], sum, temp;
	double (*kernal)( double[], int );

	switch( param ){
	case EPANECHNIKOV_KERNAL: 
		kernal = &Epanechnikov_g; 
		break;
	case NORMAL_KERNAL: 
		kernal = &Normal_g; 
		break;
	}
	
	// initialize output matrix
	*output = new int* [h];
	for( i = 0; i < h; i++ ) {
		(*output)[i] = new int [w];
		for( j = 0; j < w; j++ ) (*output)[i][j] = -1;		// default value
	}
	// initialize cluster vector
	c = new double** [h];
	for( i = 0; i < h; i++ ) {
		c[i] = new double* [w];
		for( j = 0; j < w; j++ ) c[i][j] = new double [5];
	}

	for( x = 0; x < w; x++ ){
		for( y = 0; y < h; y++ ){
			// start mean-shift for I(x, y)
			c[y][x][0] = 
				( param & LUV_COLOR ) ? ( Pixel( img, x, y, L ) / 100.0 ) : (double)Pixel( img, x, y, B ) / 255.0;
			c[y][x][1] = 
				( param & LUV_COLOR ) ? ( Pixel( img, x, y, U ) / 200.0 + 0.5 ) : (double)Pixel( img, x, y, G ) / 255.0;
			c[y][x][2] =
				( param & LUV_COLOR ) ? ( Pixel( img, x, y, V ) / 200.0 + 0.5 ) : (double)Pixel( img, x, y, R ) / 255.0;
			c[y][x][3] = (double)x / w;	c[y][x][4] = (double)y / h;
			while( 1 ){
				if( !( param & UNIFORM_KERNAL ) ){
					sum = 0.0;
					m[0] = m[1] = m[2] = m[3] = m[4] = 0.0;
					for( i = 0; i < h; i++ ){
						for( j = 0; j < w; j++ ){
							// generate data point x
							xi[0] = 
								( param & LUV_COLOR ) ? ( Pixel( img, j, i, L ) / 100.0 ) : (double)Pixel( img, j, i, B ) / 255.0;
							xi[1] = 
								( param & LUV_COLOR ) ? ( Pixel( img, j, i, U ) / 200.0 + 0.5 ) : (double)Pixel( img, j, i, G ) / 255.0;
							xi[2] =
								( param & LUV_COLOR ) ? ( Pixel( img, j, i, V ) / 200.0 + 0.5 ) : (double)Pixel( img, j, i, R ) / 255.0;
							xi[3] = (double)j / w; xi[4] = (double)i / h;
							// generate input of kernal
							for( k = 0; k < 5; k++ ) buf[k] = ( xi[k] - c[y][x][k] ) / bw;
							// calculate the kernal
							temp = Epanechnikov_g( buf, 5 );
							//  xi * g 
							for( k = 0; k < 5; k++ ) m[k] += xi[k] * temp;
							sum += temp;
						}
					}
					for( k = 0; k < 5; k++ ) m[k] = m[k] / sum - c[y][x][k];
					temp = 0.0;
					for( k = 0; k < 5; k++ ) temp = m[k] * m[k];
				}
				else{
					temp = 0.0;
					for( k = 0; k < 5; k++ ) temp = c[y][x][k] * c[y][x][k];
				}
				if( temp < ERR_BOUND ) break;
				else for( k = 0; k < 5; k++ ) c[y][x][k] += m[k];
			}
		}
	}

	// clustering
	*n = 0;
	for( sum = h * w; sum > 0.0; (*n)++ ){
		x = 0; y = 0;
		// select a cluster
		for( i = 0; i < h; i++ ){
			for( j = 0; j < w; j++ ){
				if( (*output)[i][j] == -1 ){		// if dosn't clustering
					x = j; y = i;					// mark this data point
					(*output)[i][j] = *n;
					sum--;
					break;
				}
			}
			if( j != w ) break;
		}
		// clustering
		for( i = 0; i < h; i++ ){
			for( j = 0; j < w; j++ ){
				if( (*output)[i][j] == -1 ){		// if this data point dosn't clustered yet
					if( fabs( c[i][j][0] - c[y][x][0] ) < ERR_BOUND && 
						fabs( c[i][j][1] - c[y][x][1] ) < ERR_BOUND &&  
						fabs( c[i][j][2] - c[y][x][2] ) < ERR_BOUND &&  
						fabs( c[i][j][3] - c[y][x][3] ) < ERR_BOUND &&  
						fabs( c[i][j][4] - c[y][x][4] ) < ERR_BOUND  ) {
						(*output)[i][j] = *n;
						sum--;
					}
				}
			}
		}
	}
	for( i = 0; i < h; i++ ) {
		for( j = 0; j < w; j++ ) delete c[i][j];
		delete c[i];
	}
	delete c;
	return;
}

void MeanShift( BITMAP &img, int ***output, int *n, double bw, int space ){
	if( space & SPACE ) MeanShift_5D( img, output, n, bw, space );
	else MeanShift_Color( img, output, n, bw, space );
}
