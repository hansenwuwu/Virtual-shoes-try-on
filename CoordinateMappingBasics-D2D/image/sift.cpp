
#ifdef DONT_COMPILE_THIS

#include "bitmap.h"
#include "imgproc.h"
#include "D:\Rice\workspace\matrixLib\matrix.h"
#include <math.h>
#include <stdio.h>
#include <vector>
#include <stdexcept>
using namespace std;

#define M_PI 3.14159265358979323846
#define bit( b, x, y ) ( b[(y) * stride + (x)] )
#define EXTREMA_THRES 0.01
#define CONTRAST_TRHES 0.03
#define THRESHOLD 10
#define Tr( h ) ( (h)[0][0] + (h)[1][1] )
#define Det( h ) ( (h)[0][0] * (h)[1][1] - (h)[0][1] * (h)[1][0] )
#define sqr( x ) ( (x)*(x) )
#define multi( a, b ) (a).x *= (b); (a).y *= (b);
#define CHECK_EXTREMA( buf, k, x, y, comp, thres ) \
	(	\
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[((k)-1)][((y)-1) * stride + ((x)-1) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[((k)-1)][((y)-1) * stride + (x) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[((k)-1)][((y)-1) * stride + ((x)+1) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[((k)-1)][(y) * stride + ((x)-1) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[((k)-1)][(y) * stride + (x) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[((k)-1)][(y) * stride + ((x)+1) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[((k)-1)][((y)+1) * stride + ((x)-1) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[((k)-1)][((y)+1) * stride + (x) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[((k)-1)][((y)+1) * stride + ((x)+1) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[(k)][((y)-1) * stride + ((x)-1) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[(k)][((y)-1) * stride + (x) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[(k)][((y)-1) * stride + ((x)+1) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[(k)][(y) * stride + ((x)-1) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[(k)][(y) * stride + ((x)+1) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[(k)][((y)+1) * stride + ((x)-1) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[(k)][((y)+1) * stride + (x) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[(k)][((y)+1) * stride + ((x)+1) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[((k)+1)][((y)-1) * stride + ((x)-1) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[((k)+1)][((y)-1) * stride + (x) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[((k)+1)][((y)-1) * stride + ((x)+1) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[((k)+1)][(y) * stride + ((x)-1) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[((k)+1)][(y) * stride + (x) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[((k)+1)][(y) * stride + ((x)+1) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[((k)+1)][((y)+1) * stride + ((x)-1) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[((k)+1)][((y)+1) * stride + (x) ] && \
	buf[(k)][(y) * stride + (x)] + (thres) comp buf[((k)+1)][((y)+1) * stride + ((x)+1) ] \
	)

//
// Structures
//

typedef struct _Descriptor{
	double histogram[128];
} Descriptor;

typedef struct _Gradient{
	double x, y;
} Gradient;

typedef struct _Keypoint{
	double x;
	double y;		// x, y position
	double scale;
	int octave;
	double orien;
	double magni;
	Descriptor desc;
} Keypoint;

//
// inline functions
//

inline double Gaussian( double sigma2, double x, double y ){
	return 0.5 / ( M_PI * sigma2 ) *
		exp( -0.5 * ( x * x + y * y ) / sigma2 ); 
}

inline double LoG( double sigma2, double x2, double y2 ){
	return ( x2 + y2 - 2 * sigma2 ) / sqr( sigma2 ) *
		exp( -0.5 * ( x2 + y2 ) / sigma2 );
}

inline double DoG( double sigma, double x, double y, double k ){
	return ( k - 1.0 ) * sqr( sigma ) * LoG( sqr( sigma ), sqr( x ), sqr( y ) );
}

inline int Pixel( BITMAP &bmp, int x, int y, int rgb ){
	BYTE *bits = ( BYTE* )bmp.bmBits;
	return bits[ y * bmp.bmWidthBytes + x * 3 + rgb ];
}

inline double Pixel( BITMAP &bmp, int x, int y ){
	double *bits = ( double* )bmp.bmBits;
	return bits[ y * bmp.bmWidthBytes + x ];
}

inline double Pixel( double *bits, int stride, double x, double y ){
	int _x = (int)x;
	int _y = (int)y;
	double u = x - _x;
	double v = y - _y;
	return 
		( 1.0 - u ) * ( 1.0 - v ) * bit( bits, _x, _y ) +
		( u ) * ( 1.0 - v ) * bit( bits, _x + 1, _y ) + 
		( 1.0 - u ) * ( v ) * bit( bits, _x, _y + 1 ) +
		u * v * bit( bits, _x + 1, _y + 1 );
}

inline int pow2( int x ){
	int i, ans = 1;
	for( i = 0; i < x; i++ ) ans *= 2;
	return ans;
}

inline int angle90( double x, double y ){
	if( x == 0 && y > 0 ) return 90;	// +y
	else if( x == 0 && y < 0 ) return 270;	// -y
	else if( x > 0 && y == 0 ) return 0;	// +x
	else if( x < 0 && y == 0 ) return 180;	// -x
	else if( x > 0 && y > 0 ) return (int)( atan( y / x ) * 180.0 / M_PI );	// 1
	else if( x < 0 && y > 0 ) return (int)( atan( y / x ) * 180.0 / M_PI ) + 180;	// 2
	else if( x < 0 && y < 0 ) return (int)( atan( y / x ) * 180.0 / M_PI ) + 180;	// 3
	else if( x > 0 && y < 0 ) return (int)( atan( y / x ) * 180.0 / M_PI ) + 360;	// 4
	else return 0;
}

// Convolution

void convolution(double* &img, int h, int w, int stride, int size, double *kernel, int x_low, int y_low, int x_up, int y_up ){
	double* output=(double*)malloc(h*stride*sizeof(double));
	double sum;
	int i, j, x, y;

	for(y = y_low; y < y_up; y++){		
		for(x = x_low; x < x_up; x++){
			if(y<=(size-1)/2-1 || y>=h-(size-1)/2 || x<=(size-1)/2-1 || x>=w-(size-1)/2) {
				output[y*stride+x]=0.0;
				continue;
			}
			else{
				for(i=0, sum=0; i<size; i++){
					for(j=0; j<size; j++){
						sum+=kernel[i*size+j]*img[(y+i-(size-1)/2)*stride+(x+j-(size-1)/2)];
					}
				}
				output[y*stride+x]=sum;
			}
		}
	}
	free(img);
	img=output;
}

void convolution(double* &img, int h, int w, int stride, int size, double *kernel){
	convolution( img, h, w, stride, size, kernel, 0, 0, w, h );
}

//
// SIFT
//

//
// SIFT per octave
//

void sift_octave( BITMAP &img, int oct, int s, double sig, vector<Keypoint> &key, vector<Descriptor> &des, int *pn ){
	int &n = *pn;
	int x, y, xx, yy, i, it;
	int h = img.bmHeight;
	int w = img.bmWidth;
	int stride, size;
	double dx, dy, ds, k, sigma;
	double **bits, dog[400], H[2][2];
	vec D1( 3 );
	mat D2( 3, 3 );
	
	// initialize
	bits = ( double** )calloc( ( s + 2), sizeof( double* ) );
	stride = ( img.bmWidth % 4 == 0 ) ? img.bmWidth : ( img.bmWidth + 4 - ( img.bmWidth % 4 ) );

	// copy a 8-bytes bitmap
	for( i = 0; i < s + 2; i++ ){
		bits[i] = ( double* )calloc( stride * h, sizeof(double) );	
		for( y = 0; y < h; y++ ){
			for( x = 0; x < w; x++ )
				bit( bits[i], x, y ) = Pixel( img, x, y );
		}
	}

	// construct dog
	k = pow( 2.0, 1.0 / (double)s );
	sigma = sig / k;			// sigma/k, 1, 1 sigma, 2 sigma, ... , s sigma
	for( i = 0; i < s + 2; i++, sigma *= k ){
		size = (int)( 6.0 * sigma + 1.5 );
		if( size % 2 == 0 ) size -= 1;
		for( x = 0; x < size; x++ ){
			for( y = x; y < size - x; y++ ) {
				dog[ y * size + x ] = 
				dog[ x * size + y ] =
				dog[ ( size - 1 - y ) * size + ( size - 1 - x ) ] =
				dog[ ( size - 1 - x ) * size + ( size - 1 - y ) ] =
				DoG( sigma, (double)x - 0.5 * ( size - 1 ), (double)y - 0.5 * ( size - 1 ), k );
			}
		}
		convolution( bits[i], h, w, stride, size, dog );
	}

	// find local extrema
	sigma = sig;
	n = 0;
	size = key.size();
	for( i = 1; i < s + 1; i++, sigma *= k ){
		for( y = 16; y < h - 16; y++ ){
			for( x = 16; x < w - 16; x++ ){
				// Check local extrema
				if( !CHECK_EXTREMA( bits, i, x, y, <, EXTREMA_THRES ) ){
					if( !CHECK_EXTREMA( bits, i, x, y, >, -EXTREMA_THRES ) ) continue;
				}

				// Accurate keypoint localization
				dx = dy = 0.0;
				xx = x; yy = y;
				for( it = 0; ( dx == 0.0 || dy == 0.0 ) && it < 5; it++ ){
					// D1
					D1[0] = 0.5 * ( bit( bits[i], xx + 1, yy ) - bit( bits[i], xx - 1, yy ) );
					D1[1] = 0.5 * ( bit( bits[i], xx, yy + 1 ) - bit( bits[i], xx, yy - 1 ) );
					D1[2] = 0.5 * ( bit( bits[i + 1], xx, yy ) - bit( bits[i - 1], xx, yy ) );
					
					// D2
					D2[0][0] = H[0][0] = 
						bit( bits[i], xx + 1, yy ) + bit( bits[i], xx - 1, yy ) - 2.0 * bit( bits[i], xx, yy );
					D2[1][1] = H[1][1] = 
						bit( bits[i], xx, yy + 1 ) + bit( bits[i], xx, yy - 1 ) - 2.0 * bit( bits[i], xx, yy );
					D2[2][2] = 
						bit( bits[i + 1], xx, yy ) + bit( bits[i - 1], xx, yy ) - 2.0 * bit( bits[i], xx, yy );
					D2[0][1] = D2[1][0] = H[0][1] = H[1][0] = 0.25 * (
						bit( bits[i], xx + 1, yy + 1 ) + bit( bits[i], xx - 1, yy - 1 ) -
						bit( bits[i], xx + 1, yy - 1 ) - bit( bits[i], xx - 1, yy + 1 ) );
					D2[0][2] = D2[2][0] = 0.25 * (
						bit( bits[i + 1], xx + 1, yy ) + bit( bits[i - 1], xx - 1, yy ) -
						bit( bits[i - 1], xx + 1, yy ) - bit( bits[i + 1], xx - 1, yy ) );
					D2[1][2] = D2[2][1] = 0.25 * (
						bit( bits[i + 1], xx, yy + 1 ) + bit( bits[i - 1], xx, yy - 1 ) -
						bit( bits[i - 1], xx, yy + 1 ) - bit( bits[i + 1], xx, yy - 1 ) );
					
					try{
						D1 = -1.0 * D2.inverse() * D1;
						dx = D1[0];
						dy = D1[1];
						ds = D1[2];
					}catch( ... ){
						it = 5;
						break;
					}

					if( ds < -0.5 || ds > 0.5 ) { ds = 0.0; break; }		// other octave is better

					// check x
					if( dx > 0.5 && xx < w - 2 ) { dx = 0.0; xx++; }
					else if( dx < -0.5 && xx > 1 ) { dx = 0.0; xx--; }

					// check y
					if( dy > 0.5 && yy < h - 2 ) { dy = 0.0; yy++; }
					else if( dy < -0.5 && yy > 1 ) { dy = 0.0; yy--; }
				}

				if( ds == 0.0 ) continue;
				if( it == 5 ) continue;

				// Check minimum contrast
				if( fabs( bit( bits[i], xx, yy ) + 0.5 * ( D1[0] * dx + D1[1] * dy + D1[2] * ds ) ) < CONTRAST_TRHES ) continue;
				dx += (double)xx;
				dy += (double)yy;

				// Check Edge Effect
				if( sqr( Tr(H) ) / ( Det(H) ) > THRESHOLD ) continue;

				// Save Keypoints
				key.push_back( Keypoint() );
				key[n + size].x = dx; 
				key[n + size].y = dy; 
				key[n + size].scale = sigma;
				key[n + size].octave = oct;
				n++;
			}
		}
	}

	// free memory	
	for( i = 0; i < s + 2; i++ ) free( bits[i] );
	free( bits );
}

//
// Compute Descriptor
//

#define DegToArc( x ) ( 0.01745329251994329576923690768489 * (x) )
#define Rotate( u, v ) Pixel( buf, stride, (u) * cosine - (v) * sine + key.x, (u) * sine + (v) * cosine + key.y )
// Set Width of descriptor size = 4 * 3 * scale = 12 * scale
// sigma = 0.5 * width = 6 * scale
void ComputeDesc( double buf[], int stride, Keypoint &key ){
	//		| cos	-sin |
	// R =	|			 |
	//		| sin	 cos |
/*	double cosine = cos( DegToArc( key.orien ) );
	double sine = sin( DegToArc( key.orien ) );
	double w[8][8];
	double du, dv, grad_x, grad_y, bin[8];
	Descriptor desc;
	int i, j, k, offset;

	// initial weighting table, w[0][0] is closest to origin
	for( i = 0; i < 12; i++ ){
		for( j = i; j < 12; j++ ){
			w[i][j] = w[j][i] = Gaussian( 6 * sqr( key.scale ), 0.5 + (double)i, 0.5 + (double)j );
		}
	}

	// gradient
	for( i = 0; i < 12 * scale; i += 4 ){
		for( j = 0; j < 12 * scale; j += 4 ){
			for( k = 0; k < 8; k++ ) bin[k] = 0.0;	// initialize the bin of orientation
			for( du = (double)i - 7.5; du < (double)i - 4.0; du += 1.0 ){
				for( dv = (double)j - 7.5; dv < (double)j - 4.0; dv += 1.0 ){
					grad_x = Rotate( du + 1.0, dv ) - Rotate( du, dv );
					grad_y = Rotate( du, dv + 1.0 ) - Rotate( du, dv );
					bin[ angle90( grad_x, grad_y ) / 45 ] +=
						w[ abs( (int)du ) ][ abs( (int)dv ) ] * sqrt( sqr( grad_x ) + sqr( grad_y ) );
				}
			}
			offset = i * 8 + j * 2;
			for( k = 0; k < 8; k++ ) desc.histogram[ offset + k ] = bin[k];
		}
	}

	// Normalize
	bin[0] = 0.0;
	for( i = 0; i < 128; i++ ) bin[0] += sqr( desc.histogram[i] );
	bin[0] = sqrt( bin[0] );
	for( i = offset = 0; i < 128; i++ ) {
		desc.histogram[i] /= bin[0];
		if( desc.histogram[i] > 0.2 ){	// handle the element > 0.2
			desc.histogram[i] = 0.2;
			offset++;
		}
	}

	// Re-Normalize
	if( offset > 0 ){
		bin[0] = 0.0;
		for( i = 0; i < 128; i++ ) bin[0] += sqr( desc.histogram[i] );
		bin[0] = sqrt( bin[0] );
		for( i = 0; i < 128; i++ ) desc.histogram[i] /= bin[0];
	}

	key.desc = desc;*/
}

//
// Compute Orientation and magnitude per keypoint
//

void ComputeOrien( BITMAP img[], vector<Keypoint> &key, int n ){
	int i, j;
	int stride = img[ key[n].octave ].bmWidthBytes;
	int x = (int)( key[n].x + 0.5 );
	int y = (int)( key[n].y + 0.5 );
	int r, size = 6 * key[n].scale + 1.5;
	int x_low = max( x - size, 0 );
	int x_up = min( x + size, img[ key[n].octave ].bmWidth );
	int y_low = max( y - size, 0 );
	int y_up = min( y + size, img[ key[n].octave ].bmHeight );
	int win_x_low, win_x_up, win_y_low, win_y_up, theta;
	double bin[36] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	double grad_x, grad_y, kernal[400];
	double *bits = ( double* )img[ key[n].octave ].bmBits;
	double *buf = ( double* )calloc( img[ key[n].octave ].bmHeight * img[ key[n].octave ].bmWidthBytes, sizeof( double ) );
	int flag;

	if( size % 2 == 0 ) size--;
	r = size / 2;

	// copy image
	win_x_low = max( x_low - r, 0 );
	win_x_up = min( x_up + r, img[ key[n].octave ].bmWidth );
	win_y_low = max( y_low - r, 0 );
	win_y_up = min( y_up + r, img[ key[n].octave ].bmHeight );
	for( i = win_x_low; i < win_x_up; i++ ){
		for( j = win_y_low; j < win_y_up; j++ ) bit( buf, i, j ) = bit( bits, i, j );
	}
	
	// construct gaussian kernal and smooth image
	for( i = 0; i < size; i++ ){
		for( j = i; j < size - i; j++ ){
			kernal[ j * size + i ] = 
			kernal[ i * size + j ] = 
			kernal[ ( size - 1 - j ) * size + ( size - 1 - i ) ] =
			kernal[ ( size - 1 - i ) * size + ( size - 1 - j ) ] =
			Gaussian( sqr( key[n].scale ), (double)i - 0.5 * ( size - 1 ), (double)j - 0.5 * ( size - 1 ) );
		}
	}
	convolution( buf, img[ key[n].octave ].bmHeight, img[ key[n].octave ].bmWidth, 
		stride, size, kernal, x_low, y_low, x_up, y_up );

	// voting for orientation
	for( i = x_low; i < x_up; i++ ){
		for( j = y_low; j < y_up; j++ ){
			if( i + 1 < x_up ) grad_x = bit( buf, i + 1, j ) - bit( buf, i, j );
			else grad_x = bit( buf, i, j ) - bit( buf, i - 1, j );
			if( j + 1 < y_up ) grad_y = bit( buf, i, j + 1 ) - bit( buf, i, j );
			else grad_y = bit( buf, i, j ) - bit( buf, i, j - 1 );
			theta = angle90( grad_x, grad_y );
			bin[ theta / 36 ] += 
				Gaussian( sqr( 1.5 * key[n].scale ), i - x, j - y ) * sqrt( sqr( grad_x ) + sqr( grad_y ) );
		}
	}

	// find max orientation
	flag =0;
	for( i = 1; i < 36; i++ ) if( bin[i] > bin[flag] ) flag = i;
	key[n].orien = flag * 10;
	key[n].magni = bin[flag];

	ComputeDesc( buf, stride, key[n] );

	// find the other orientation ( > 80% )
	for( i = 0; i < 36; i++ ){
		if( i == flag ) continue;
		if( bin[i] > 0.8 * bin[flag] ){
			Keypoint newkey = key[n];
			newkey.orien = i * 10;
			newkey.magni = bin[i];
			ComputeDesc( buf, stride, newkey );
			key.push_back( newkey );
		}
	}

	free( buf );
}

//
// SIFT driver
//

void sift( BITMAP &img, int oct, int s, double sig, double **pkey, double **pdesc, int *pn ){
	vector< Keypoint > Key;
	vector< Descriptor> Desc;
	BITMAP *img8;				// image pyramid
	double *bits, *old_bits, k;
	int x, y, i, n;
	int h = img.bmHeight, w = img.bmWidth;
	int stride;

	// initialize image array
	img8 = ( BITMAP* ) calloc( oct, sizeof( BITMAP ) );
	img8[0] = img;

	// initialize
	img8[0].bmWidthBytes = ( w % 4 == 0 ) ? w : ( w + 4 - ( w % 4 ) );
	stride = img8[0].bmWidthBytes;
	img8[0].bmBits = calloc( stride * h, sizeof(double) );
	bits = old_bits = ( double* )img8[0].bmBits;

	// Setting gray images
	for( y = 0; y < h; y++ ){
		for( x = 0; x < w; x++ ){
			k = Pixel( img, x, y, 0 ) + Pixel( img, x, y, 1 ) + Pixel( img, x, y, 2 ); 
			bits[ y * stride + x ] = k / 765;
		}
	}

	n = 0;
	for( i = 0; i < oct; i++ ){
		// Detect Keypoints per octave
		sift_octave( img8[i], i, s, sig, Key, Desc, pn );
		n += *pn;
		
		if( i == oct - 1 ) break;

		// construct next octave image
		img8[i + 1].bmHeight = h = h / 2;
		img8[i + 1].bmWidth = w = w / 2;
		img8[i + 1].bmWidthBytes = ( w % 4 == 0 ) ? w : ( w + 4 - ( w % 4 ) );
		img8[i + 1].bmBits = calloc( stride * h, sizeof(double) );
		bits = ( double* )img8[i + 1].bmBits;
		
		// scaling
		for( y = 0; y < h; y++ ){
			for( x = 0; x < w; x++ ){
				bits[ y * img8[i + 1].bmWidthBytes + x ] = 0.25 * (
					old_bits[ y * 2 * stride + x * 2 ] + 
					old_bits[ ( y * 2 + 1 ) * stride + x * 2 ] + 
					old_bits[ y * 2 * stride + ( x * 2 + 1 ) ] + 
					old_bits[ ( y * 2 + 1 ) * stride + ( x * 2 + 1 ) ]
					);
//					old_bits[ y * 2 * stride + x * 2 ];
			}
		}
		stride = img8[i + 1].bmWidthBytes;
		old_bits = bits;
	}

	// generate keypoint descriptor
	for( i = 0; i < n; i++ ){
		ComputeOrien( img8, Key, i );
	}

	*pn = n = Key.size();
	IMAGE test( "im32.bmp" );
	for( x = 0; x < n; x++ ){
		for( i = 0; i <= Key[x].magni * pow2( Key[x].octave ) * 1000; i++ ){
			double dx = i * cos( M_PI / 180 * Key[x].orien );
			double dy = i * sin( M_PI / 180 * Key[x].orien );
			test.SetPixel( Key[x].x * pow2( Key[x].octave ) + dx, Key[x].y * pow2( Key[x].octave ) + dy, 0x0000FF ); 
		}
	}
	WriteBMP( "test.bmp", test.BMP() );
	printf( "testing\n" );

	*pkey = new double [2 * n];
	*pdesc = new double [128 * n];
	for( i = 0; i < n; i++ ){
		(*pkey)[i * 2 + 0] = Key[i].x * pow2( Key[i].octave );
		(*pkey)[i * 2 + 1] = Key[i].y * pow2( Key[i].octave );
		for( x = 0; x < 128; x++ ) (*pdesc)[i * 128 + x] = Key[i].desc.histogram[x];
	}

	// free memory
	for( i = 0; i < oct; i++ ) free( img8[i].bmBits );
	free( img8 );
}


#endif