//  
//  Image process
//

#include "imgproc.h"
#include <math.h>
#include "linear.h"
#include "algo.h"
#include <string.h>

//
// define and macro
//

#define M_PI	3.1415926535897932384626433832795

//
// inline and template
//

template< typename T > T SQR( T x ){
	return x * x;
}

inline float Gauss( float x, float a ){
	return ( float )exp( -0.5f * x * x / a / a ); 
}

inline float Gauss( double x, double a ){
	return exp( -0.5f * x * x / a / a ); 
}

inline void GenGaussian( float G[], int w ){
	int a = w / 2, i, j;
	float x, b = 0.5f / ( float )SQR( a ), s = b / ( float )M_PI;
	for( i = -a; i <= a; i++ ){
		for( j = i; j <= a; j++ ){
			x = ( float )( SQR( i ) + SQR( j ) );
			G[ ( a + i ) * w + a + j ] = 
			G[ ( a + j ) * w + a + i ] = ( float )( exp( -x * b ) * s );
		}
	}
}

inline void GenGaussian( double G[], int w ){
	int a = w / 2, i, j;
	double x, b = 0.5 / ( double )SQR( a ), s = b / ( double )M_PI;
	for( i = -a; i <= a; i++ ){
		for( j = i; j <= a; j++ ){
			x = ( double )( SQR( i ) + SQR( j ) );
			G[ ( a + i ) * w + a + j ] = 
			G[ ( a + j ) * w + a + i ] = exp( -x * b ) * s;
		}
	}
}

inline BYTE CheckColor( const double x ){
	if( x <= 0.5 ) return 0x00;
	else if( x > 254.5 ) return 0xff;
	else return ( BYTE )x;
}

inline BYTE CheckColor( const float x ){
	if( x <= 0.5f ) return 0x00;
	else if( x > 254.5f ) return 0xff;
	else return ( BYTE )x;
}

inline BYTE CheckColor( const int x ){
	if( x < 0 ) return 0x00;
	else if( x > 255 ) return 0xff;
	else return ( BYTE )x;
}

// function declaration

//
// average
//

void AveImage( float x[], float ave[], int col, int row ){
	float *pa = ave, *px = x;
	int i, j, w = col - 1, h = row - 1;

	// ( 0, 0 )
	*pa = *( x ) + *( x + 1 ) + *( x + col ) + *( x + col + 1 );
	pa++;	px++;

	// ( 0, 1 )
	*pa = *( pa - 1 ) + *( px + 1 ) + *( px + col + 1 );
	pa++;	px++;

	// ( 0, 2 ) ~ ( 0, col - 2 )
	for( i = 2; i < w; i++, pa++, px++ ){
		*pa = *( pa - 1 ) - *( px - 2 ) - *( px + col - 2 ) + *( px + 1 ) + *( px + col + 1 );
	}

	// ( 0, col - 1 )
	*pa = *( pa - 1 ) - *( px - 2 ) - *( px + col - 2 );
	pa++;	px++;

	// ( 1, 0 )
	*pa = *( pa - col ) + *( px + col ) + *( px + col + 1 );
	pa++;	px++;

	// ( 1, 1 )
	*pa = *( pa - 1 ) + *( px - col + 1 ) + *( px + 1 ) + *( px + col + 1 );
	pa++;	px++;

	// ( 1, 2 ) ~ ( 1, w - 2 )
	for( i = 2; i < w; i++, pa++, px++ ){
		*pa = *( pa - 1 ) - *( px - col - 2 ) - *( px - 2 ) - *( px + col - 2 ) + *( px - col + 1 ) + *( px + 1 ) + *( px + col + 1 );
	}

	// ( 1, w - 1 )
	*pa = *( pa - 1 ) - *( px - col - 2 ) - *( px - 2 ) - *( px + col - 2 );
	pa++;	px++;

	for( i = 2; i < h; i++, pa++, px++ ){
		// ( i, 0 )
		*pa = *( pa - col ) - *( px - 2 * col ) - *( px - 2 * col + 1 ) + *( px + col ) + *( px + col + 1 );
		pa++;	px++;

		// ( i, 1 )
		*pa = *( pa - 1 ) + *( px - col + 1 ) + *( px + 1 ) + *( px + col + 1 );
		pa++;	px++;

		// ( i, 2 ) ~ ( i, w - 2 )
		for( j = 2; j < w; j++, pa++, px++ ){
			*pa = *( pa - 1 ) - *( px - col - 2 ) - *( px - 2 ) - *( px + col - 2 ) + *( px - col + 1 ) + *( px + 1 ) + *( px + col + 1 );
		}

		// ( i, w - 1 )
		*pa = *( pa - 1 ) - *( px - col - 2 ) - *( px - 2 ) - *( px + col - 2 );
	}

	// ( h - 1, 0 )
	*pa = *( pa - col ) - *( px - 2 * col ) - *( px - 2 * col + 1 );
	pa++;	px++;

	// ( h - 1, 1 )
	*pa = *( pa - 1 ) + *( px + 1 ) + *( px - col + 1 );
	pa++;	px++;

	// ( h - 1, 2 ) ~ ( h - 1, w - 2 )
	for( i = 2; i < w; i++, pa++, px++ ){
		*pa = *( pa - 1 ) - *( px - col - 2 ) - *( px - 2 ) + *( px - col + 1 ) + *( px + 1 );
	}

	// ( h - 1, w - 1 )
	*pa = *( pa - 1 ) - *( px - 2 ) - *( px - col - 2 );
	
	// average
	pa = ave;
	*( pa++ ) /= 4.0;
	for( i = 1; i < w; i++, pa++ ) *pa /= 6.0;
	*( pa++ ) /= 4.0;
	for( i = 1; i < w; i++ ){
		*( pa++ ) /= 6.0;
		for( j = 1; j < h; j++, pa++ ) *pa /= 9.0;
		*( pa++ ) /= 6.0;
	}
	*( pa++ ) /= 4.0;
	for( i = 0; i < w; i++, pa++ ) *pa /= 6.0;
	*( pa++ ) /= 4.0;
}

//
// E[x^2]
//

void Ave2Image( float x[], float ave[], int col, int row ){
	float *pa = ave, *px = x;
	int i, j, w = col - 1, h = row - 1;

	// ( 0, 0 )
	*pa = SQR( *( x ) ) + SQR( *( x + 1 ) ) + SQR( *( x + col ) ) + SQR( *( x + col + 1 ) );
	pa++;	px++;

	// ( 0, 1 )
	*pa = *( pa - 1 ) + SQR( *( px + 1 ) ) + SQR( *( px + col + 1 ) );
	pa++;	px++;

	// ( 0, 2 ) ~ ( 0, col - 2 )
	for( i = 2; i < w; i++, pa++, px++ ){
		*pa = *( pa - 1 ) - SQR( *( px - 2 ) ) - SQR( *( px + col - 2 ) ) + SQR( *( px + 1 ) ) + SQR( *( px + col + 1 ) );
	}

	// ( 0, col - 1 )
	*pa = *( pa - 1 ) - SQR( *( px - 2 ) ) - SQR( *( px + col - 2 ) );
	pa++;	px++;

	// ( 1, 0 )
	*pa = *( pa - col ) + SQR( *( px + col ) ) + SQR( *( px + col + 1 ) );
	pa++;	px++;

	// ( 1, 1 )
	*pa = *( pa - 1 ) + SQR( *( px - col + 1 ) ) + SQR( *( px + 1 ) ) + SQR( *( px + col + 1 ) );
	pa++;	px++;

	// ( 1, 2 ) ~ ( 1, w - 2 )
	for( i = 2; i < w; i++, pa++, px++ ){
		*pa = *( pa - 1 ) - SQR( *( px - col - 2 ) ) - SQR( *( px - 2 ) ) - SQR( *( px + col - 2 ) )
			+ SQR( *( px - col + 1 ) ) + SQR( *( px + 1 ) ) + SQR( *( px + col + 1 ) );
	}

	// ( 1, w - 1 )
	*pa = *( pa - 1 ) - SQR( *( px - col - 2 ) ) - SQR( *( px - 2 ) ) - SQR( *( px + col - 2 ) );
	pa++;	px++;

	for( i = 2; i < h; i++, pa++, px++ ){
		// ( i, 0 )
		*pa = *( pa - col ) - SQR( *( px - 2 * col ) ) - SQR( *( px - 2 * col + 1 ) ) + SQR( *( px + col ) ) + SQR( *( px + col + 1 ) );
		pa++;	px++;

		// ( i, 1 )
		*pa = *( pa - 1 ) + SQR( *( px - col + 1 ) ) + SQR( *( px + 1 ) ) + SQR( *( px + col + 1 ) );
		pa++;	px++;

		// ( i, 2 ) ~ ( i, w - 2 )
		for( j = 2; j < w; j++, pa++, px++ ){
			*pa = *( pa - 1 ) - SQR( *( px - col - 2 ) ) - SQR( *( px - 2 ) ) - SQR( *( px + col - 2 ) ) + 
				SQR( *( px - col + 1 ) ) + SQR( *( px + 1 ) ) + SQR( *( px + col + 1 ) );
		}

		// ( i, w - 1 )
		*pa = *( pa - 1 ) - SQR( *( px - col - 2 ) ) - SQR( *( px - 2 ) ) - SQR( *( px + col - 2 ) );
	}

	// ( h - 1, 0 )
	*pa = *( pa - col ) - SQR( *( px - 2 * col ) ) - SQR( *( px - 2 * col + 1 ) );
	pa++;	px++;

	// ( h - 1, 1 )
	*pa = *( pa - 1 ) + SQR( *( px + 1 ) ) + SQR( *( px - col + 1 ) );
	pa++;	px++;

	// ( h - 1, 2 ) ~ ( h - 1, w - 2 )
	for( i = 2; i < w; i++, pa++, px++ ){
		*pa = *( pa - 1 ) - SQR( *( px - col - 2 ) ) - SQR( *( px - 2 ) ) + SQR( *( px - col + 1 ) ) + SQR( *( px + 1 ) );
	}

	// ( h - 1, w - 1 )
	*pa = *( pa - 1 ) - SQR( *( px - 2 ) ) - SQR( *( px - col - 2 ) );
	
	// average get E[x^2]
	pa = ave;
	*( pa++ ) /= 4.0;
	for( i = 1; i < w; i++, pa++ ) *pa /= 6.0;
	*( pa++ ) /= 4.0;
	for( i = 1; i < w; i++ ){
		*( pa++ ) /= 6.0;
		for( j = 1; j < h; j++, pa++ ) *pa /= 9.0;
		*( pa++ ) /= 6.0;
	}
	*( pa++ ) /= 4.0;
	for( i = 0; i < w; i++, pa++ ) *pa /= 6.0;
	*( pa++ ) /= 4.0;
}

void VarImage( float ave[], float ave2[], float var[], int col, int row ){
	float *pa = ave, *pa2 = ave2;
	int i, n = row * col;

	for( i = 0; i < n; i++, pa++, pa2++ ) var[i] = *pa2 - SQR( *pa );
}

//
// Draw line
//

void DrawLine( BYTE *img, int w, int h, int stride, 
			  int x1, int y1, int x2, int y2, int r, int g, int b ){
	int i, prev, next, 
		x = abs( x2 - x1 ),
		y = abs( y2 - y1 );
	float m, round;
	BYTE *bit;

	// line from p1 to p2
	if( x > y ){							// go alone x
		// let x2 > x1
		if( x1 > x2 ){
			next = x1;	x1 = x2;	x2 = next;	// swap x
			next = y1;	y1 = y2;	y2 = next;	// swap y
		}
		bit = img + y1 * stride + x1 * 3;	// p1
		m = ( float )( y2 - y1 ) / ( float )( x2 - x1 );		// slope
		round = ( m > 0.0f ) ? 0.5f : -0.5f;
		y = prev = y1;						// previous y value
		for( i = 0; i <= x; i++, bit += 3 ){
			bit[0] = b;
			bit[1] = g;
			bit[2] = r;
			// determine y value
			next = ( int )( m * ( float )i + round ) + y1;
			if( next > prev ) {
				y = next;
				bit += stride;
			}
			else if( next < prev ){
				y = next;
				bit -= stride;
			}
			prev = next;
		}
	}
	else{									// go alone y
		// let y2 > y1
		if( y1 > y2 ){
			next = x1;	x1 = x2;	x2 = next;	// swap x
			next = y1;	y1 = y2;	y2 = next;	// swap y
		}
		bit = img + y1 * stride + x1 * 3;	// p1
		m = ( float )( x2 - x1 ) / ( float )( y2 - y1 );		// 1 / slope
		round = ( m > 0.0f ) ? 0.5f : -0.5f;
		x = prev = x1;						// previous y value
		for( i = 0; i <= y; i++, bit += stride ){
			bit[0] = b;
			bit[1] = g;
			bit[2] = r;
			// determine y value
			next = ( int )( m * ( float )i + round ) + x1;
			if( next > prev ) {
				x = next;
				bit += 3;
			}
			else if( next < prev ){
				x = next;
				bit -= 3;
			}
			prev = next;
		}
	}
	bit = img + y2 * stride + x2 * 3;
	bit[0] = b;
	bit[1] = g; 
	bit[2] = r;
}

void DrawLine( BITMAP *img, int x1, int y1, int x2, int y2, int r, int g, int b ){
	DrawLine( ( BYTE* )img->bmBits, img->bmWidth, img->bmHeight, img->bmWidthBytes,
		x1, y1, x2, y2, r, g, b );
}

void DrawLine( char *img, int w, int h, int x1, int y1, int x2, int y2, char v ){
	int i, prev, next, 
		x = abs( x2 - x1 ),
		y = abs( y2 - y1 );
	float m, r;
	char *bit;

	// line from p1 to p2
	if( x > y ){							// go alone x
		// let x2 > x1
		if( x1 > x2 ){
			next = x1;	x1 = x2;	x2 = next;	// swap x
			next = y1;	y1 = y2;	y2 = next;	// swap y
		}
		bit = img + y1 * w + x1;	// p1
		m = ( float )( y2 - y1 ) / ( float )( x2 - x1 );		// slope
		r = ( m > 0.0f ) ? 0.5f : -0.5f;
		y = prev = y1;						// previous y value
		for( i = 0; i <= x; i++, bit ++ ){
			bit[0] = v;
			// determine y value
			next = ( int )( m * ( float )i + r ) + y1;
			if( next > prev ) {
				y = next;
				bit += w;
			}
			else if( next < prev ){
				y = next;
				bit -= w;
			}
			prev = next;
		}
	}
	else{									// go alone y
		// let y2 > y1
		if( y1 > y2 ){
			next = x1;	x1 = x2;	x2 = next;	// swap x
			next = y1;	y1 = y2;	y2 = next;	// swap y
		}
		bit = img + y1 * w + x1;	// p1
		m = ( float )( x2 - x1 ) / ( float )( y2 - y1 );		// 1 / slope
		r = ( m > 0.0f ) ? 0.5f : -0.5f;
		x = prev = x1;						// previous y value
		for( i = 0; i <= y; i++, bit += w ){
			bit[0] = v;
			// determine y value
			next = ( int )( m * ( float )i + r ) + x1;
			if( next > prev ) {
				x = next;
				bit ++;
			}
			else if( next < prev ){
				x = next;
				bit --;
			}
			prev = next;
		}
	}
	// endpoint setting
	bit = img + y2 * w + x2;
	bit[0] = v;
}

//
// Dither
//

void Dither( BYTE* img, int w, int h, int stride, double* threshold, int m_w, int m_h ){
	int i, j, x, y;
	BYTE *bit;

	for( i = 0; i < h; i++ ){
		y = i % m_h;
		bit = img + i * stride;
		for( j = 0; j < w; j++, bit += 3 ){
			x = j % m_w;
			bit[0] = ( BYTE )( ( ( int )bit[0] + ( int )bit[1] + ( int )bit[2] ) / 3 );
			if( bit[0] > threshold[ y * m_w + x ] * 255 )
				bit[0] = bit[1] = bit[2] = 0xFF;	//white
			else
				bit[0] = bit[1] = bit[2] = 0x00;	//black
		}
	}
}

void Dither( BITMAP *img, double* threshold, int m_w, int m_h ){
	Dither( ( BYTE* )img->bmBits, img->bmHeight, img->bmWidth, img->bmWidthBytes, 
		threshold, m_w, m_h );
}

//
// Contrast
//

void Contrast( BYTE* img, int w, int h, int stride, BYTE threshold, double ratio ){
	int i, j;
	BYTE *bit;
	for( i = 0; i < h; i++ ){
		bit = img + i * stride;
		for( j = 0; j < w; j++, bit++ ){
			if( *bit > threshold )
				*bit += BYTE( ( 0xFF - *bit ) * ratio );	//light
			else
				*bit -= BYTE( *bit * ratio );	//dark
		}
	}
}

void Contrast( BITMAP *img, BYTE threshold, double ratio ){
	Contrast( ( BYTE* )img->bmBits, img->bmWidth, img->bmWidthBytes, img->bmWidthBytes, 
		threshold, ratio );
}

//
// Convolution
//

#define k( i, j ) kernel[ ( (i) + ( size / 2 ) ) * size + (j) + (size / 2 )]

void Convolution( BYTE *img, int w, int h, int stride, int size, double *kernel ){
	BYTE *output = img,
		*input = ( BYTE* )malloc( h * stride * sizeof( BYTE ) );
	double Rsum, Gsum, Bsum, ksum, *k;
	BYTE *inbit, *outbit;
	int i, j, x, y, size2;

	size2 = size / 2;

	// initialize kernal
	k = kernel;
	for( i = 0, ksum = 0; i < size; i++ ){
		for( j = 0; j < size; j++ ) {
			ksum += *k;
			k++;
		}
	}

	// initialize input image
	x = h * stride;
	for( i = 0; i < x; i++ ) input[i] = output[i];

	// convolution
	for( y = 0; y < h; y++ ){
		inbit = input + y * stride;
		outbit = output + y * stride;
		for( x = 0; x < w; x++, inbit += 3, outbit += 3 ){
			if( y < size2 ||
				y >= h - size2 ||
				x < size2 ||
				x >= w - size2 ) {
				// don't change color
				outbit[0] = inbit[0];	outbit[1] = inbit[1];	outbit[2] = inbit[2];
			}
			else{
				for( i = -size2, Rsum = Gsum = Bsum = 0; i <= size2; i++ ){
					for( j = -size2; j <= size2; j++){
						Bsum += k( i, j ) * input[ ( y + i ) * stride + ( x + j ) * 3 + 0 ];
						Gsum += k( i, j ) * input[ ( y + i ) * stride + ( x + j ) * 3 + 1 ];
						Rsum += k( i, j ) * input[ ( y + i ) * stride + ( x + j ) * 3 + 2 ];
					}
				}
				outbit[0] = ( BYTE )( Bsum/ksum );
				outbit[1] = ( BYTE )( Gsum/ksum );
				outbit[2] = ( BYTE )( Rsum/ksum );
			}
		}
	}
//	img = output;
	free( input );
}

void Convolution( BYTE *img, int w, int h, int stride, int size, float *kernel ){
	BYTE *output = img,
		*input = ( BYTE* )malloc( h * stride * sizeof( BYTE ) );
	float Rsum, Gsum, Bsum, ksum, *k;
	BYTE *inbit, *outbit;
	int i, j, x, y, size2;

	size2 = size / 2;

	// initialize kernal
	k = kernel;
	for( i = 0, ksum = 0; i < size; i++ ){
		for( j = 0; j < size; j++ ) {
			ksum += *k;
			k++;
		}
	}

	// initialize input image
	x = h * stride;
	for( i = 0; i < x; i++ ) input[i] = output[i];

	// convolution
	for( y = 0; y < h; y++ ){
		inbit = input + y * stride;
		outbit = output + y * stride;
		for( x = 0; x < w; x++, inbit += 3, outbit += 3 ){
			if( y < size2 ||
				y >= h - size2 ||
				x < size2 ||
				x >= w - size2 ) {
				// don't change color
				outbit[0] = inbit[0];	outbit[1] = inbit[1];	outbit[2] = inbit[2];
			}
			else{
				for( i = -size2, Rsum = Gsum = Bsum = 0; i <= size2; i++ ){
					for( j = -size2; j <= size2; j++){
						Bsum += k( i, j ) * input[ ( y + i ) * stride + ( x + j ) * 3 + 0 ];
						Gsum += k( i, j ) * input[ ( y + i ) * stride + ( x + j ) * 3 + 1 ];
						Rsum += k( i, j ) * input[ ( y + i ) * stride + ( x + j ) * 3 + 2 ];
					}
				}
				outbit[0] = ( BYTE )( Bsum/ksum );
				outbit[1] = ( BYTE )( Gsum/ksum );
				outbit[2] = ( BYTE )( Rsum/ksum );
			}
		}
	}
//	img = output;
	free( input );
}

void Convolution(BITMAP *img, int size, double *kernel){
	Convolution( ( BYTE*)img->bmBits, img->bmWidth, img->bmHeight, 
		img->bmWidthBytes, size, kernel);
}

void Convolution( BITMAP *img, int size, float *kernel ){
	Convolution( ( BYTE*)img->bmBits, img->bmWidth, img->bmHeight, 
		img->bmWidthBytes, size, kernel);
}

#undef k

//
// Gaussian Filter
//

void GaussFilter( BYTE *img, int w, int h, int stride, int win ){
	float *g;
	win = ( win / 2 ) * 2 + 1;
	g = ( float* )malloc( SQR( win ) * sizeof( float ) );
	GenGaussian( g, win );
	Convolution( img, w, h, stride, win, g );
	free( g );
}

void GaussFilter( BITMAP *img, int win ){
	GaussFilter( ( BYTE* )img->bmBits, img->bmWidth, img->bmHeight, img->bmWidthBytes, win );
}

//
// Gradient Image
//

void GradientX( BYTE *img, int w, int h, int stride, int *g ){
	int i, j, *dx,
		n = w * h, w2 = w - 1, h2 = h - 1;
	BYTE *bit;

	// initialize
	for( i = 0; i < n; i++ ) g[i] = 0;

	// calculate gradient
	for( i = 1; i < h2; i++ ){
		bit = img + i * stride + 3;
		dx = g + i * w + 1;
		for( j = 1; j < w2; j++, bit += 3, dx++ ){
			*dx = 2 * ( ( int )bit[ 3 + 0 ] + ( int )bit[ 3 + 1 ] + ( int )bit[ 3 + 2 ] ) + 
				( int )bit[ 3 - stride + 0 ] + ( int )bit[ 3 - stride + 1 ] + ( int )bit[ 3 - stride + 2 ] + 
				( int )bit[ 3 + stride + 0 ] + ( int )bit[ 3 + stride + 1 ] + ( int )bit[ 3 + stride + 2 ];
			*dx -= 2 * ( ( int )bit[ -3 + 0 ] + ( int )bit[ -3 + 1 ] + ( int )bit[ -3 + 2 ] ) + 
				( int )bit[ -3 - stride + 0 ] + ( int )bit[ -3 - stride + 1 ] + ( int )bit[ -3 - stride + 2 ] + 
				( int )bit[ -3 + stride + 0 ] + ( int )bit[ -3 + stride + 1 ] + ( int )bit[ -3 + stride + 2 ];
			*dx /= 24;
		}
	}
}

void GradientY( BYTE *img, int w, int h, int stride, int *g ){
	int i, j, *dy,
		n = w * h, w2 = w - 1, h2 = h - 1;
	BYTE *bit;

	// initialize
	for( i = 0; i < n; i++ ) g[i] = 0;

	// calculate gradient
	for( i = 1; i < h2; i++ ){
		bit = img + i * stride + 3;
		dy = g + i * w + 1;
		for( j = 1; j < w2; j++, bit += 3, dy++ ){
			*dy = 2 * ( ( int )bit[ stride + 0 ] + ( int )bit[ stride + 1 ] + ( int )bit[ stride + 2 ] ) +
				( int )bit[ stride + 3 + 0 ] + ( int )bit[ stride + 3 + 1 ] + ( int )bit[ stride + 3 + 2 ] +
				( int )bit[ stride - 3 + 0 ] + ( int )bit[ stride - 3 + 1 ] + ( int )bit[ stride - 3 + 2 ];
			*dy -= 2 * ( ( int )bit[ -stride + 0 ] + ( int )bit[ -stride + 1 ] + ( int )bit[ -stride + 2 ] ) +
				( int )bit[ -stride + 3 + 0 ] + ( int )bit[ -stride + 3 + 1 ] + ( int )bit[ -stride + 3 + 2 ] +
				( int )bit[ -stride - 3 + 0 ] + ( int )bit[ -stride - 3 + 1 ] + ( int )bit[ -stride - 3 + 2 ];
			*dy /= 24;
		}
	}
}

//
// Bilateral Filtering
//

void BilateralFilterGrey( BYTE *img, int w, int h, int win, int n ){
	BilateralFilterGrey( img, img, w, h, win, n );
}

void BilateralFilterGrey( BYTE *img1, BYTE *img2, int w, int h, int win, int n ){
	int i, j, k, b, t, x, y, p, w2, h2,
		win2 = win / 2,
		height = 256 / n, layb = height / 2;
	float *G, *g, *layer, *layer2, kp, temp, fg, u, v,
		alpha = ( float )50.0f;
	BYTE *bit;
	void *mem;

	// initialize and memory pool
	win = win2 * 2 + 1;
	mem = malloc( ( win * win + n * w * h ) * sizeof( float ) );
	G = ( float* )mem;					// win * win
	g = G + win * win;					// n * w * h

	GenGaussian( G, win );

	// generate intensity kernel
	layer = g;
	for( k = 0; k < n; k++ ){
		b = k * height;					// bottom of layer
		t = b + height;					// top of layer
		p = b + layb;					// mid of layer
		bit = img2;
		j = w * h;
		for( i = 0; i < j; i++, bit++, layer++ ) {
			// || Ip - Iq ||
			temp = ( float )abs( ( p - ( int )( *bit ) ) );	
			*layer = Gauss( temp, alpha );
		}
	}

	// bilateral filtering
	h2 = h - win2;
	w2 = w - win2;
	for( y = win2; y < h2; y++ ){
		bit = img1 + y * w + win2;
		for( x = win2; x < w2; x++, bit++ ){
			p = *bit;					// current pixel intensity
			temp = 0.0f;				// smoothed pixel internsity
			kp = 0.0f;					// normalization factor
			// select intensity layer
			for( i = k = 0; i < n; i++, k += height ){
				if( p >= k && p <= k + height ) break;
			}
			if( i == n ) i--;
			layer = g + i * w * h;		// i-th layer
			// interpolate
			if( i == n - 1 || i == 0 ){	// topest or lowest layer, can't interpolate
				for( i = -win2; i <= win2; i++ ){
					for( j = -win2; j <= win2; j++ ){
						fg = G[ ( i + win2 ) * win + j + win2 ] * layer[ ( y + i ) * w + x + j ];
						kp += fg;
						temp += img1[ ( y + i ) * w + x + j ] * fg;
					}
				}
			}
			else{						// interpolate i and i+1 layer
				k += height / 2;		// mid intensity of layer
				if( p < k ){			// interpolate with lower layer
					// u:i , v:i-1
					layer2 = layer - w * h;		// i-1 -th layer
					v = ( float )( ( k - p ) % height ) / ( float )height;
					u = 1.0f - v;
				}
				else{
					layer2 = layer + w * h;		// i+1 -th layer
					v = ( float )( ( p - k ) % height ) / ( float )height;
					u = 1.0f - v;
				}
				for( i = -win2; i <= win2; i++ ){
					for( j = -win2; j <= win2; j++ ){
						b = ( y + i ) * w + x + j;
						fg = u * layer[b] + v * layer2[b];
						fg *= G[ ( i + win2 ) * win + j + win2 ];
						kp += fg;
						temp += img1[b] * fg;
					}
				}
			}
			*bit = ( BYTE )( temp / kp );
		}
	}

	// free memory pool
	free( mem );
}

void BilateralFilterGrey( float *img, float *img1, float *img2, int w, int h, int win, int n ){
	int i, j, k, x, y, w2, h2, win2 = win / 2;
	float *G, *g, *layer, *layer2, 
		kp, temp, fg, u, v, p, b, t,
		height = 1.0f / ( float )n, layb = height / 2.0f,
		alpha = ( float ).25f;
	float *bit;
	void *mem;

	// initialize and memory pool
	win = win2 * 2 + 1;
	mem = malloc( ( win * win + n * w * h ) * sizeof( float ) );
	G = ( float* )mem;					// win * win
	g = G + win * win;					// n * w * h

	GenGaussian( G, win );

	// generate intensity kernel
	layer = g;
	for( k = 0; k < n; k++ ){
		b = ( float )k * height;		// bottom of layer
		t = b + height;					// top of layer
		p = b + layb;					// mid of layer
		bit = img2;
		j = w * h;
		for( i = 0; i < j; i++, bit++, layer++ ) {
			// || Ip - Iq ||
			temp = ( float )fabs( p - *bit );	
			*layer = Gauss( temp, alpha );
		}
	}

	// bilateral filtering
	h2 = h - win2;
	w2 = w - win2;
	for( y = win2; y < h2; y++ ){
		bit = img1 + y * w + win2;
		for( x = win2; x < w2; x++, bit++ ){
			p = *bit;					// current pixel intensity
			temp = 0.0f;				// smoothed pixel internsity
			kp = 0.0f;					// normalization factor
			// select intensity layer
			for( i = 0, b = 0.0f; i < n; i++, b += height ){
				if( p >= b && p <= b + height ) break;
			}
			if( i == n ) i--;
			layer = g + i * w * h;		// i-th layer
			// interpolate
			if( i == n - 1 || i == 0 ){	// topest or lowest layer, can't interpolate
				for( i = -win2; i <= win2; i++ ){
					for( j = -win2; j <= win2; j++ ){
						fg = G[ ( i + win2 ) * win + j + win2 ] * layer[ ( y + i ) * w + x + j ];
						kp += fg;
						temp += img[ ( y + i ) * w + x + j ] * fg;
					}
				}
			}
			else{						// interpolate i and i+1 layer
				b += height / 2.0;		// mid intensity of layer
				if( p < b ){			// interpolate with lower layer
					// u:i , v:i-1
					layer2 = layer - w * h;		// i-1 -th layer
					v = ( float )( b - p ) / ( float )height;
					u = 1.0f - v;
				}
				else{
					layer2 = layer + w * h;		// i+1 -th layer
					v = ( float )( p - b ) / ( float )height;
					u = 1.0f - v;
				}
				for( i = -win2; i <= win2; i++ ){
					for( j = -win2; j <= win2; j++ ){
						k = ( y + i ) * w + x + j;
						fg = u * layer[k] + v * layer2[k];
						fg *= G[ ( i + win2 ) * win + j + win2 ];
						kp += fg;
						temp += img[k] * fg;
					}
				}
			}
			*bit = temp / kp;
		}
	}

	// free memory pool
	free( mem );
}

void BilateralFilterGrey( double *img, double *img1, double *img2, int w, int h, int win, int n ){
	int i, j, k, x, y, w2, h2, win2 = win / 2;
	double *G, *g, *layer, *layer2, 
		kp, temp, fg, u, v, p, b, t,
		height = 1.0f / ( double )n, layb = height / 2.0,
		alpha = ( double ).25;
	double *bit;
	void *mem;

	// initialize and memory pool
	win = win2 * 2 + 1;
	mem = malloc( ( win * win + n * w * h ) * sizeof( double ) );
	G = ( double* )mem;					// win * win
	g = G + win * win;					// n * w * h

	GenGaussian( G, win );

	// generate intensity kernel
	layer = g;
	for( k = 0; k < n; k++ ){
		b = ( double )k * height;		// bottom of layer
		t = b + height;					// top of layer
		p = b + layb;					// mid of layer
		bit = img2;
		j = w * h;
		for( i = 0; i < j; i++, bit++, layer++ ) {
			// || Ip - Iq ||
			temp = fabs( p - *bit );	
			*layer = Gauss( temp, alpha );
		}
	}

	// bilateral filtering
	h2 = h - win2;
	w2 = w - win2;
	for( y = win2; y < h2; y++ ){
		bit = img1 + y * w + win2;
		for( x = win2; x < w2; x++, bit++ ){
			p = *bit;					// current pixel intensity
			temp = 0.0f;				// smoothed pixel internsity
			kp = 0.0f;					// normalization factor
			// select intensity layer
			for( i = 0, b = 0.0f; i < n; i++, b += height ){
				if( p >= b && p <= b + height ) break;
			}
			if( i == n ) i--;
			layer = g + i * w * h;		// i-th layer
			// interpolate
			if( i == n - 1 || i == 0 ){	// topest or lowest layer, can't interpolate
				for( i = -win2; i <= win2; i++ ){
					for( j = -win2; j <= win2; j++ ){
						fg = G[ ( i + win2 ) * win + j + win2 ] * layer[ ( y + i ) * w + x + j ];
						kp += fg;
						temp += img[ ( y + i ) * w + x + j ] * fg;
					}
				}
			}
			else{						// interpolate i and i+1 layer
				b += height / 2.0;		// mid intensity of layer
				if( p < b ){			// interpolate with lower layer
					// u:i , v:i-1
					layer2 = layer - w * h;		// i-1 -th layer
					v = ( double )( b - p ) / ( double )height;
					u = 1.0f - v;
				}
				else{
					layer2 = layer + w * h;		// i+1 -th layer
					v = ( double )( p - b ) / ( double )height;
					u = 1.0f - v;
				}
				for( i = -win2; i <= win2; i++ ){
					for( j = -win2; j <= win2; j++ ){
						k = ( y + i ) * w + x + j;
						fg = u * layer[k] + v * layer2[k];
						fg *= G[ ( i + win2 ) * win + j + win2 ];
						kp += fg;
						temp += img[k] * fg;
					}
				}
			}
			*bit = temp / kp;
		}
	}

	// free memory pool
	free( mem );
}

void BilateralFilter( BYTE *img, int w, int h, int stride, int win, int segment ){
	int i, j, k;
	BYTE *rgb = ( BYTE* )malloc( h * w * sizeof( BYTE ) ),
		*RGB, *pixel;

	for( k = 0; k < 3; k++ ){		// RGB
		RGB = rgb;
		for( i = 0; i < h; i++ ){
			pixel = img + i * stride;
			for( j = 0; j < w; j++, RGB++, pixel += 3 ){
				*RGB = pixel[k];	// only red
			}
		}
		BilateralFilterGrey( rgb, w, h, win, segment );
		RGB = rgb;
		for( i = 0; i < h; i++ ){
			pixel = img + i * stride;
			for( j = 0; j < w; j++, RGB++, pixel += 3 ){
				pixel[k] = *RGB;	// only red
			}
		}
	}

	free( rgb );
}


void BilateralFilter( BYTE *src, BYTE *tar, int w, int h, int stride, int win, int segment ){
	int i, j, k;
	BYTE *rgb1 = ( BYTE* )malloc( 2 * h * w * sizeof( BYTE ) ),
		*rgb2 = rgb1 + h * w,
		*RGB1, *RGB2, *bit1, *bit2;

	for( k = 0; k < 3; k++ ){		// RGB
		RGB1 = rgb1;
		RGB2 = rgb2;
		for( i = 0; i < h; i++ ){
			bit1 = src + i * stride;
			bit2 = tar + i * stride;
			for( j = 0; j < w; j++, RGB1++, RGB2++, bit1 += 3, bit2 += 3 ){
				*RGB1 = bit1[k];	// only red
				*RGB2 = bit2[k];
			}
		}
		BilateralFilterGrey( rgb1, rgb2, w, h, win, segment );
		RGB1 = rgb1;
		RGB2 = rgb2;
		for( i = 0; i < h; i++ ){
			bit1 = src + i * stride;
			bit2 = tar + i * stride;
			for( j = 0; j < w; j++, RGB1++, RGB2++, bit1 += 3, bit2 += 3 ){
				bit1[k] = *RGB1;	// only red
				bit2[k] = *RGB2;
			}
		}
	}
	free( rgb1 );
}

void BilateralFilter( BITMAP *img, int win, int segment ){
	BilateralFilter( 
		( BYTE* )img->bmBits, img->bmWidth, img->bmHeight, img->bmWidthBytes, win, segment );
}

int BilateralFilter( BITMAP *img1, const BITMAP &img2, int win, int segment ){
	if( img1->bmHeight == img2.bmHeight && img1->bmWidth == img2.bmWidth ){
		BilateralFilter(
			( BYTE* )img1->bmBits, ( BYTE* )img2.bmBits, 
			img1->bmWidth, img1->bmHeight, img1->bmWidthBytes, win, segment );
		return 0;
	}
	else return 1;
}

//
// Warping
//

void Warping( BITMAP *img, double H[9], char p ){
	Warping( ( BYTE* )img->bmBits, img->bmWidth, img->bmHeight, img->bmWidthBytes, H, p );
}

void Warping( BYTE *img, int w, int h, int stride, double H[9], char param ){
	int i, j, k, u, v;
	BYTE *bits, *pixel;
	double x, y, z;
	double Hinv[9], temp[9];

	bits = img;
	pixel = new BYTE [stride * h];
	memcpy( Hinv, H, 9 * sizeof( double ) );
//	inverse( H, Hinv, 3, temp );
	inverse( Hinv, 3, temp );

	for( i = 0; i < h; i++ ){
		for( j = 0; j < w; j++ ){
			x = Hinv[0] * j + Hinv[1] * i + Hinv[2];
			y = Hinv[3] * j + Hinv[4] * i + Hinv[5];
			z = Hinv[6] * j + Hinv[7] * i + Hinv[8];
			x /= z; y /= z;
			u = (int)x; v = (int)y;

			if( x < 0 || y < 0 ) continue;
			if( param == IMG_NEAREST ){
				//	2 3
				//	0 1
				k = 0;
				if( x - (double)u > 0.5 ) k++;
				if( y - (double)v > 0.5 ) k += 2;
				u += k % 2;
				v += k / 2;
				if( u < 0 || v < 0 || u >= w || v >= h ) continue;
				for( k = 0; k < 3; k++ ) 
					pixel[ i * stride + j * 3 + k ] = bits[ v * stride + u * 3 + k ];
			}
			else if( param == IMG_BILINEAR ){
				if( u < 0 || v < 0 || u >= w - 1 || v >= h - 1 ) continue;
				x = x - (double)u;
				y = y - (double)v;
				for( k = 0; k < 3; k++ ){
					pixel[ i * stride + j * 3 + k ] = ( BYTE )(
						( 1.0 - x ) * ( 1.0 - y ) * bits[ v * stride + u * 3 + k ] +
						x * ( 1.0 - y ) * bits[ v * stride + ( u + 1 ) * 3 + k ] +
						( 1.0 - x ) * y * bits[ ( v + 1 ) * stride + u * 3 + k ] +
						x * y * bits[ ( v + 1 ) * stride + ( u + 1 ) * 3 + k ]
						);
				}
			}
			else if( param == IMG_HERMITE ){
				if( u < 0 || v < 0 || u >= w - 1 || v >= h - 1 ) continue;
				x = x - (double)u;
				y = y - (double)v;
				x = 3.0 * x * x - 2.0 * x * x * x;
				y = 3.0 * y * y - 2.0 * y * y * y;
				for( k = 0; k < 3; k++ ){
					pixel[ i * stride + j * 3 + k ] = ( BYTE )(
						( 1.0 - x ) * ( 1.0 - y ) * bits[ v * stride + u * 3 + k ] +
						x * ( 1.0 - y ) * bits[ v * stride + ( u + 1 ) * 3 + k ] +
						( 1.0 - x ) * y * bits[ ( v + 1 ) * stride + u * 3 + k ] +
						x * y * bits[ ( v + 1 ) * stride + ( u + 1 ) * 3 + k ]
						);
				}
			}
		}
	}
	for( i = 0; i < h * stride; i++ ) bits[i] = pixel[i];
	delete pixel;
}

//
// Select Region
//

#define color( x, y, i )	bits[ (y) * stride + (x) * 3 + (i) ]
#define region( x, y )		select[ (y) * w + (x) ]

void Select( BITMAP *bmp, int x, int y ){
	Select( ( BYTE* )bmp->bmBits, bmp->bmWidth, bmp->bmHeight, bmp->bmWidthBytes, x, y );
}

void Select( BYTE *bits, int w, int h, int stride, int x, int y ){
	int n, i, j;
	bool *select;
	BYTE value[3] = { color( x, y, 0 ), color( x, y, 1 ), color( x, y, 2 ) };

	// initialize region image
	select = ( bool* )calloc( h * w, sizeof( bool ) );	// all element is 0
	region( x, y ) = 1;		// initial point

	do{
		n = 0;
		for( i = 0; i < w; i++ ){
			for( j = 0; j < h; j++ ){
				if( !region( i, j ) &&						// not selected yet
					color( i, j, 0 ) == value[0] && 
					color( i, j, 1 ) == value[1] &&
					color( i, j, 2 ) == value[2] &&	(		// same color
					( i > 0 && region( i - 1, j ) ) ||		// left neighbor
					( i < w - 1 && region( i + 1, j ) ) ||	// right neighbor
					( j > 0 && region( i, j - 1 ) ) ||		// down neighbor
					( j < h - 1 && region( i, j + 1 ) )		// top neighbor 
					) ){
					region( i, j ) = 1;
					n++;
				}
			}
		}
		for( i = w - 1; i >= 0; i-- ){
			for( j = h - 1; j >= 0; j-- ){
				if( !region( i, j ) &&						// not selected yet
					color( i, j, 0 ) == value[0] && 
					color( i, j, 1 ) == value[1] &&
					color( i, j, 2 ) == value[2] &&	(		// same color
					( i > 0 && region( i - 1, j ) ) ||		// left neighbor
					( i < w - 1 && region( i + 1, j ) ) ||	// right neighbor
					( j > 0 && region( i, j - 1 ) ) ||		// down neighbor
					( j < h - 1 && region( i, j + 1 ) )		// top neighbor 
					) ){
					region( i, j ) = 1;
					n++;
				}
			}
		}
	} while( n > 0 );	// until no extended region

	// set black outside the region
	j = stride * h;
	for( i = 0; i < j; i++ ) bits[i] = 0x00;
	for( i = 0; i < w; i++ ){
		for( j = 0; j < h; j++ ){
			if( region( i, j ) ){
				color( i, j, 0 ) = //value[0];
				color( i, j, 1 ) = //value[1];
				color( i, j, 2 ) = 0xff;//value[2];
			}
		}
	}

	free( select );
}
#undef color
#undef region

void Select( BYTE *img, char *region, int w, int h, int stride, int p[], int n ){
	int i, j, *p1, *p2;
	char *r;
	BYTE *bit;

	for( i = 0, j = w * h; i < j; i++ ) region[i] = 0;	// initialize
	
	// draw boundary region
	p1 = p;
	p2 = p + 2;
	for( i = 1; i < n; i++ ) {
		DrawLine( region, w, h, p1[0], p1[1], p2[0], p2[1], -1 );
		p1 = p2;
		p2 += 2;
	}
	p2 = p;
	DrawLine( region, w, h, p1[0], p1[1], p2[0], p2[1], -1 );

	// initialize unselect point
	n = 0;
	r = region;
	for( i = 0; i < w; i++, r++ ) if( *r == 0 ) { *r = 1; n++; break; }
	r = region + ( h - 1 ) * w;
	for( i = 0; i < w; i++, r++ ) if( *r == 0 ) { *r = 1; n++; break; }
	r = region;
	for( i = 0; i < h; i++, r += w ) if( *r == 0 ) { *r = 1; n++; break; }
	r = region + w - 1;
	for( i = 0; i < h; i++, r += w ) if( *r == 0 ) { *r = 1; n++; break; }
	if( n == 0 ){		// all selecting case
		for( i = 0, j = h * w; i < j; i++ ){
			if( region[i] == 0 ) region[i] = 1;
		}
		return ;
	}

	// furfill unbounded region
	do{
		n = 0;			// counting
		// forword
		r = region;
		for( i = 0; i < h; i++ ){
			for( j = 0; j < w; j++, r++ ){
				if( r[0] == 0 && (
					( j > 0 && r[-1] == 1 ) ||
					( j < w - 1 && r[1] == 1 ) ||
					( i > 0 && r[-w] == 1 ) ||
					( i < h - 1 && r[w] == 1 ) ) ){
					r[0] = 1;
					n++;
				}
			}
		}
		// backword
		r--;
		for( i = 0; i < h; i++ ){
			for( j = 0; j < w; j++, r-- ){
				if( r[0] == 0 && (
					( j > 0 && r[-1] == 1 ) ||
					( j < w - 1 && r[1] == 1 ) ||
					( i > 0 && r[-w] == 1 ) ||
					( i < h - 1 && r[w] == 1 ) ) ){
					r[0] = 1;
					n++;
				}
			}
		}
	} while( n );	// until no region added
	
	// black unselected region
	r = region;
	for( i = 0; i < h; i++ ){
		bit = img + i * stride;
		for( j = 0; j < w; j++, r++, bit += 3 ){
			if( r[0] > 0 ) r[0] = 0;		// unselected region
			else if( r[0] == 0 ) r[0] = 1;	// selected region
			// edge = -1
		}
	}
}

void Select( BITMAP *img, char *region, int p[], int n ){
	Select( ( BYTE* )img->bmBits, region, img->bmWidth, img->bmHeight, img->bmWidthBytes, p, n );
}

void Select( BYTE *img, int w, int h, int stride, int p[], int n ){
	int i, j;
	BYTE *bit;
	char *r, *region;

	r = region = ( char* )malloc( w * h * sizeof( char ) );
	Select( img, region, w, h, stride, p, n );
	for( i = 0; i < h; i++ ){
		bit = img + i * stride;
		for( j = 0; j < w; j++, bit += 3, r++ ){
			if( r[0] == 0 ) bit[0] = bit[1] = bit[2] = 0x00;
		}
	}
	free( region );
}

void Select( BITMAP *img, int p[], int n ){
	Select( ( BYTE* )img->bmBits, img->bmWidth, img->bmHeight, img->bmWidthBytes, p, n );
}

//
// Scaling
//

void Scale( BYTE *img1, int w1, int h1, int stride1, BYTE *img2, int w2, int h2, int stride2 ){
	int y1, y2, x1, x2;
	float v1, v2, u1, u2;
	BYTE *bit1, *bit2;

	// inverse warping
	for( y2 = 0; y2 < h2; y2++ ){
		bit2 = img2 + y2 * stride2;
		v2 = ( float )y2 / ( float )h2 * ( float )h1;	
		y1 = (int)v2;					// the old pixel
		v2 = v2 - ( float )y1;			// the ratio of interpolation
		v1 = 1.0 - v2;
		for( x2 = 0; x2 < w2; x2++, bit2 += 3 ){
			u2 = (float)x2 / ( float )w2 * ( float )w1;
			x1 = (int)u2;				// the old pixel
			u2 = u2 - ( float )x1;		// the ratio of interpolation
			u1 = 1.0 - u2;
			bit1 = img1 + y1 * stride1 + x1 * 3;
			// bilinear
			bit2[0] = u1 * v1 * bit1[0];
			bit2[1] = u1 * v1 * bit1[1];
			bit2[2] = u1 * v1 * bit1[2];
			if( y1 + 1 < h1 ){
				bit2[0] += u1 * v2 * bit1[ stride1 + 0 ];
				bit2[1] += u1 * v2 * bit1[ stride1 + 1 ];
				bit2[2] += u1 * v2 * bit1[ stride1 + 2 ];
			}
			if( x1 + 1 < w1 ){
				bit2[0] += u2 * v1 * bit1[ 3 + 0 ];
				bit2[1] += u2 * v1 * bit1[ 3 + 1 ];
				bit2[2] += u2 * v1 * bit1[ 3 + 2 ];
			}
			if( x1 + 1 < w1 && y1 + 1 < h1 ){
				bit2[0] += u2 * v2 * bit1[ 3 + stride1 + 0 ] ;
				bit2[1] += u2 * v2 * bit1[ 3 + stride1 + 1 ] ;
				bit2[2] += u2 * v2 * bit1[ 3 + stride1 + 2 ] ;
			}
		}
	}
}

void Scale8( float *img1, int w1, int h1, float *img2, int w2, int h2 ){
	int y1, y2, x1, x2;
	float v1, v2, u1, u2, *bit1, *bit2;

	// inverse warping
	for( y2 = 0; y2 < h2; y2++ ){
		bit2 = img2 + y2 * w2;
		v2 = ( float )y2 / ( float )h2 * ( float )h1;	
		y1 = (int)v2;					// the old pixel
		v2 = v2 - ( float )y1;			// the ratio of interpolation
		v1 = 1.0 - v2;
		for( x2 = 0; x2 < w2; x2++, bit2++ ){
			u2 = (float)x2 / ( float )w2 * ( float )w1;
			x1 = (int)u2;				// the old pixel
			u2 = u2 - ( float )x1;		// the ratio of interpolation
			u1 = 1.0 - u2;
			bit1 = img1 + y1 * w1 + x1;
			// bilinear
			bit2[0] = u1 * v1 * bit1[0];
			if( y1 + 1 < h1 )
				bit2[0] += u1 * v2 * bit1[ w1 ];
			if( x1 + 1 < w1 )
				bit2[0] += u2 * v1 * bit1[ 1 ];
			if( x1 + 1 < w1 && y1 + 1 < h1 )
				bit2[0] += u2 * v2 * bit1[ 1 + w1 ] ;
		}
	}
}

void Scale8( double *img1, int w1, int h1, double *img2, int w2, int h2 ){
	int y1, y2, x1, x2;
	double v1, v2, u1, u2, *bit1, *bit2;

	// inverse warping
	for( y2 = 0; y2 < h2; y2++ ){
		bit2 = img2 + y2 * w2;
		v2 = ( double )y2 / ( double )h2 * ( double )h1;	
		y1 = (int)v2;					// the old pixel
		v2 = v2 - ( double )y1;			// the ratio of interpolation
		v1 = 1.0 - v2;
		for( x2 = 0; x2 < w2; x2++, bit2++ ){
			u2 = (double)x2 / ( double )w2 * ( double )w1;
			x1 = (int)u2;				// the old pixel
			u2 = u2 - ( double )x1;		// the ratio of interpolation
			u1 = 1.0 - u2;
			bit1 = img1 + y1 * w1 + x1;
			// bilinear
			bit2[0] = u1 * v1 * bit1[0];
			if( y1 + 1 < h1 )
				bit2[0] += u1 * v2 * bit1[ w1 ];
			if( x1 + 1 < w1 )
				bit2[0] += u2 * v1 * bit1[ 1 ];
			if( x1 + 1 < w1 && y1 + 1 < h1 )
				bit2[0] += u2 * v2 * bit1[ 1 + w1 ] ;
		}
	}
}

/**/
