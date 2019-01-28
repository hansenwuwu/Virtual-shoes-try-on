#include "imgproc.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

// 
// Macro
//

#define M_PI		3.1415926535897932384626433832795
#define SQR( x )	( (x) * (x) )
#define MIN( a, b )		( ( (a) < (b) ? (a) : (b) )
#define MAX( a, b )		( ( (a) > (b) ? (a) : (b) )
#define MIN3( a, b, c ) ( ( (a) < (b) ) ? MIN( (a), (c) ) : MIN( (b), (c) ) )
#define MAX3( a, b, c ) ( ( (a) > (b) ) ? MAX( (a), (c) ) : MAX( (b), (c) ) )

//
// inline function
//

inline int InImage( const ImageBuffer &buf, const int x, const int y ) {
	if( x < 0 || y < 0 || x >= buf.w || y >= buf.h ) return 0;
	else return 1;
}

inline int InImage( const Selected &buf, const int x, const int y ) {
	if( x < 0 || y < 0 || x >= buf.w || y >= buf.h ) return 0;
	else return 1;
}

inline int SelectImage( const Selected &buf, const int x, const int y ){
	if( InImage( buf, x, y ) == 0 ) return 0;
	buf.img[ y * buf.w + x ] = ( char )0xff;
	return 1;
}

inline void GenGaussian( float *g, int w ){
	int a = w / 2, i, j;
	const float 
		b = 0.5f / ( float )SQR( ( float )a * 0.33f ), 
		s = b / ( float )M_PI;
	float x;

	if( w == 3 ){
		g[0] = g[2] = g[6] = g[8] = 0.05882352941176470588235294117647f;
		g[1] = g[3] = g[5] = g[7] = 0.11764705882352941176470588235294f;
		g[4] = 0.29411764705882352941176470588235f;
		return;
	}

	for( i = -a; i <= a; i++ ){
		for( j = i; j <= a; j++ ){
			x = ( float )( SQR( i ) + SQR( j ) );
			g[ ( a + i ) * w + a + j ] = 
			g[ ( a + j ) * w + a + i ] = ( float )( exp( -x * b ) * s );
		}
	}
}

//               //
// Type Template //
//               //

//
// Select Region
//

template < typename T >
int Select1( ImageBuffer buf, Selected sel, ImagePoint p, float thres, T* bits ){
	ImagePoint *stack;
	T *bit;
	char *sel_bit;
	int i = 0;

	if( buf.img == 0x00 || sel.img == 0x00 ) return 0;	// check image
	memset( sel.img, 0x00, sel.w * sel.h );				// select nothing
	sel.img[ p.y * sel.w + p.x ] = ( char )0xff;		// starting point
	if( InImage( buf, p.x, p.y ) == 0 ) return 0;		// if starting point is feasible
	stack = ( ImagePoint* )malloc( sel.w * sel.h * sizeof( ImagePoint ) );

#define push_stack( a, b ) i++; stack[i].x = (a); stack[i].y = (b);
	stack[i] = p;				// push stack
	while( i > -1 ){			// if stack is not empty
		p = stack[i];			// top of stack
		i--;					// pop stack
		bit = bits + p.y * buf.stride + p.x;
		sel_bit = sel.img + p.y * sel.w + p.x;
		if( InImage( buf, p.x + 1, p.y ) &&						// right
			sel_bit[ 1 ] == 0x00 &&			// non-selected
			fabs( ( float )bit[1] - ( float )bit[0] ) < thres ){	// is neighbor
			SelectImage( sel, p.x + 1, p.y );
			push_stack( p.x + 1, p.y );
		}
		if( InImage( buf, p.x - 1, p.y ) &&						// left
			sel_bit[ -1 ] == 0x00 &&			// non-selected
			fabs( ( float )bit[-1] - ( float )bit[0] ) < thres ){	// is neighbor
			SelectImage( sel, p.x - 1, p.y );
			push_stack( p.x - 1, p.y );
		}
		if( InImage( buf, p.x, p.y + 1 ) &&						// top
			sel_bit[ sel.w ] == 0x00 &&		// non-selected
			fabs( ( float )bit[ buf.stride ] - ( float )bit[0] ) < thres ){	// is neighbor
			SelectImage( sel, p.x, p.y + 1 );
			push_stack( p.x, p.y + 1 );
		}
		if( InImage( buf, p.x, p.y - 1 ) &&						// down
			sel_bit[ -sel.w ] == 0x00 &&			// non-selected
			fabs( ( float )bit[ -buf.stride ] - ( float )bit[0] ) < thres ){	// is neighbor
			SelectImage( sel, p.x, p.y - 1 );
			push_stack( p.x, p.y - 1 );
		}
	}
#undef push_stack
	free( stack );
	return 1;
}

template < typename T >
int Select3( ImageBuffer buf, Selected sel, ImagePoint p, float thres, T* bits ){
	ImagePoint *stack;
	T *bit;
	char *sel_bit;
	int i = 0;

	if( buf.img == 0x00 || sel.img == 0x00 ) return 0;	// check image
	memset( sel.img, 0x00, sel.w * sel.h );				// select nothing
	sel.img[ p.y * sel.w + p.x ] = ( char )0xff;		// starting point
	if( InImage( buf, p.x, p.y ) == 0 ) return 0;		// if starting point is feasible
	stack = ( ImagePoint* )malloc( sel.w * sel.h * sizeof( ImagePoint ) );

#define push_stack( a, b ) i++; stack[i].x = (a); stack[i].y = (b);
	stack[i] = p;				// push stack
	while( i > -1 ){			// if stack is not empty
		p = stack[i];			// top of stack
		i--;					// pop stack
		bit = bits + p.y * buf.stride + p.x * 3;
		sel_bit = sel.img + p.y * sel.w + p.x;
		if( InImage( buf, p.x + 1, p.y ) &&						// right
			sel_bit[ 1 ] == 0x00 &&			// non-selected
			.33f * fabs( ( float )bit[3] + ( float )bit[4] + ( float )bit[5]
			- ( float )bit[0] - ( float )bit[1] - ( float )bit[2] ) < thres ){	// is neighbor
			SelectImage( sel, p.x + 1, p.y );
			push_stack( p.x + 1, p.y );
		}
		if( InImage( buf, p.x - 1, p.y ) &&						// left
			sel_bit[ -1 ] == 0x00 &&			// non-selected
			.33f * fabs( ( float )bit[-3] + ( float )bit[-2] + ( float )bit[-1] 
			- ( float )bit[0] - ( float )bit[1] - ( float )bit[2] ) < thres ){	// is neighbor
			SelectImage( sel, p.x - 1, p.y );
			push_stack( p.x - 1, p.y );
		}
		if( InImage( buf, p.x, p.y + 1 ) &&						// top
			sel_bit[ sel.w ] == 0x00 &&		// non-selected
			.33f * fabs( ( float )bit[ buf.stride ] + ( float )bit[ buf.stride + 1 ] + ( float )bit[ buf.stride + 2 ]
			- ( float )bit[0] - ( float )bit[1] - ( float )bit[2] ) < thres ){	// is neighbor
			SelectImage( sel, p.x, p.y + 1 );
			push_stack( p.x, p.y + 1 );
		}
		if( InImage( buf, p.x, p.y - 1 ) &&						// down
			sel_bit[ -sel.w ] == 0x00 &&			// non-selected
			.33f * fabs( ( float )bit[ -buf.stride ] + ( float )bit[ -buf.stride + 1 ] + ( float )bit[ -buf.stride + 2 ]
			- ( float )bit[0] - ( float )bit[1] - ( float )bit[2] ) < thres ){	// is neighbor
			SelectImage( sel, p.x, p.y - 1 );
			push_stack( p.x, p.y - 1 );
		}
	}
#undef push_stack
	free( stack );
	return 1;
}

template < typename T >
int Select4( ImageBuffer buf, Selected sel, ImagePoint p, float thres, T* bits ){
	ImagePoint *stack;
	T *bit;
	char *sel_bit;
	int i = 0;

	if( buf.img == 0x00 || sel.img == 0x00 ) return 0;	// check image
	memset( sel.img, 0x00, sel.w * sel.h );				// select nothing
	sel.img[ p.y * sel.w + p.x ] = ( char )0xff;		// starting point
	if( InImage( buf, p.x, p.y ) == 0 ) return 0;		// if starting point is feasible
	stack = ( ImagePoint* )malloc( sel.w * sel.h * sizeof( ImagePoint ) );

#define push_stack( a, b ) i++; stack[i].x = (a); stack[i].y = (b);
	stack[i] = p;				// push stack
	while( i > -1 ){			// if stack is not empty
		p = stack[i];			// top of stack
		i--;					// pop stack
		bit = bits + p.y * buf.stride + p.x * 4;
		sel_bit = sel.img + p.y * sel.w + p.x;
		if( InImage( buf, p.x + 1, p.y ) &&						// right
			sel_bit[ 1 ] == 0x00 &&			// non-selected
			.33f * fabs( ( float )bit[4] + ( float )bit[5] + ( float )bit[6]
			- ( float )bit[0] - ( float )bit[1] - ( float )bit[2] ) < thres ){	// is neighbor
			SelectImage( sel, p.x + 1, p.y );
			push_stack( p.x + 1, p.y );
		}
		if( InImage( buf, p.x - 1, p.y ) &&						// left
			sel_bit[ -1 ] == 0x00 &&			// non-selected
			.33f * fabs( ( float )bit[-4] + ( float )bit[-3] + ( float )bit[-2] 
			- ( float )bit[0] - ( float )bit[1] - ( float )bit[2] ) < thres ){	// is neighbor
			SelectImage( sel, p.x - 1, p.y );
			push_stack( p.x - 1, p.y );
		}
		if( InImage( buf, p.x, p.y + 1 ) &&						// top
			sel_bit[ sel.w ] == 0x00 &&		// non-selected
			.33f * fabs( ( float )bit[ buf.stride ] + ( float )bit[ buf.stride + 1 ] + ( float )bit[ buf.stride + 2 ]
			- ( float )bit[0] - ( float )bit[1] - ( float )bit[2] ) < thres ){	// is neighbor
			SelectImage( sel, p.x, p.y + 1 );
			push_stack( p.x, p.y + 1 );
		}
		if( InImage( buf, p.x, p.y - 1 ) &&						// down
			sel_bit[ -sel.w ] == 0x00 &&			// non-selected
			.33f * fabs( ( float )bit[ -buf.stride ] + ( float )bit[ -buf.stride + 1 ] + ( float )bit[ -buf.stride + 2 ]
			- ( float )bit[0] - ( float )bit[1] - ( float )bit[2] ) < thres ){	// is neighbor
			SelectImage( sel, p.x, p.y - 1 );
			push_stack( p.x, p.y - 1 );
		}
	}
#undef push_stack
	free( stack );
	return 1;
}

//
// Convolution
//

template < typename T >
int DrawLine1( ImageBuffer buf, int x1, int y1, int x2, int y2, T gray ){
	int i, prev, next,
		x = abs( x2 - x1 ), 
		y = abs( y2 - y1 );
	float m, round;
	T *bit;

	if( InImage( buf, x1, y1 ) == 0 ) return 0;
	if( InImage( buf, x2, y2 ) == 0 ) return 0;

	// line from p1 to p2
	if( x > y ){							// go alone x
		// let x2 > x1
		if( x1 > x2 ){
			next = x1;	x1 = x2;	x2 = next;	// swap x
			next = y1;	y1 = y2;	y2 = next;	// swap y
		}
		bit = ( T* )buf.img + y1 * buf.stride + x1;	// p1
		m = ( float )( y2 - y1 ) / ( float )( x2 - x1 );		// slope
		round = ( m > 0.0f ) ? 0.5f : -0.5f;
		y = prev = y1;						// previous y value
		for( i = 0; i <= x; i++, bit++ ){
			bit[0] = gray;
			// determine y value
			next = ( int )( m * ( float )i + round ) + y1;
			if( next > prev ) {
				y = next;
				bit += buf.stride;
			}
			else if( next < prev ){
				y = next;
				bit -= buf.stride;
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
		bit = ( T* )buf.img + y1 * buf.stride + x1;	// p1
		m = ( float )( x2 - x1 ) / ( float )( y2 - y1 );		// 1 / slope
		round = ( m > 0.0f ) ? 0.5f : -0.5f;
		x = prev = x1;						// previous y value
		for( i = 0; i <= y; i++, bit += buf.stride ){
			bit[0] = gray;
			// determine y value
			next = ( int )( m * ( float )i + round ) + x1;
			if( next > prev ) {
				x = next;
				bit++;
			}
			else if( next < prev ){
				x = next;
				bit++;
			}
			prev = next;
		}
	}
	bit = ( T* )buf.img + y2 * buf.stride + x2;
	bit[0] = gray;
	return 1;
}

template < typename T >
int DrawLine3( ImageBuffer buf, int x1, int y1, int x2, int y2, T red, T green, T blue ){
	int i, prev, next,
		x = abs( x2 - x1 ), 
		y = abs( y2 - y1 );
	float m, round;
	T *bit;

	if( InImage( buf, x1, y1 ) == 0 ) return 0;
	if( InImage( buf, x2, y2 ) == 0 ) return 0;

	// line from p1 to p2
	if( x > y ){							// go alone x
		// let x2 > x1
		if( x1 > x2 ){
			next = x1;	x1 = x2;	x2 = next;	// swap x
			next = y1;	y1 = y2;	y2 = next;	// swap y
		}
		bit = ( T* )buf.img + y1 * buf.stride + x1 * 3;	// p1
		m = ( float )( y2 - y1 ) / ( float )( x2 - x1 );		// slope
		round = ( m > 0.0f ) ? 0.5f : -0.5f;
		y = prev = y1;						// previous y value
		for( i = 0; i <= x; i++, bit += 3 ){
			bit[0] = red;
			bit[1] = green;
			bit[2] = blue;
			// determine y value
			next = ( int )( m * ( float )i + round ) + y1;
			if( next > prev ) {
				y = next;
				bit += buf.stride;
			}
			else if( next < prev ){
				y = next;
				bit -= buf.stride;
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
		bit = ( T* )buf.img + y1 * buf.stride + x1 * 3;	// p1
		m = ( float )( x2 - x1 ) / ( float )( y2 - y1 );		// 1 / slope
		round = ( m > 0.0f ) ? 0.5f : -0.5f;
		x = prev = x1;						// previous y value
		for( i = 0; i <= y; i++, bit += buf.stride ){
			bit[0] = red;
			bit[1] = green;
			bit[2] = blue;
			// determine y value
			next = ( int )( m * ( float )i + round ) + x1;
			if( next > prev ) {
				x = next;
				bit += 3;
			}
			else if( next < prev ){
				x = next;
				bit += 3;
			}
			prev = next;
		}
	}
	bit = ( T* )buf.img + y2 * buf.stride + x2 * 3;
	bit[0] = red;
	bit[1] = green; 
	bit[2] = blue;
	return 1;
}

template < typename T >
int DrawLine4( ImageBuffer buf, int x1, int y1, int x2, int y2, T red, T green, T blue, T alpha ){
	int i, prev, next,
		x = abs( x2 - x1 ), 
		y = abs( y2 - y1 );
	float m, round;
	T *bit;

	if( InImage( buf, x1, y1 ) == 0 ) return 0;
	if( InImage( buf, x2, y2 ) == 0 ) return 0;

	// line from p1 to p2
	if( x > y ){							// go alone x
		// let x2 > x1
		if( x1 > x2 ){
			next = x1;	x1 = x2;	x2 = next;	// swap x
			next = y1;	y1 = y2;	y2 = next;	// swap y
		}
		bit = ( T* )buf.img + y1 * buf.stride + x1 * 4;	// p1
		m = ( float )( y2 - y1 ) / ( float )( x2 - x1 );		// slope
		round = ( m > 0.0f ) ? 0.5f : -0.5f;
		y = prev = y1;						// previous y value
		for( i = 0; i <= x; i++, bit += 4 ){
			bit[0] = red;
			bit[1] = green;
			bit[2] = blue;
			bit[3] = alpha;
			// determine y value
			next = ( int )( m * ( float )i + round ) + y1;
			if( next > prev ) {
				y = next;
				bit += buf.stride;
			}
			else if( next < prev ){
				y = next;
				bit -= buf.stride;
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
		bit = ( T* )buf.img + y1 * buf.stride + x1 * 4;	// p1
		m = ( float )( x2 - x1 ) / ( float )( y2 - y1 );		// 1 / slope
		round = ( m > 0.0f ) ? 0.5f : -0.5f;
		x = prev = x1;						// previous y value
		for( i = 0; i <= y; i++, bit += buf.stride ){
			bit[0] = red;
			bit[1] = green;
			bit[2] = blue;
			bit[3] = alpha;
			// determine y value
			next = ( int )( m * ( float )i + round ) + x1;
			if( next > prev ) {
				x = next;
				bit += 4;
			}
			else if( next < prev ){
				x = next;
				bit += 4;
			}
			prev = next;
		}
	}
	bit = ( T* )buf.img + y2 * buf.stride + x2 * 4;
	bit[0] = red;
	bit[1] = green; 
	bit[2] = blue;
	bit[3] = alpha;
	return 1;
}

template < typename T >
int DrawTriangleLtoR1( ImageBuffer buf, ImagePoint p1, ImagePoint p2, ImagePoint p3, T x1, T x2, T x3 ){
	T *bits = buf.img, *bit1, *bit2;
	int i, j;
	return 1;
}

//
// Filling
//

template < typename T, typename C >
int Filling( ImageBuffer buf, Region r, T *bits, C color ){
	T *bit;
	char *rbit;
	int i, j;

	for( i = 0; i < buf.h; i++ ){
		bit = bits + i * buf.stride;
		rbit = r.img + i * r.w;
		for( j = 0; j < buf.w; j++, bit += buf.color, rbit++ ){
			if( rbit[0] ) memcpy( bit, &color, sizeof( color ) );
		}
	}
	return 1;
}

//
// Convolution
//

template < typename T >
int Convolution( ImageBuffer buf, int size, float *f, T *bits ){
	T *bit, *output, *outbit;
	float temp[4], *pf;
	const int w = size / 2,
		left = w, right = buf.w - w, top = buf.h - w, bottom = w,
		offset = w * ( buf.stride + buf.color ),
		bitpixel = buf.color * sizeof( float );
	int i, j, k, x, y; 

	if( size % 2 == 0 ) return 0;

	// copy image
	output = ( T* )malloc( buf.h * buf.stride * sizeof( T ) );
	if( output == NULL ) return 0;				// fail to allocate
	memcpy( output, buf.img, buf.h * buf.stride * sizeof( T ) );

	for( i = bottom; i < top; i++ ){
		// target pixel
		outbit = output + i * buf.stride + left * buf.color;
		bits = ( T* )buf.img + i * buf.stride + left * buf.color - offset;
		for( j = left; j < right; j++, outbit += buf.color, bits += buf.color ){
			// set to zero
			memset( temp, 0x00, bitpixel );
			// convolution
			pf = f;
			for( x = 0; x < size; x++ ){
				bit = bits + x * buf.stride;
				for( y = 0; y < size; y++, pf++, bit += buf.color ) {
					for( k = 0; k < buf.color; k++ ) temp[k] += pf[0] * ( float )bit[k];
				}
			}
			// set pixel
			for( k = 0; k < buf.color; k++ ) outbit[k] = ( T )temp[k];
		}
	}

	memcpy( buf.img, output, buf.h * buf.stride * sizeof( T ) );
	free( output );
	return 1;
}

//
// Gaussian Filtering
//

template < typename T >
int GaussFilter( ImageBuffer buf, int size, T* bits ){
	float *g;
	int temp;
	
	// check window size
	if( size % 2 == 0 ) return 0;

	// memory allocate
	g = ( float * )malloc( SQR( size ) * sizeof( float ) );
	if( g == NULL ) return 0;
	GenGaussian( g, size );
	
	temp = Convolution( buf, size, g, bits );
	free( g );
	return temp;
}

//
// Scaling
//

template < typename T >
int Scale( ImageBuffer buf, T *img2, int w2, int h2, int stride2 ){
	const int w1 = buf.w, h1 = buf.h,
		stride1 = buf.stride;
	int y1, y2, x1, x2, i;
	float v1, v2, u1, u2, temp[4];
	T *img1 = ( T* )buf.img, *bit1, *bit2;

	// inverse warping
	for( y2 = 0; y2 < h2; y2++ ){
		bit2 = img2 + y2 * stride2;
		v2 = ( float )y2 / ( float )h2 * ( float )h1;
		y1 = ( int )v2;					// the old pixel
		v2 = v2 - ( float )y1;			// the ratio of interpolation
		v1 = 1.0f - v2;
		for( x2 = 0; x2 < w2; x2++, bit2 += buf.color ){
			u2 = ( float )x2 / ( float )w2 * ( float )w1;
			x1 = ( int )u2;				// the old pixel
			u2 = u2 - ( float )x1;		// the ratio of interpolation
			u1 = 1.0f - u2;
			bit1 = img1 + y1 * stride1 + x1 * buf.color;
			// bilinear
			for( i = 0; i < buf.color; i++ ) 
				temp[i] = u1 * v1 * ( float )bit1[i];
			if( y1 + 1 < h1 ){
				for( i = 0; i < buf.color; i++ ) 
					temp[i] += u1 * v2 * ( float )bit1[ stride1 + i ];
			}
			if( x1 + 1 < w1 ){
				for( i = 0; i < buf.color; i++ ) 
					temp[i] += u2 * v1 * ( float )bit1[ buf.color + i ];
			}
			if( x1 + 1 < w1 && y1 + 1 < h1 ){
				for( i = 0; i < buf.color; i++ ) 
					temp[i] += u2 * v2 * ( float )bit1[ buf.color + stride1 + i ] ;
			}
			for( i = 0; i < buf.color; i++ ) bit2[i] = ( T )temp[i];
		}
	}
	return 1;
}

//
// Dithering
//

template < typename T >
void Dither( ImageBuffer buf, T *threshold, T max, T min ){
	int i, j, k;
	T *bit, *bits = ( T* )buf.img;

	for( i = 0; i < buf.h; i++ ){
		bit = bits + i * buf.stride;
		for( j = 0; j < buf.w; j++, bit += buf.color ){
			for( k = 0; k < buf.color; k++ ) bit[k] = ( bit[k] < threshold[k] ) ? min : max;
		}
	}
}

//
// Export Functions
//

//
// DrawLine
//

int DrawLine1b( ImageBuffer buf, ImagePoint p1, ImagePoint p2, Color1b c ){
	switch( buf.color ){
	case 1:	return DrawLine1( buf, p1.x, p1.y, p2.x, p2.y, c.x ); 
	case 3: return DrawLine3( buf, p1.x, p1.y, p2.x, p2.y, c.x, c.x, c.x ); 
	case 4: return DrawLine4( buf, p1.x, p1.y, p2.x, p2.y, c.x, c.x, c.x, ( BYTE )0xff ); 
	default: return 0;
	}
}
int DrawLine3b( ImageBuffer buf, ImagePoint p1, ImagePoint p2, Color3b c ){
	switch( buf.color ){
	case 1:	return 0;
	case 3: return DrawLine3( buf, p1.x, p1.y, p2.x, p2.y, c.r, c.g, c.b );
	case 4: return DrawLine4( buf, p1.x, p1.y, p2.x, p2.y, c.r, c.g, c.b, ( BYTE )0xff );
	default: return 0;
	}
}
int DrawLine4b( ImageBuffer buf, ImagePoint p1, ImagePoint p2, Color4b c ){
	switch( buf.color ){
	case 1:	
	case 3: return 0;
	case 4: return DrawLine4( buf, p1.x, p1.y, p2.x, p2.y, c.r, c.b, c.g, c.a );
	default: return 0;
	}
}
int DrawLine1i( ImageBuffer buf, ImagePoint p1, ImagePoint p2, Color1i c ){
	switch( buf.color ){
	case 1:	return DrawLine1( buf, p1.x, p1.y, p2.x, p2.y, c.x ); 
	case 3: return DrawLine3( buf, p1.x, p1.y, p2.x, p2.y, c.x, c.x, c.x ); 
	case 4: return DrawLine4( buf, p1.x, p1.y, p2.x, p2.y, c.x, c.x, c.x, 255 ); 
	default: return 0;
	}
}
int DrawLine3i( ImageBuffer buf, ImagePoint p1, ImagePoint p2, Color3i c ){
	switch( buf.color ){
	case 1:	return 0;
	case 3: return DrawLine3( buf, p1.x, p1.y, p2.x, p2.y, c.r, c.g, c.b );
	case 4: return DrawLine4( buf, p1.x, p1.y, p2.x, p2.y, c.r, c.g, c.b, 255 );
	default: return 0;
	}
}
int DrawLine4i( ImageBuffer buf, ImagePoint p1, ImagePoint p2, Color4i c ){
	switch( buf.color ){
	case 1:	
	case 3: return 0;
	case 4: return DrawLine4( buf, p1.x, p1.y, p2.x, p2.y, c.r, c.b, c.g, c.a );
	default: return 0;
	}
}
int DrawLine1f( ImageBuffer buf, ImagePoint p1, ImagePoint p2, Color1f c ){
	switch( buf.color ){
	case 1:	return DrawLine1( buf, p1.x, p1.y, p2.x, p2.y, c.x ); 
	case 3: return DrawLine3( buf, p1.x, p1.y, p2.x, p2.y, c.x, c.x, c.x ); 
	case 4: return DrawLine4( buf, p1.x, p1.y, p2.x, p2.y, c.x, c.x, c.x, 255.0f ); 
	default: return 0;
	}
}
int DrawLine3f( ImageBuffer buf, ImagePoint p1, ImagePoint p2, Color3f c ){
	switch( buf.color ){
	case 1:	return 0;
	case 3: return DrawLine3( buf, p1.x, p1.y, p2.x, p2.y, c.r, c.g, c.b );
	case 4: return DrawLine4( buf, p1.x, p1.y, p2.x, p2.y, c.r, c.g, c.b, 255.0f );
	default: return 0;
	}
}
int DrawLine4f( ImageBuffer buf, ImagePoint p1, ImagePoint p2, Color4f c ){
	switch( buf.color ){
	case 1:	
	case 3: return 0;
	case 4: return DrawLine4( buf, p1.x, p1.y, p2.x, p2.y, c.r, c.b, c.g, c.a );
	default: return 0;
	}
}
int DrawLine1d( ImageBuffer buf, ImagePoint p1, ImagePoint p2, Color1d c ){
	switch( buf.color ){
	case 1:	return DrawLine1( buf, p1.x, p1.y, p2.x, p2.y, c.x ); 
	case 3: return DrawLine3( buf, p1.x, p1.y, p2.x, p2.y, c.x, c.x, c.x ); 
	case 4: return DrawLine4( buf, p1.x, p1.y, p2.x, p2.y, c.x, c.x, c.x, 255.0 ); 
	default: return 0;
	}
}
int DrawLine3d( ImageBuffer buf, ImagePoint p1, ImagePoint p2, Color3d c ){
	switch( buf.color ){
	case 1:	return 0;
	case 3: return DrawLine3( buf, p1.x, p1.y, p2.x, p2.y, c.r, c.g, c.b );
	case 4: return DrawLine4( buf, p1.x, p1.y, p2.x, p2.y, c.r, c.g, c.b, 255.0 );
	default: return 0;
	}
}
int DrawLine4d( ImageBuffer buf, ImagePoint p1, ImagePoint p2, Color4d c ){
	switch( buf.color ){
	case 1:	
	case 3: return 0;
	case 4: return DrawLine4( buf, p1.x, p1.y, p2.x, p2.y, c.r, c.b, c.g, c.a );
	default: return 0;
	}
}

//
// Filling
//

int Filling1b( ImageBuffer buf, Region r, Color1b c ){
	return Filling( buf, r, ( unsigned char * )buf.img, c );
}

int Filling3b( ImageBuffer buf, Region r, Color3b c ){
	return Filling( buf, r, ( unsigned char * )buf.img, c );
}

int Filling4b( ImageBuffer buf, Region r, Color4b c ){
	return Filling( buf, r, ( unsigned char * )buf.img, c );
}

int Filling1i( ImageBuffer buf, Region r, Color1i c ){
	return Filling( buf, r, ( int * )buf.img, c );
}

int Filling3i( ImageBuffer buf, Region r, Color3i c ){
	return Filling( buf, r, ( int * )buf.img, c );
}

int Filling4i( ImageBuffer buf, Region r, Color4i c ){
	return Filling( buf, r, ( int * )buf.img, c );
}

int Filling1f( ImageBuffer buf, Region r, Color1f c ){
	return Filling( buf, r, ( float * )buf.img, c );
}

int Filling3f( ImageBuffer buf, Region r, Color3f c ){
	return Filling( buf, r, ( float * )buf.img, c );
}

int Filling4f( ImageBuffer buf, Region r, Color4f c ){
	return Filling( buf, r, ( float * )buf.img, c );
}

int Filling1d( ImageBuffer buf, Region r, Color1d c ){
	return Filling( buf, r, ( double * )buf.img, c );
}

int Filling3d( ImageBuffer buf, Region r, Color3d c ){
	return Filling( buf, r, ( double * )buf.img, c );
}

int Filling4d( ImageBuffer buf, Region r, Color4d c ){
	return Filling( buf, r, ( double * )buf.img, c );
}

//
// Iamge Selection
//

int Selectb( ImageBuffer buf, Selected sel, ImagePoint p, float thres ){
	switch( buf.color ){
	case 1: return Select1( buf, sel, p, thres, ( unsigned char * )buf.img );
	case 3: return Select3( buf, sel, p, thres, ( unsigned char * )buf.img );
	case 4: return Select4( buf, sel, p, thres, ( unsigned char * )buf.img );
	default: return 0;
	}
}

int Selecti( ImageBuffer buf, Selected sel, ImagePoint p, float thres ){
	switch( buf.color ){
	case 1: return Select1( buf, sel, p, thres, ( int * )buf.img );
	case 3: return Select3( buf, sel, p, thres, ( int * )buf.img );
	case 4: return Select4( buf, sel, p, thres, ( int * )buf.img );
	default: return 0;
	}
}

int Selectf( ImageBuffer buf, Selected sel, ImagePoint p, float thres ){
	switch( buf.color ){
	case 1: return Select1( buf, sel, p, thres, ( float * )buf.img );
	case 3: return Select3( buf, sel, p, thres, ( float * )buf.img );
	case 4: return Select4( buf, sel, p, thres, ( float * )buf.img );
	default: return 0;
	}
}

int Selectd( ImageBuffer buf, Selected sel, ImagePoint p, float thres ){
	switch( buf.color ){
	case 1: return Select1( buf, sel, p, thres, ( double * )buf.img );
	case 3: return Select3( buf, sel, p, thres, ( double * )buf.img );
	case 4: return Select4( buf, sel, p, thres, ( double * )buf.img );
	default: return 0;
	}
}

//
// Convolution
//

int Convolutionb( ImageBuffer buf, int size, float *f ){
	return Convolution( buf, size, f, ( BYTE * )buf.img );
}

int Convolutioni( ImageBuffer buf, int size, float *f ){
	return Convolution( buf, size, f, ( int * )buf.img );
}

int Convolutionf( ImageBuffer buf, int size, float *f ){
	return Convolution( buf, size, f, ( float * )buf.img );
}

int Convolutiond( ImageBuffer buf, int size, float *f ){
	return Convolution( buf, size, f, ( double * )buf.img );
}

//
// Gaussian Filtering
//

int GaussFilterb( ImageBuffer buf, int win ){
	return GaussFilter( buf, win, ( BYTE * )buf.img );
}

int GaussFilteri( ImageBuffer buf, int win ){
	return GaussFilter( buf, win, ( int * )buf.img );
}

int GaussFilterf( ImageBuffer buf, int win ){
	return GaussFilter( buf, win, ( float * )buf.img );
}

int GaussFilterd( ImageBuffer buf, int win ){
	return GaussFilter( buf, win, ( double * )buf.img );
}

//
// Scaling
//

int Scaleb( ImageBuffer src, ImageBuffer tar ){
	if( tar.img == NULL ) return 0;
	return Scale( src, ( BYTE* )tar.img, tar.w, tar.h, tar.stride );
}

int Scalei( ImageBuffer src, ImageBuffer tar ){
	if( tar.img == NULL ) return 0;
	return Scale( src, ( int* )tar.img, tar.w, tar.h, tar.stride );
}

int Scalef( ImageBuffer src, ImageBuffer tar ){
	if( tar.img == NULL ) return 0;
	return Scale( src, ( float* )tar.img, tar.w, tar.h, tar.stride );
}

int Scaled( ImageBuffer src, ImageBuffer tar ){
	if( tar.img == NULL ) return 0;
	return Scale( src, ( double* )tar.img, tar.w, tar.h, tar.stride );
}

//
// Dithering
//

void Ditherb( ImageBuffer src, unsigned char *t ){
	Dither( src, t, ( unsigned char )0xff, ( unsigned char )0x00 );
}

void Ditheri( ImageBuffer src, int *t, int max, int min ){
	Dither( src, t, max, min );
}

void Ditherf( ImageBuffer src, float *t, float max, float min ){
	Dither( src, t, max, min );
}

void Ditherd( ImageBuffer src, double *t, double max, double min ){
	Dither( src, t, max, min );
}

//
// Color Converting
//

int RGBtoNTSCi( ImageBuffer img ){
	static const int m[3][3] = {
		{ 299, 587, 114 },
		{ 596, -274, -322 },
		{ 211, -523, 312 }
	};
	int temp[3], *bit, *bits, i, j;
	if( img.color < 3 ) return 0;

	bits = ( int* )img.img;
	for( i = 0; i < img.h; i++ ){
		bit = ( int* )bits + i * img.stride;
		for( j = 0; j < img.w; j++, bit += img.color ){
			temp[0] = m[0][0] * bit[0] + m[0][1] * bit[1] + m[0][2] * bit[2];
			temp[1] = m[1][0] * bit[0] + m[1][1] * bit[1] + m[1][2] * bit[2];
			temp[2] = m[2][0] * bit[0] + m[2][1] * bit[1] + m[2][2] * bit[2];
			bit[0] = temp[0] / 1000;
			bit[1] = temp[1] / 1000;
			bit[2] = temp[2] / 1000;
		}
	}
	return 1;
}

int RGBtoNTSCf( ImageBuffer img ){
	static const float m[3][3] = {
		{ 0.299f, 0.587f, 0.114f },
		{ 0.596f, -0.274f, -0.322f },
		{ 0.211f, -0.523f, 0.312f }
	};
	float temp[3], *bit, *bits;
	int i, j;
	if( img.color < 3 ) return 0;

	bits = ( float* )img.img;
	for( i = 0; i < img.h; i++ ){
		bit = ( float* )bits + i * img.stride;
		for( j = 0; j < img.w; j++, bit += img.color ){
			temp[0] = m[0][0] * bit[0] + m[0][1] * bit[1] + m[0][2] * bit[2];
			temp[1] = m[1][0] * bit[0] + m[1][1] * bit[1] + m[1][2] * bit[2];
			temp[2] = m[2][0] * bit[0] + m[2][1] * bit[1] + m[2][2] * bit[2];
			bit[0] = temp[0];
			bit[1] = temp[1];
			bit[2] = temp[2];
		}
	}
	return 1;
}

int RGBtoNTSCd( ImageBuffer img ){
	static const double m[3][3] = {
		{ 0.299, 0.587, 0.114 },
		{ 0.596, -0.274, -0.322 },
		{ 0.211, -0.523, 0.312 }
	};
	double temp[3], *bit, *bits;
	int i, j;
	if( img.color < 3 ) return 0;

	bits = ( double* )img.img;
	for( i = 0; i < img.h; i++ ){
		bit = ( double* )bits + i * img.stride;
		for( j = 0; j < img.w; j++, bit += img.color ){
			temp[0] = m[0][0] * bit[0] + m[0][1] * bit[1] + m[0][2] * bit[2];
			temp[1] = m[1][0] * bit[0] + m[1][1] * bit[1] + m[1][2] * bit[2];
			temp[2] = m[2][0] * bit[0] + m[2][1] * bit[1] + m[2][2] * bit[2];
			bit[0] = temp[0];
			bit[1] = temp[1];
			bit[2] = temp[2];
		}
	}
	return 1;
}

int NTSCtoRGBi( ImageBuffer img ){
	static const int m[3][3] = {
		{ 1000, 956, 621 },
		{ 1000, -272, -647 },
		{ 1000, -1106, 1703 }
	};
	int temp[3], *bit, *bits, i, j;
	if( img.color < 3 ) return 0;

	bits = ( int* )img.img;
	for( i = 0; i < img.h; i++ ){
		bit = ( int* )bits + i * img.stride;
		for( j = 0; j < img.w; j++, bit += img.color ){
			temp[0] = m[0][0] * bit[0] + m[0][1] * bit[1] + m[0][2] * bit[2];
			temp[1] = m[1][0] * bit[0] + m[1][1] * bit[1] + m[1][2] * bit[2];
			temp[2] = m[2][0] * bit[0] + m[2][1] * bit[1] + m[2][2] * bit[2];
			bit[0] = temp[0] / 1000;
			bit[1] = temp[1] / 1000;
			bit[2] = temp[2] / 1000;
		}
	}
	return 1;
}

int NTSCtoRGBf( ImageBuffer img ){
	static const float m[3][3] = {
		{ 1.000f, 0.956f, 0.621f },
		{ 1.000f, -0.272f, -0.647f },
		{ 1.000f, -1.106f, 1.703f }
	};
	float temp[3], *bit, *bits;
	int i, j;
	if( img.color < 3 ) return 0;

	bits = ( float* )img.img;
	for( i = 0; i < img.h; i++ ){
		bit = ( float* )bits + i * img.stride;
		for( j = 0; j < img.w; j++, bit += img.color ){
			temp[0] = m[0][0] * bit[0] + m[0][1] * bit[1] + m[0][2] * bit[2];
			temp[1] = m[1][0] * bit[0] + m[1][1] * bit[1] + m[1][2] * bit[2];
			temp[2] = m[2][0] * bit[0] + m[2][1] * bit[1] + m[2][2] * bit[2];
			bit[0] = temp[0];
			bit[1] = temp[1];
			bit[2] = temp[2];
		}
	}
	return 1;
}

int NTSCtoRGBd( ImageBuffer img ){
	static const double m[3][3] = {
		{ 1.000, 0.956, 0.621 },
		{ 1.000, -0.272, -0.647 },
		{ 1.000, -1.106, 1.703 }
	};
	double temp[3], *bit, *bits;
	int i, j;
	if( img.color < 3 ) return 0;

	bits = ( double* )img.img;
	for( i = 0; i < img.h; i++ ){
		bit = ( double* )bits + i * img.stride;
		for( j = 0; j < img.w; j++, bit += img.color ){
			temp[0] = m[0][0] * bit[0] + m[0][1] * bit[1] + m[0][2] * bit[2];
			temp[1] = m[1][0] * bit[0] + m[1][1] * bit[1] + m[1][2] * bit[2];
			temp[2] = m[2][0] * bit[0] + m[2][1] * bit[1] + m[2][2] * bit[2];
			bit[0] = temp[0];
			bit[1] = temp[1];
			bit[2] = temp[2];
		}
	}
	return 1;
}