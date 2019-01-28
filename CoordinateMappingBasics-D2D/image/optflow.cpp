#include "vision.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "linear2.h"

//
// macro
//

#define SQR( x )			( (x) * (x) )
#define IsSeq( x )			( (x).n > 1 )

//
// algorithm template
//

template< typename T >
void HornOptFlow( OpticalFlow *opt, T *t0, T *t1, int w, int h, int color, int stride ){
	const int w1 = w - 1, h1 = h - 1;
	const float wt1 = 0.1666666666666666667f;	// 1 / 6
	const float wt2 = 0.0833333333333333333f;	// 1 / 12
	const float wt3 = 0.25f / ( float )color;	// 1 / 4 / color
	const float a2 = SQR( opt->alpha );			// a2
	float temp[2][2][2];
	float *u, *v, *pu, *pv, *pu0, *pv0, *ex, *ey, *et, *pex, *pey, *pet, 
		u0, v0;
	T *pt0, *pt1;
	int i, j, k, jmp;

	// initialize
	ex = opt->buf;
	ey = ex + w * h;
	et = ey + w * h;
	u = et + w * h;
	v = u + w * h;

	// compute partial derivatives
	for( i = 1; i < h1; i++ ){
		jmp = i * w + 1;
		pt0 = t0 + i * stride + color;	// pointer to t0 image
		pt1 = t1 + i * stride + color;	// pointer to t1 image
		pex = ex + jmp;			// pointer to Ex
		pey = ey + jmp;			// pointer to Ey
		pet = et + jmp;			// pointer to Et
		for( j = 1; j < w1; j++, pex++, pey++, pet++ ){
			pex[0] = pey[0] = pet[0] = 0.0f;
			for( k = 0; k < color; k++, pt0++, pt1++ ){
				temp[0][0][0] = ( float )pt0[0];
				temp[1][0][0] = ( float )pt0[ color ];
				temp[0][1][0] = ( float )pt0[ stride ];
				temp[1][1][0] = ( float )pt0[ stride + color ];
				temp[0][0][1] = ( float )pt1[0];
				temp[1][0][1] = ( float )pt1[ color ];
				temp[0][1][1] = ( float )pt1[ stride ];
				temp[1][1][1] = ( float )pt1[ stride + color ];
				*pex += temp[1][0][0] + temp[1][0][1] + temp[1][1][0] + temp[1][1][1] 
					- temp[0][0][0] - temp[0][0][1] - temp[0][1][0] - temp[0][1][1];
				*pey += temp[0][1][0] + temp[0][1][1] + temp[1][1][0] + temp[1][1][1] 
					- temp[0][0][0] - temp[0][0][1] - temp[1][0][0] - temp[1][0][1];
				*pet += temp[0][0][1] + temp[0][1][1] + temp[1][0][1] + temp[1][1][1] 
					- temp[0][0][0] - temp[0][1][0] - temp[1][0][0] - temp[1][1][0];
			}
			pex[0] *= wt3; 
			pey[0] *= wt3; 
			pet[0] *= wt3; 
		}
	}
	// compute optical flow
	for( k = 0; k < opt->iter; k++ ){
		for( i = 1; i < h1; i++ ){
			jmp = i * w + 1;
			pu = u + jmp;
			pv = v + jmp;
			pu0 = opt->u + jmp;		// pointer to optical flow u
			pv0 = opt->v + jmp;		// pointer to optical flow v
			pex = ex + jmp;			// pointer to Ex
			pey = ey + jmp;			// pointer to Ey
			pet = et + jmp;			// pointer to Et
			for( j = 1; j < w1; j++, pu++, pv++, pu0++, pv0++, pex++, pey++, pet++ ){
				u0 = wt1 * ( pu0[ -1 ] + pu0[ 1 ] + pu0[ -w ] + pu0[ w ] ) +
					wt2 * ( pu0[ -w - 1 ] + pu0[ w - 1 ] + pu0[ -w + 1 ] + pu0[ w + 1 ] );
				v0 = wt1 * ( pv0[ -1 ] + pv0[ 1 ] + pv0[ -w ] + pv0[ w ] ) +
					wt2 * ( pv0[ -w - 1 ] + pv0[ w - 1 ] + pv0[ -w + 1 ] + pv0[ w + 1 ] );
				***temp = ( pex[0] * u0 + pey[0] * v0 + pet[0] ) / ( a2 + SQR( pex[0] ) + SQR( pey[0] ) );
				pu[0] = u0 - pex[0] * ***temp;
				pv[0] = v0 - pey[0] * ***temp;
			}
			pu[0] = pu[-1];	u[ jmp - 1 ] = u[ jmp ];
			pv[0] = pv[-1];	v[ jmp - 1 ] = v[ jmp ];
		}
		memcpy( u, u + w, w * sizeof( float ) );
		memcpy( u + i * w, u + ( i - 1 ) * w, w * sizeof( float ) );
		memcpy( v, v + w, w * sizeof( float ) );
		memcpy( v + i * w, v + ( i - 1 ) * w, w * sizeof( float ) );
		memcpy( opt->u, u, w * h * sizeof( float ) );
		memcpy( opt->v, v, w * h * sizeof( float ) );
	}
}

// Lucas / Kanade's Algorithm

template < typename T >
void LKOptFlow( OpticalFlow *opt, T *t1, T *t0, int w, int h, int color, int stride ){
	const int w1 = w - 1;
	const int h1 = h - 1;
	const float e2 = SQR( opt->alpha );
	float temp[2][3], p[6], dp[6], he[6][6], 
		*sd, *psd, *dx, *dy, *pdx, *pdy, *pu, *pv; 
	float x, y, u1, v1, u2, v2, e;
	T *bit1, *bit2;
	int index[6], 
		i, j, k, t, jmp, xi, yi, n;

	// initialize memory pool
	dx = opt->buf;				// h * w * color
	dy = dx + color * w * h;	// h * w * color
	sd = dy + color * w * h;	// steepest descent( color * 6 )
	memset( p, 0x00, 6 * sizeof( float ) );

	// compute gradient image by Sobel operation
	for( i = 1; i < h1; i++ ){
		// jump
		jmp = ( i * w + 1 ) * color;
		bit1 = t0 + i * stride + color;	
		pdx = dx + jmp;		pdy = dy + jmp;
		for( j = 1; j < w1; j++ ){
			for( k = 0; k < color; k++, bit1++, pdx++, pdy++ ){
				// x-sobel
				temp[0][0] = ( float )bit1[ stride + color ];
				temp[0][1] = ( float )bit1[ color ];
				temp[0][2] = ( float )bit1[ -stride + color ];
				temp[1][0] = ( float )bit1[ stride - color ];
				temp[1][1] = ( float )bit1[ -color ];
				temp[1][2] = ( float )bit1[ -stride - color ];
				*pdx = 0.125f * ( temp[0][0] + temp[0][1] + temp[0][1] + temp[0][2] - 
					temp[1][0] - temp[1][1] - temp[1][1] - temp[1][2] );
				// y-sobel
				temp[0][1] = ( float )bit1[ stride ];
				temp[1][1] = ( float )bit1[ -stride ];
				*pdy = 0.125f * ( temp[0][0] + temp[0][1] + temp[0][1] + temp[1][0] - 
					temp[0][2] - temp[1][1] - temp[1][1] - temp[1][2] );
			}
		}
	}
	// LK-Algotithm
	for( t = 0; t < opt->iter; t++ ){
		// zerolize
		memset( he[0], 0x00, 36 * sizeof( float ) );
		memset( dp, 0x00, 6 * sizeof( float ) );
		e = 0.0f;	n = 0;
		// compute steepest descent
		for( i = 1; i < h1; i++ ){
			bit2 = t1 + i * stride + color;	
			for( j = 1; j < w1; j++, bit2 += color ){
				// compute coordinate after warping
				x = ( 1.0f + p[0] ) * ( float )j + p[2] * ( float )i + p[4];
				y = p[1] * ( float )j + ( 1.0f + p[3] ) * ( float )i + p[5];
				xi = ( int )x;	yi = ( int )y;
				// check point after warping
				if( xi < 1 || yi < 1 || xi >= w1 - 1 || yi >= h1 - 1 ) continue;
				u2 = ( float )fmod( x, 1.0f );	u1 = 1.0f - u2;	// weights of bilinear
				v2 = ( float )fmod( y, 1.0f );	v1 = 1.0f - v2;	// weights of bilinear
				psd = sd;
				jmp = ( yi * w + xi ) * color;					// jump of derivation
				for( k = 0; k < color; k++, psd += 6, jmp++ ){
					// compute steepest descent
					temp[0][0] = dx[ jmp ];
					temp[0][1] = dy[ jmp ];
					psd[0] = temp[0][0] * j;	psd[1] = temp[0][1] * j;
					psd[2] = temp[0][0] * i;	psd[3] = temp[0][1] * i;
					psd[4] = temp[0][0];		psd[5] = temp[0][1];
					// compute increment p by bilinear
					**temp = ( float )bit2[k] - (
						u1 * v1 * ( float )t0[ yi * stride + xi * color + k ] +
						u2 * v1 * ( float )t0[ yi * stride + ( xi + 1 ) * color + k ] +
						u1 * v2 * ( float )t0[ ( yi + 1 ) * stride + xi * color + k ] +
						u2 * v2 * ( float )t0[ ( yi + 1 ) * stride + ( xi + 1 ) * color + k ] ); 
					dp[0] += **temp * psd[0];	dp[1] += **temp * psd[1];	
					dp[2] += **temp * psd[2];	dp[3] += **temp * psd[3];
					dp[4] += **temp * psd[4];	dp[5] += **temp * psd[5];
					// compute Hessian matrix
					he[0][0] += psd[0] * psd[0];	he[0][1] += psd[0] * psd[1];	he[0][2] += psd[0] * psd[2];	
					he[0][3] += psd[0] * psd[3];	he[0][4] += psd[0] * psd[4];	he[0][5] += psd[0] * psd[5];
					he[1][1] += psd[1] * psd[1];	he[1][2] += psd[1] * psd[2];	he[1][3] += psd[1] * psd[3];	
					he[1][4] += psd[1] * psd[4];	he[1][5] += psd[1] * psd[5];	
					he[2][2] += psd[2] * psd[2];	he[2][3] += psd[2] * psd[3];	
					he[2][4] += psd[2] * psd[4];	he[2][5] += psd[2] * psd[5];
					he[3][3] += psd[3] * psd[3];	he[3][4] += psd[3] * psd[4];	he[3][5] += psd[3] * psd[5];
					he[4][4] += psd[4] * psd[4];	he[4][5] += psd[4] * psd[5];
					he[5][5] += psd[5] * psd[5];
					// adding error
					e += SQR( **temp );
					n++;
				}
			}
		}
		e /= ( float )n;
		if( e < e2 ) break;
		// complete hessian matrix
		he[1][0] = he[0][1];	he[2][0] = he[0][2];	he[3][0] = he[0][3];	he[4][0] = he[0][4];	he[5][0] = he[0][5];
		he[2][1] = he[1][2];	he[3][1] = he[1][3];	he[4][1] = he[1][4];	he[5][1] = he[1][5];
		he[3][2] = he[2][3];	he[4][2] = he[2][4];	he[5][2] = he[2][5];
		he[4][3] = he[3][4];	he[5][3] = he[3][5];
		he[5][4] = he[4][5];
		// LU-factorize
		LUfactf( he[0], index, 6, temp[0] );
		SolveLUf( he[0], index, dp, 6 );
		// update parameter
		p[0] += dp[0];
		p[1] += dp[1];
		p[2] += dp[2];
		p[3] += dp[3];
		p[4] += dp[4];
		p[5] += dp[5];
	}
	// save the optical flow
	pu = opt->u;	pv = opt->v;
	for( i = 0; i < h; i++ ){
		for( j = 0; j < w; j++, pu++, pv++ ){
			pu[0] = p[0] * ( float )j + p[2] * ( float )i + p[4];
			pv[0] = p[1] * ( float )j + p[3] * ( float )i + p[5];
		}
	}
}

// 
// functions
// 

// Memory

int NewOpticalFlow( OpticalFlow *opt, const SequenceBuffer seq, OptFlowAlgo algo ){
	switch( ( int )algo ){
	case OF_HORN_ALGO:
		// velocity image
		opt->u = ( float * )calloc( seq.w * seq.h * 2, sizeof( float ) );
		opt->v = opt->u + seq.w * seq.h;
		if( opt->u == NULL ) return 0;
		// calculating buffer
		opt->buf = ( float * )malloc( seq.w * seq.h * 5 * sizeof( float ) );
		if( opt->buf == NULL ) return 0;
		break;
	case OF_LK_ALGO:
		// velocity image
		opt->u = ( float * )calloc( seq.w * seq.h * 2, sizeof( float ) );
		opt->v = opt->u + seq.w * seq.h;
		if( opt->u == NULL ) return 0;
		// calculating buffer( gradient image )
		opt->buf = ( float * )malloc( seq.color * ( seq.w * seq.h * 2 + 6 ) * sizeof( float ) );
		break;
	default:
		return 0;
	}	
	return 1;
}

int DeleteOpticalFlow( OpticalFlow *opt ){
	if( opt->u == NULL || opt->v == NULL ) return 0;
	free( opt->u );
	if( opt->buf ) free( opt->buf );
	opt->u = opt->v = opt->buf = NULL;
	return 1;
}

// Horn Optical Flow

int HornOptFlowb( OpticalFlow *opt, const SequenceBuffer s ){
	if( !IsSeq( s ) ) return 0;
	HornOptFlow( opt, ( unsigned char * )s.seq[0], ( unsigned char * )s.seq[1], s.w, s.h, s.color, s.stride );
	return 1;
}

int HornOptFlowi( OpticalFlow *opt, const SequenceBuffer s ){
	if( !IsSeq( s ) ) return 0;
	HornOptFlow( opt, ( int * )s.seq[0], ( int * )s.seq[1], s.w, s.h, s.color, s.stride );
	return 1;
}

int HornOptFlowf( OpticalFlow *opt, const SequenceBuffer s ){
	if( !IsSeq( s ) ) return 0;
	HornOptFlow( opt, ( float * )s.seq[0], ( float * )s.seq[1], s.w, s.h, s.color, s.stride );
	return 1;
}

int HornOptFlowd( OpticalFlow *opt, const SequenceBuffer s ){
	if( !IsSeq( s ) ) return 0;
	HornOptFlow( opt, ( double * )s.seq[0], ( double * )s.seq[1], s.w, s.h, s.color, s.stride );
	return 1;
}

// Lucas Kanade's Optical Flow

int LKOptFlowb( OpticalFlow *opt, const SequenceBuffer s ){
	if( !IsSeq( s ) ) return 0;
	LKOptFlow( opt, ( unsigned char * )s.seq[0], ( unsigned char * )s.seq[1], s.w, s.h, s.color, s.stride );
	return 1;
}

int LKOptFlowi( OpticalFlow *opt, const SequenceBuffer s ){
	if( !IsSeq( s ) ) return 0;
	LKOptFlow( opt, ( int * )s.seq[0], ( int * )s.seq[1], s.w, s.h, s.color, s.stride );
	return 1;
}

int LKOptFlowf( OpticalFlow *opt, const SequenceBuffer s ){
	if( !IsSeq( s ) ) return 0;
	LKOptFlow( opt, ( float * )s.seq[0], ( float * )s.seq[1], s.w, s.h, s.color, s.stride );
	return 1;
}

int LKOptFlowd( OpticalFlow *opt, const SequenceBuffer s ){
	if( !IsSeq( s ) ) return 0;
	LKOptFlow( opt, ( double * )s.seq[0], ( double * )s.seq[1], s.w, s.h, s.color, s.stride );
	return 1;
}