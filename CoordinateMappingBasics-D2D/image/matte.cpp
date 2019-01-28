#include "imgproc.h"
#include "linear.h"
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <stdio.h>

//
// define and macro
//

#define UNKNOWN		DBL_MAX
#define IsCorner( i, j, w, h )	\
	( ( i == 0 && j == 0 ) || ( i == 0 && j == w - 1 ) || ( i == h - 1 && j == w - 1 ) || ( i == h - 1 && j == 0 ) )
#define IsBoundary( i, j, w, h )	\
	( i == 0 || j == 0 || i == h - 1 || j == w - 1 )

//
// inline functions
//

template< typename T > T SQR( T x ){
	return x * x;
}

template< typename T > BYTE CheckColor( T x ){
	if( x > 254.0 ) return 0xff;
	else if( x < 0.5 ) return 0x00;
	else return ( BYTE )x;
}

//
// sparse matrix structure
//

typedef struct _node{
	double x;
	int i, j;
} Node;

//
// function decleration
//

void CompMeanVar( double *u, double *s, BYTE *img, int w, int h, int stride, int i, int j );
void CompMeanVar( double *u, double *s, BYTE *img, int w, int h, int stride );
int BuildL( Node *L, double *u, double *s, BYTE *img, int w, int h, int stride, double e );

//
// functions definition
//

#define AddMean	u[0] += temp[0];	u[1] += temp[1];	u[2] += temp[2];
#define AddVar	\
	s[0] += temp[0] * temp[0];	s[1] += temp[0] * temp[1];	s[2] += temp[0] * temp[2];	\
	s[3] += temp[1] * temp[0];	s[4] += temp[1] * temp[1];	s[5] += temp[1] * temp[2];	\
	s[6] += temp[2] * temp[0];	s[7] += temp[2] * temp[1];	s[8] += temp[2] * temp[2];

void CompMeanVar( double *u, double *s, BYTE *bit, int w, int h, int stride, int i, int j ){
	int n = 0;
	double temp[3], p;

	// mid
	temp[0] = ( double )bit[0];
	temp[1] = ( double )bit[1];
	temp[2] = ( double )bit[2];
	AddMean;	AddVar;	n++;
	if( i > 0 ){		// down
		temp[0] = ( double )bit[ -stride + 0 ];
		temp[1] = ( double )bit[ -stride + 1 ];
		temp[2] = ( double )bit[ -stride + 2 ];
		AddMean;	AddVar;	n++;
	}
	if( i < h - 1 ){	// top
		temp[0] = ( double )bit[ stride + 0 ];
		temp[1] = ( double )bit[ stride + 1 ];
		temp[2] = ( double )bit[ stride + 2 ];
		AddMean;	AddVar;	n++;
	}
	if( j > 0 ){			// left
		temp[0] = ( double )bit[ -3 + 0 ];
		temp[1] = ( double )bit[ -3 + 1 ];
		temp[2] = ( double )bit[ -3 + 2 ];
		AddMean;	AddVar;	n++;
		if( i > 0 ){		// down
			temp[0] = ( double )bit[ -stride - 3 + 0 ];
			temp[1] = ( double )bit[ -stride - 3 + 1 ];
			temp[2] = ( double )bit[ -stride - 3 + 2 ];
			AddMean;	AddVar;	n++;
		}
		if( i < h - 1 ){	// top
			temp[0] = ( double )bit[ stride - 3 + 0 ];
			temp[1] = ( double )bit[ stride - 3 + 1 ];
			temp[2] = ( double )bit[ stride - 3 + 2 ];
			AddMean;	AddVar;	n++;
		}
	}
	if( j < w - 1 ){		// right
		temp[0] = ( double )bit[ 3 + 0 ];
		temp[1] = ( double )bit[ 3 + 1 ];
		temp[2] = ( double )bit[ 3 + 2 ];
		AddMean;	AddVar;	n++;
		if( i > 0 ){		// down
			temp[0] = ( double )bit[ -stride + 3 + 0 ];
			temp[1] = ( double )bit[ -stride + 3 + 1 ];
			temp[2] = ( double )bit[ -stride + 3 + 2 ];
			AddMean;	AddVar;	n++;
		}
		if( i < h - 1 ){	// top
			temp[0] = ( double )bit[ stride + 3 + 0 ];
			temp[1] = ( double )bit[ stride + 3 + 1 ];
			temp[2] = ( double )bit[ stride + 3 + 2 ];
			AddMean;	AddVar;	n++;
		}
	}
	// normalize
	p = 1.0 / ( double )n;
	u[0] *= p;	u[1] *= p;	u[2] *= p;
	s[0] *= p;	s[1] *= p;	s[2] *= p;	s[3] *= p;	s[4] *= p;	s[5] *= p;	s[6] *= p;	s[7] *= p;	s[8] *= p;	
	s[0] -= u[0] * u[0];	s[1] -= u[0] * u[1];	s[2] -= u[0] * u[2];
	s[3] -= u[1] * u[0];	s[4] -= u[1] * u[1];	s[5] -= u[1] * u[2];
	s[6] -= u[2] * u[0];	s[7] -= u[2] * u[1];	s[8] -= u[2] * u[2];
}

void CompMeanVar( double *u, double *v, BYTE *img, int w, int h, int stride ){
	double *pu, *pv;
	BYTE *bit;
	int i, j;

	pu = u;	pv = v;
	for( i = 0; i < h; i++ ){
		bit = img + i * stride;
		for( j = 0; j < w; j++, bit += 3, pu += 3, pv += 9 ){
			CompMeanVar( pu, pv, bit, w, h, stride, i, j );
		}
	}
}

int BuildL( Node *L, double *u, double *v, BYTE *img, int w, int h, int stride, double e ){
	double *pu, *pv, temp[12];
	BYTE *biti, *bitj;
	int i, k, ii, ij, ji, jj, ki, kj, n, bottom, top, left, right;
	void *mem = malloc( 9 * sizeof( double ) );

	k = 0;
	// node in row (i)
	for( ii = 0; ii < h; ii++ ){
		biti = img + ii * stride;
		for( ij = 0; ij < w; ij++, biti += 3 ){
			// node in col (j) that must has a window k include i and j
			for( ji = max( ii - 2, 0 ); ji < h && ji - ii < 3; ji++ ){
				bitj = img + ji * stride + ( max( ij - 2, 0 ) ) * 3;
				for( jj = max( ij - 2, 0 ); jj < w && jj - ij < 3; jj++, bitj += 3 ){
					if( ii * w + ij > ji * w + jj ) continue;	// only compute upper triangle of L
					// set L
					L[k].i = ii * w + ij;
					L[k].j = ji * w + jj;
					L[k].x = 0.0;
					// find possible position of k
					bottom = max( min( ii, ji ) - 1, 0 );
					top = min( max( ii, ji ) + 1, h - 1 );
					left = max( min( ij, jj ) - 1, 0 );
					right = min( max( ij, jj ) + 1, w - 1 );
					// for all possible window k
					for( ki = bottom; ki <= top; ki++ ){
						pu = u + ( ki * w + left ) * 3;			// mean of window k
						pv = v + ( ki * w + left ) * 9;			// var of window k
						for( kj = left; kj <= right; kj++, pu += 3, pv += 9 ){
							// if k include i and j
							if( abs( ki - ii ) < 2 && abs( ki - ji ) < 2 && abs( kj - ij ) < 2 && abs( kj - jj ) < 2 ){
								if( IsCorner( ki, kj, w, h ) ) n = 4;
								else if( IsBoundary( ki, kj, w, h ) ) n = 6;
								else n = 9;
								temp[9] = e / ( double )n;
								for( i = 0; i < 9; i++ ) temp[i] = pv[i];
								temp[0] += temp[9];
								temp[4] += temp[9];
								temp[8] += temp[9];
								inverse( temp, 3, mem );
								temp[9] = ( double )bitj[0] - pu[0];
								temp[10] = ( double )bitj[1] - pu[1];
								temp[11] = ( double )bitj[2] - pu[2];
								temp[0] = temp[0] * temp[9] + temp[1] * temp[10] + temp[2] * temp[11];
								temp[1] = temp[3] * temp[9] + temp[4] * temp[10] + temp[5] * temp[11];
								temp[2] = temp[6] * temp[9] + temp[7] * temp[10] + temp[8] * temp[11];
								temp[9] = ( double )biti[0] - pu[0];
								temp[10] = ( double )biti[1] - pu[1];
								temp[11] = ( double )biti[2] - pu[2];
								temp[0] = temp[0] * temp[9] + temp[1] * temp[10] + temp[2] * temp[11];
								L[k].x -= ( 1.0 + temp[0] ) / ( double )n;
								if( L[k].i == L[k].j ) L[k].x += 1.0;
							}
						}
					}
					k++;	// next element of L
				}
			}
		}
	}
	free( mem );
	return k;
}

int Matte( BYTE *img, BYTE *alpha, int w, int h, int stride, int f[], int nf, int b[], int nb ){
	double *u, *v, *c, *v2, *pv, *pc, d, e, temp, temp2, r = 0.2;
	Node *L;
	BYTE *bit;
	int n, i, j, k, *p;
	void *mem;

	L = ( Node* )malloc( 13 * w * h * sizeof( Node ) );
	mem = calloc( 12 * w * h, sizeof( double ) );
	u = ( double* )mem;		// 3wh
	v = u + 3 * w * h;		// 9wh
	c = u;
	v2 = v + w * h;

	if( L == NULL || mem == NULL ) return 1;

	CompMeanVar( u, v, img, w, h, stride );
	n = BuildL( L, u, v, img, w, h, stride, 1.0 );

	// constraint
	j = w * h;
	for( i = 0; i < j; i++ ) c[i] = UNKNOWN;
	p = f;
	for( i = 0; i < nf; i++, p += 2 ) c[ p[1] * w + p[0] ] = 255.0;
	p = b;
	for( i = 0; i < nb; i++, p += 2 ) c[ p[1] * w + p[0] ] = 0.0;
	// initial guess
	pv = v;	pc = c;
	for( i = 0; i < h; i++ ) {
		for( j = 0; j < w; j++, pv++, pc++ ){
			if( pc[0] != UNKNOWN ) pv[0] = pc[0];
			else{
				// defult to foreground
				temp = DBL_MAX;	pv[0] = 255.0;
				// nearest distance to foreground
				p = f;
				for( k = 0; k < nf; k++, p += 2 ){
					e = p[1] - i;	d = p[0] - j;
					d = SQR( e ) + SQR( d );
					if( d < temp ) temp = d;
				}
				// nearest distance to background
				p = b;
				for( k = 0; k < nb; k++, p += 2 ){
					e = p[1] - i;	d = p[0] - j;
					d = SQR( e ) + SQR( d );
					if( d < temp ) {
						// nearest to background
						pv[0] = 0.0;	break;
					}
				}
			}
		}
	}
/**/
	// solve normalized cut
	do{
		// initialize
		d = e = 0.0;
		for( i = 0, j = w * h; i < j; i++ ) v2[i] = 0.0;
		// v2 = Av
		for( i = 0; i < n; i++ ){
			if( L[i].i != L[i].j ) {
				v2[ L[i].j ] += L[i].x * v[ L[i].i ];
				v2[ L[i].i ] += L[i].x * v[ L[i].j ];
			}
			else v2[ L[i].i ] += L[i].x * v[ L[i].j ];
		}
		for( i = 0, j = w * h; i < j; i++ ){
			e += v[i] * v2[i];						// energy
			if( c[i] != UNKNOWN ) continue;
			temp = v[i] - r * v2[i];				// gradient descent
			if( temp < 0.0 ) temp = 0.0;
			else if( temp > 255.0 ) temp = 255.0;
			temp2 = SQR( temp - v[i] );
			v[i] = temp;							// update
			if( d < temp2 ) d = temp2;				// max difference
		}
		printf( "%lf\t%lf\n", e, d );
	} while( d > 0.01 );
/**/
	// output
	pv = v;
	for( i = 0; i < h; i++ ){
		bit = alpha + i * stride;
		for( j = 0; j < w; j++, bit += 3, pv++ ){
			bit[0] = bit[1] = bit[2] = CheckColor( pv[0] );
		}
	}

	free( mem );
	free( L );
	return 0;
}

int Matte( BYTE *img, BYTE *scrib, BYTE *alpha, int w, int h, int stride ){
	int i, j, nf, nb, *f, *b, *pf, *pb;
	BYTE *bit;

	// counting
	nb = nf = 0;
	for( i = 0; i < h; i++ ){
		bit = scrib + i * stride;
		for( j = 0; j < w; j++, bit += 3 ){
			if( bit[0] == 0xff && bit[1] == 0x00 && bit[2] == 0x00 ) nb++;
			else if( bit[0] == 0x00 && bit[1] == 0x00 && bit[2] == 0xff ) nf++;
		}
	}

	// memory pool
	pf = f = ( int* )malloc( nf * 2 * sizeof( int ) );
	pb = b = ( int* )malloc( nb * 2 * sizeof( int ) );

	// setting
	for( i = 0; i < h; i++ ){
		bit = scrib + i * stride;
		for( j = 0; j < w; j++, bit += 3 ){
			if( bit[0] == 0xff && bit[1] == 0x00 && bit[2] == 0x00 ) {
				pb[0] = j; pb[1] = i; pb += 2;
				bit[0] = bit[1] = bit[2] = 0x00;
			}
			else if( bit[0] == 0x00 && bit[1] == 0x00 && bit[2] == 0xff ) {
				pf[0] = j; pf[1] = i; pf += 2;
				bit[0] = bit[1] = bit[2] = 0xff;
			}
			else 
				bit[0] = bit[1] = bit[2] = 0x77;
		}
	}
	// Matting
	i = Matte( img, scrib, w, h, stride, f, nf, b, nb );

	free( f );
	free( b );
	return i;
}

/**/