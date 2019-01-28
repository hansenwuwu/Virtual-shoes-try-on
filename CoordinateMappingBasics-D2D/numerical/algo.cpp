#include "algo.h"
#include "linear.h"
#include <limits.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

// macro

#ifndef M_PI
#define M_PI	3.1415926535897932384626433832795
#endif
#define ERR		0.00001
#define SQR(x)	( (x) * (x) )

//
// K-Mean Clustering
//

void Kmean( int *data, int *cluster, int *domain, int *out, int n, int d , int c, 
		   bool init, void *mem ){
	int i, j, k, dist, min, it, *center, *cn;
	void *memory;

	memory = malloc( ( c + c * d ) * sizeof( int ) );
	center = ( int* )memory;
	cn = center + c * d;

	// initialize cluster center
	if( init ){
		for( j = 0; j < d; j++ ) {
			min = domain[j * 2 + 0];
			for( i = 0; i < c; i++ ) center[i * d + j] = rand() % ( domain[j * 2 + 1] - min ) + min;
		}
	}
	else for( i = 0; i < c * d; i++ ) center[i] = cluster[i];

	// clustering
	it = 0;
	do{
		it++;
		for( i = 0; i < c * d; i++ ) {
			cluster[i] = center[i];		// update cluster center
			center[i] = 0;				// initialize new cluster center
		}
		for( i = 0; i < c; i++ ) cn[i] = 0;				// initialize cluster counter
		for( i = 0; i < n; i++ ){		// for every data point
			// determine cluster
			min = INT_MAX;
			for( j = 0; j < c; j++ ){	// for every clustering center
				dist = 0;
				for( k = 0; k < d; k++ ) dist += SQR( data[i * d + k ] - cluster[j * d + k] );
				if( dist < min ) {
					min = dist;
					out[i] = j;			// determine cluster
				}
			}
			k = out[i];					// cluster
			for( j = 0; j < d; j++ ) center[ k * d + j ] += data[i * d + j];
			cn[k]++;
		}
		// average new center
		for( i = 0; i < c; i++ ){
			if( cn[i] == 0 ) {			// bad clustering center
				// re-randomize cluster center
				min = domain[j * 2 + 0];
				for( j = 0; j < d; j++ ) center[i * d + j] = rand() % ( domain[j * 2 + 1] - min ) + min;
			}
			// average
			else for( j = 0; j < d; j++ ) center[i * d + j] /= cn[i];
		}
		// calculate distance between new and old cluster centers
		for( i = 0; i < c; i++ ){	// for every cluster center
			if( abs( cluster[i * d + j] - center[i * d + j] ) > 1 ) break;
		}
	} while( i < c && it < 50 );

	free( memory );
}

void Kmean( double *data, double *cluster, double *domain, int *out, int n, int d, int c, 
		   bool init, void *mem ){
	int i, j, k, it, *cn;
	double dist, min, *center;
	void *memory;

	memory = malloc( c * d * sizeof( double ) + c * sizeof( int ) );
	center = ( double* )memory;
	cn = ( int* )( center + c * d );

	// generate cluster center
	if( init ){
		for( j = 0; j < d; j++ ) {
			min = domain[j * 2 + 0];
			for( i = 0; i < c; i++ ) 
				center[i * d + j] = ( double )rand() / ( double )RAND_MAX * ( domain[j * 2 + 1] - min ) + min;
		}
	}
	else for( i = 0; i < c * d; i++ ) center[i] = cluster[i];

	// clustering
	it = 0;
	do{
		it++;
		for( i = 0; i < c * d; i++ ) {
			cluster[i] = center[i];		// update cluster center
			center[i] = 0.0;				// initialize new cluster center
		}
		for( i = 0; i < c; i++ ) cn[i] = 0;				// initialize cluster counter
		for( i = 0; i < n; i++ ){		// for every data point
			// determine cluster
			min = DBL_MAX;
			for( j = 0; j < c; j++ ){	// for every clustering center
				dist = 0.0;
				for( k = 0; k < d; k++ ) dist += SQR( data[i * d + k ] - cluster[j * d + k] );
				if( dist < min ) {
					min = dist;
					out[i] = j;			// determine cluster
				}
			}
			k = out[i];					// cluster
			for( j = 0; j < d; j++ ) center[ k * d + j ] += data[i * d + j];
			cn[k]++;
		}
		// average new center
		for( i = 0; i < c; i++ ){
			if( cn[i] == 0 ) {			// bad clustering center
				// re-randomize cluster center
				min = domain[j * 2 + 0];
				for( j = 0; j < d; j++ ) 
					center[i * d + j] = ( double )rand() / ( double )RAND_MAX * ( domain[j * 2 + 1] - min ) + min;
			}
			// average
			else for( j = 0; j < d; j++ ) center[i * d + j] /= cn[i];
		}
		// calculate distance between new and old cluster centers
		for( i = 0; i < c; i++ ){	// for every cluster center
			if( fabs( cluster[i * d + j] - center[i * d + j] ) > ERR ) break;
		}
	} while( i < c && it < 50 );

	free( memory );
}

//
// Histogram Equalization
//

void Histogram( int *data, int b, int t, int n, int *H ){
	int i, size = t - b + 1;

	// memory setting
	if( H == NULL )	H = ( int* )calloc( size, sizeof( int ) );
	else for( i = 0; i < size; i++ ) H[i] = 0;

	for( i = 0; i < n; i++ ) H[ data[i] ]++;
	for( i = 1; i < size; i++ ) H[i] = H[ i - 1 ] + H[i];
	for( i = 0; i < n; i++ ) 
		data[i] = ( int )( ( double )t / ( double )n * ( double )H[ data[i] ] ) + b;
	
	if( H == NULL ) free( H );
}

void Histogram( int *data1, int *data2, int b, int t, 
			   int n1, int n2, int *H1, int *H2 ){
	int *H, i, j, low, up, jmp, temp, size = t - b + 1;
	float r;

	// memory allocate
	if( H1 == NULL ) H1 = ( int* )calloc( size, sizeof( int ) );
	else for( i = 0; i < size; i++ ) H1[i] = 0;
	if( H2 == NULL ) H2 = ( int* )calloc( size, sizeof( int ) );
	else for( i = 0; i < size; i++ ) H2[i] = 0;

	// calculate pdf
	for( i = 0; i < n1; i++ ) H1[ data1[i] - b ]++;
	for( i = 0; i < n2; i++ ) H2[ data2[i] - b ]++;

	// threshold	n * 0.2
	temp = ( int )( 0.2f * ( float )n1 );
	for( i = 0; i < size; i++ ) if( H1[i] > temp ) H1[i] = temp;
	temp = ( int )( 0.2f * ( float )n2 );
	for( i = 0; i < size; i++ ) if( H2[i] > temp ) H2[i] = temp;
	
	// smooth
	temp = size - 2;
	for( i = 2; i < temp; i++ ) {
		H1[i] = ( H1[i - 2] + 2 * H1[i - 1] + 4 * H1[i] + 2 * H1[i + 1] + H1[i + 2] ) / 9;
		H2[i] = ( H2[i - 2] + 2 * H2[i - 1] + 4 * H2[i] + 2 * H2[i + 1] + H2[i + 2] ) / 9;
	}
	// trnaslate pdf to cdf
	for( i = 1; i < size; i++ ) {
		H1[i] = H1[ i - 1 ] + H1[i];
		H2[i] = H2[ i - 1 ] + H2[i];
	}
	// ratio
	r = ( float )H1[ size - 1 ] / ( float )H2[ size - 1 ];
	// equal sample size
	if( n1 != n2 )
		for( i = 0; i < size; i++ ) H2[i] = ( int )( ( float )H2[i] * r );
	// domain of H2: [ low, up ]
	for( low = 0; H2[ low + 1 ] == 0 ; low++ );
	for( up = size - 1; H2[ up - 1 ] == H2[up]; up-- );
	H = H2 + low;
	size = up - low;
		
	// histogram transform
	for( i = 0; i < n1; i++ ){
		// binary search for H2 = H1
		temp = H1[ data1[i] - b ];
		jmp = j = size / 2;
		do{
			if( H[ j ] <= temp && H[ j + 1 ] >= temp )	// equal
				break;
			jmp = ( jmp + 1 ) / 2;
			if( H[j] < temp ) j += jmp;						// move right
			else j -= jmp;										// move left
		} while( jmp > 1 );
		// save value
		data1[i] = j + b + low;
	}
}

void Histogram( int *source, int *target, int *b, int *t, int d, int n1, int n2, int time ){
	int *H1, *H2, *xr1, *xr2, *x1, *x2, *p, 
		jmp, size, i, j, k, a, it, min, max;
	float f[2], shift = 2.0f * ( float )M_PI / ( float )time, angle, c, s;
	void *mem;

	// calculate size
	size = 0;
	for( i = 0; i < d; i++ ) {
		j = t[i] - b[i];
		size += SQR( j );
	}
	size = ( int )sqrt( ( double )size ) * 2;

	// memory pool
	mem = malloc( ( 2 * size + ( n1 + n2 ) * ( d + 1 ) ) * sizeof( int ) );
	H1 = ( int* )mem;				// size
	H2 = H1 + size;					// size
	xr1 = H2 + size;				// n1 * d
	xr2 = xr1 + n1 * d;				// n2 * d
	x1 = xr2 + n2 * d;				// n1
	x2 = x1 + n1;					// n2

	// iteration
	angle = 0.0;			// intialize
	for( it = 0; it < time; it++, angle += shift ){	
		// rotate at plain i(0) and a
		a = it % ( d - 1 ) + 1;
		c = cos( angle );	// cos
		s = sin( angle );	// sin

		// rotate source
		p = source;
		jmp = 0;
		for( i = 0; i < n1; i++ ){
			// copy
			for( j = 1; j < d; j++ ) xr1[ jmp + j ] = p[j];
			// rotate
			f[0] = ( float )p[0];
			f[1] = ( float )p[a];
			xr1[ jmp ] = ( int )( c * f[0] - s * f[1] );
			xr1[ jmp + a ] = ( int )( s * f[0] + c * f[1] );
			// next vector
			p += d;	jmp += d;
		}

		// roate target
		p = target;
		jmp = 0;
		for( i = 0; i < n2; i++ ){
			// copy
			for( j = 1; j < d; j++ ) xr2[ jmp + j ] = p[j];
			// rotate
			f[0] = ( float )p[0];
			f[1] = ( float )p[a];
			xr2[ jmp ] = ( int )( c * f[0] - s * f[1] );
			xr2[ jmp + a ] = ( int )( s * f[0] + c * f[1] );
			// next vector
			p += d;	jmp += d;
		}

		// project to each dimension
		for( k = 0; k < d; k++ ){	
			// project to first dimension
			max = 0;	min = INT_MAX;						// find domain range
			for( i = 0, j = k; i < n1; i++, j += d ) {		// source
				x1[i] = xr1[j];	
				if( x1[i] < min ) min = x1[i];
				if( x1[i] > max ) max = x1[i];
			}
			for( i = 0, j = k; i < n2; i++, j += d ) {		// target
				x2[i] = xr2[j];
				if( x2[i] < min ) min = x2[i];
				if( x2[i] > max ) max = x2[i];
			}
			// histogram transfer
			Histogram( x1, x2, min, max, n1, n2, H1, H2 );
			p = xr1 + k;
			for( i = 0; i < n1; i++, p += d ) *p = x1[i];	// save back
		}
		// rotate back
		p = source;
		jmp = 0;
		for( i = 0; i < n1; i++ ){
			for( j = 1; j < d; j++ ) p[j] = xr1[ jmp + j ] ;
			f[0] = ( float )xr1[ jmp ];
			f[1] = ( float )xr1[ jmp + a ];
			p[0] = ( int )( c * f[0] + s * f[1] + 0.5f );
			p[a] = ( int )( -s * f[0] + c * f[1] + 0.5f );
			jmp += d;
			p += d;
		}
	}

	// free memory pool
	free( mem );
}

double NormalizedCut( double *W, double *u, int n, int time ){
	double *A, *a, *v, max = 0.0, prev, temp;
	int i, j, k, it, *ind;
	void *mem;
	clock_t t;

	// memory allocate
	k = SQR( n );
	mem = malloc( ( k + n ) * sizeof( double ) + n * sizeof( int ) );
	A = ( double * )mem;		// ( n^2 + n ) / 2
	v = A + k;					// n
	ind = ( int* )( v + n );	// n

	// initialize
	for( i = 0; i < n; i++ ) u[i] = ( double )rand() / ( double )RAND_MAX;
	for( i = 0; i < k; i++ ) A[i] = W[i];
	for( i = 0; i < k; i += n + 1 ) printf( "%lf\n", A[i] );
	getchar();

	t = clock();
	for( k = 0; k < n; k++ ){
		// LU decomposite for x = A^-1y -> Ax = y, solve x
		LUfact( A, ind, n, v );
		// the i-th largest eigenvalue
		prev = DBL_MAX;
		for( it = 0; it < time; it++ ){
			if( fabs( max - prev ) < ERR ) break;
			// calculate v
/*			for( i = k; i < n; i++ ){
				v[i] = 0.0;
				a = A + i * n + k;
				for( j = k; j < n; j++ ) v[i] += u[j] * *( a++ );
			}*/
			for( i = 0; i < n; i++ ) v[i] = u[i];
			SolveLU( A, ind, v, n );
			// calculate u
			prev = max;
			max = 0.0;
			for( i = k; i < n; i++ ) {
				if( fabs( v[i] ) > fabs( max ) ) max = v[i];
			}
			if( fabs( max ) > ERR ) for( i = k; i < n; i++ ) u[i] = v[i] / max;
		}
		temp = max;
		// householder transform to A( k+1 ), Hu = e1
		LUcomb( A, ind, n );
		if( k < 1 ){
			max = 0.0;
			for( i = k; i < n; i++ ) max += SQR( u[i] );
			max = sqrt( max );
			for( i = k; i < n; i++ ) u[i] /= max;		// alpha = 1
			max = 1.0 - u[k];							// H
			u[k] -= 1.0;								// u1
			// p
			for( i = k; i < n; i++ ){
				v[i] = 0.0;
				a = A + i * n + k;
				for( j = k; j < n; j++ ) v[i] += *( a++ ) * u[j];
				v[i] /= max;	// Au / H
			}
			// K
			prev = 0.0;
			for( i = k; i < n; i++ ) prev += u[i] * v[i];
			prev /= max;		// 2K
			for( i = k; i < n; i++ ){
				a = A + i * n + k;
				for( j = k; j < n; j++ ){
					*( a++ ) += - v[i] * u[j] - v[j] * u[i] + prev * u[i] * u[j];
				}
			}
		}
		v[k] = temp;
	}
	printf( "find eigenvalue: %lf\n", ( double )( clock() - t ) / ( double )CLOCKS_PER_SEC );
	t = clock();

	j = 0; k = 0;	// j: smallest, k: second smallest
/*
	for( i = 0; i < n; i++ ) if( fabs( v[i] ) < fabs( v[j] ) ) j = i;
	for( i = 0; i < n; i++ ) if( i != j && fabs( v[i] ) < fabs( v[k] ) ) k = i;
/**/
	for( i = 0; i < n; i++ ) if( v[i] < v[j] ) j = i;
	for( i = 0; i < n; i++ ) if( i != j && v[i] < v[k] ) k = i;
/**/
	temp = v[j];	// eigenvalue

	temp = InvIter( W, u, temp, n, 50, v );

	// normalize
	max = 0.0;
//	for( i = 0; i < n; i++ ) max += SQR( u[i] );
//	max = sqrt( max );
	for( i = 0; i < n; i++ ) if( fabs( u[i] ) > fabs( max ) ) max = u[i];
	for( i = 0; i < n; i++ ) u[i] /= max;

	free( mem );
	return temp;
}