
#include "linear2.h"
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdlib.h>

// 
// Macro
//

#define SQR(x)	( (x) * (x) )

//
// Determinant
//

float detf( float *m, int n ){
	float ans = 1.0f;
	float ratio, max, temp;
	int interchange = 0;
	int i, j, k, n1 = n - 1;
	int row;

	for( i = 0; i < n1; i++ ){
		row = i;
		max = 0.0f;
		for( j = i; j < n; j++ ){
			if( fabs( m[j * n + i] ) > fabs( max ) ) {
				max = m[j * n + i];
				row = j;
			}
		}
		if( m[row * n + i] == 0.0f ) continue;	//all downside value is 0
		if( row != i ){
			for( j = 0; j < n; j++ ) {
				temp = m[i * n + j];
				m[i * n + j] = m[row * n + j];
				m[row * n + j] = temp;
			}
			interchange =! interchange;
		}
		for( j = i + 1; j < n; j++ ) {
			if( m[j * n + i] == 0.0f ) continue;
			ratio = m[j * n + i] / m[i * n + i];
			for( k = 0; k < n; k++ ) m[j * n + k] -= ( m[i * n + k] * ratio );
		}
	}
	for(i=0; i<n; i++) ans *= m[i * n + i];

	if( interchange ) return -ans;
	else return ans;
}

double detd( double *m, int n ){
	double ans = 1.0f;
	double ratio, max, temp;
	int interchange = 0;
	int i, j, k, n1 = n - 1;
	int row;

	for( i = 0; i < n1; i++ ){
		row = i;
		max = 0.0;
		for( j = i; j < n; j++ ){
			if( fabs( m[j * n + i] ) > fabs( max ) ) {
				max = m[j * n + i];
				row = j;
			}
		}
		if( m[row * n + i] == 0.0 ) continue;	//all downside value is 0
		if( row != i ){
			for( j = 0; j < n; j++ ) {
				temp = m[i * n + j];
				m[i * n + j] = m[row * n + j];
				m[row * n + j] = temp;
			}
			interchange =! interchange;
		}
		for( j = i + 1; j < n; j++ ) {
			if( m[j * n + i] == 0.0 ) continue;
			ratio = m[j * n + i] / m[i * n + i];
			for( k = 0; k < n; k++ ) m[j * n + k] -= ( m[i * n + k] * ratio );
		}
	}
	for(i=0; i<n; i++) ans *= m[i * n + i];

	if( interchange ) return -ans;
	else return ans;
}

// 
// inverse
//

int invf( float *inv, float *m, int n ){
	float max, ratio;
	int i, j, k, row;

	// initialize
	memset( inv, 0x00, n * n * sizeof( float ) );
	for( i = 0; i < n; i++ ) inv[ i * n + i ] = 1.0f;

	for( j = 0; j < n; j++ ){
		max = 0.0f;
		// find the pivot value
		for( i = j; i < n; i++ ){
			if( fabs( m[i * n + j] ) > fabs( max ) ) {
				max = m[i * n + j];
				row = i;
			}
		}
		// check if equal to 0
		if( max == 0.0f ) return 0;
		// row exchange if pivot isn't a diogonal element
		if( row != j ){
			for( k = 0; k < n; k++ ) {
				// swap
				max = m[row * n + k];
				m[row * n + k] = m[j * n + k];
				m[j * n + k] = max;
				max = inv[row * n + k];
				inv[row * n + k] = inv[j * n + k];
				inv[j * n + k] = max;
			}
		}
		// row operations
		for( i = 0; i < n; i++ ){
			if( i == j ){
				ratio = 1.0f / m[i * n + j];
				for( k = 0; k < n; k++ ){
					m[i * n + k] *= ratio;
					inv[i * n + k] *= ratio;
				}
			}
			else{
				if( m[i * n + j] == 0.0f ) continue;
				ratio = m[i * n + j] / m[j * n + j];
				for( k = 0; k < n; k++ ){
					m[i * n + k] -= ratio * m[j * n + k];
					inv[i * n + k] -= ratio * inv[j * n + k];
				}
			}
		}
	}
	return 1;
}

int invd( double *inv, double *m, int n ){
	double max, ratio;
	int i, j, k, row;

	// initialize
	memset( inv, 0x00, n * n * sizeof( double ) );
	for( i = 0; i < n; i++ ) inv[ i * n + i ] = 1.0;

	for( j = 0; j < n; j++ ){
		max = 0.0;
		// find the pivot value
		for( i = j; i < n; i++ ){
			if( fabs( m[i * n + j] ) > fabs( max ) ) {
				max = m[i * n + j];
				row = i;
			}
		}
		// check if equal to 0
		if( max == 0.0 ) return 0;
		// row exchange if pivot isn't a diogonal element
		if( row != j ){
			for( k = 0; k < n; k++ ) {
				// swap
				max = m[row * n + k];
				m[row * n + k] = m[j * n + k];
				m[j * n + k] = max;
				max = inv[row * n + k];
				inv[row * n + k] = inv[j * n + k];
				inv[j * n + k] = max;
			}
		}
		// row operations
		for( i = 0; i < n; i++ ){
			if( i == j ){
				ratio = 1.0 / m[i * n + j];
				for( k = 0; k < n; k++ ){
					m[i * n + k] *= ratio;
					inv[i * n + k] *= ratio;
				}
			}
			else{
				if( m[i * n + j] == 0.0 ) continue;
				ratio = m[i * n + j] / m[j * n + j];
				for( k = 0; k < n; k++ ){
					m[i * n + k] -= ratio * m[j * n + k];
					inv[i * n + k] -= ratio * inv[j * n + k];
				}
			}
		}
	}
	return 1;
}

//
// Psudo-Inverse, by SVD
//

int pinvf( Matf_out pinv, Matf_in x, Row row, Col col, Errf e ){
	Matf u, vT, buf, pu, pv, px;
	Vecf w;
	void *mem;
	float temp;
	int n, i, j, k;
	
	// memory pool
	n = ( row < col ) ? row : col;			// n = min( row, col )
	mem = calloc( row * col + n + 2 * SQR( n ), sizeof( float ) );
	if( mem == NULL ) return 0;

	memset( pinv, 0x00, row * col * sizeof( float ) );
	if( row <= col ){
		u = ( float * )mem;
		w = u + SQR( row );
		vT = w + n;
		buf = vT + row * col;
	}
	else{						// n = col
		// allocate memory
		u = ( float * )mem;
		w = u + row * col;
		vT = w + n;
		buf = vT + SQR( col );
	}

	// SVD
	svdf( x, u, w, vT, row, col, e, buf );

	// ( u * w-1 )T = w-1 * uT 
	pu = u;
	for( i = 0; i < row; i++ ){
		for( j = 0; j < n; j++, pu++ ){
			temp = ( w[j] < e ) ? 0.f : ( 1.f / w[j] );	// inverse value of w[i]
			*pu *= temp;
		}
	}
	// v * uT = ( u * vT )T
	pu = u;
	for( i = 0; i < row; i++ ){
		px = pinv + i;
		for( j = 0; j < col; j++ ){
			pv = vT + j;
			for( k = 0; k < n; k++, pv += col ){
				px[0] += pu[k] * pv[0];
			}
			px += row;
		}
		pu += n;
	}

	// free memory pool
	free( mem );
	return 1;
}

//
// Singular Value Decompisition
//

void svdf( Matf x, Matf u, Vecf w, Matf v, int row, int col, Errf e, Matf xx ){
	float *xi, *xj, *px, *pv, *pu, temp;
	int n, i, j, k;

	// initialize
	n = ( row < col ) ? row : col;			// n = min( row, col )
	memset( xx, 0x00, SQR( n ) * sizeof( float ) );
	
	if( row <= col ){
		// XXT
		for( i = 0; i < row; i++ ){
			px = pv = xx + i * n + i;
			for( j = i; j < row; j++, px++, pv += n ){
				xi = x + i * col;
				xj = x + j * col;
				for( k = 0; k < col; k++, xi++, xj++ ) px[0] += xi[0] * xj[0];
				pv[0] = px[0];
			}
		}
		Jacobif( xx, w, u, n, 30 );
		for( i = 0; i < n; i++ ) w[i] = ( w[i] > e ) ? ( float )sqrt( w[i] ) : 0.0f;
		// sorting
		for( i = 0; i < n - 1; i++ ){
			// pick max wi
			for( j = i + 1, k = i; j < n; j++ )	if( w[j] > w[k] ) k = j;
			if( i == k ) continue;			// don't need to swap
			// swap w
			temp = w[i];
			w[i] = w[k];
			w[k] = temp;
			// swap u
			pv = u + i;
			pu = u + k;
			for( j = 0; j < n; j++, pv += n, pu += n ) {
				temp = pu[0];
				pu[0] = pv[0];
				pv[0] = temp;
			}
		}
		// vT = uT * x
		pv = v;
		for( i = 0; i < row; i++ ){				// XT * ui
			for( j = 0; j < col; j++, pv++ ){	// u1 ... un
				pu = u + i;
				px = x + j;
				pv[0] = 0.0f;
				for( k = 0; k < row; k++, pu += n, px += col ) pv[0] += pu[0] * px[0];
				if( ( float )fabs( w[i] ) > e )	pv[0] /= w[i];			// normalize by eigenvalue
			}
		}
	}
	else{		// col < row
		// XTX
		for( i = 0; i < col; i++ ){
			px = pv = xx + i * n + i;
			for( j = i; j < col; j++, px++, pv += n ){
				xi = x + i;
				xj = x + j;
				for( k = 0; k < row; k++, xi += n, xj += n ) px[0] += xi[0] * xj[0];
				pv[0] = px[0];
			}
		}
		Jacobif( xx, w, v, n, 30 );
		for( i = 0; i < n; i++ ) w[i] = ( w[i] > e ) ? ( float )sqrt( w[i] ) : 0.0f;
		// sorting
		for( i = 0; i < n - 1; i++ ){
			// pick max wi
			for( j = i + 1, k = i; j < n; j++ )	if( w[j] > w[k] ) k = j;
			if( i == k ) continue;			// don't need to swap
			// swap w
			temp = w[i];
			w[i] = w[k];
			w[k] = temp;
			// swap v
			pv = v + i;
			pu = v + k;
			for( j = 0; j < n; j++, pv += n, pu += n ) {
				temp = pu[0];
				pu[0] = pv[0];
				pv[0] = temp;
			}
		}/**/
		// transpose v to vT
		for( i = 0; i < n; i++ ){
			pu = v + i * n + i + 1;
			pv = v + ( i + 1 ) * n + i;
			for( j = i + 1; j < n; j++, pu++, pv += n ){
				temp = pu[0];
				pu[0] = pv[0];
				pv[0] = temp;
			}
		}/**/
		// u = XvT
		pu = u;
		for( i = 0; i < row; i++ ){				// X * vi
			for( j = 0; j < col; j++, pu++ ){	// v1 ... vn
				pv = v + j * col;
				px = x + i * col;
				pu[0] = 0.0f;
				for( k = 0; k < col; k++, pv++, px++ ) pu[0] += pv[0] * px[0];
				if( ( float )fabs( w[j] ) > e )	pu[0] /= w[j];			// normalize by eigenvalue
			}
		}
	}
}

void svdd( Matd x, Matd u, Vecd w, Matd v, int row, int col, Errd e, Matd xx ){
	double *xi, *xj, *px, *pv, *pu, temp;
	int n, i, j, k;

	// initialize
	n = ( row < col ) ? row : col;			// n = min( row, col )
	memset( xx, 0x00, SQR( n ) * sizeof( double ) );
	
	if( row <= col ){
		// XXT
		for( i = 0; i < row; i++ ){
			px = pv = xx + i * n + i;
			for( j = i; j < row; j++, px++, pv += n ){
				xi = x + i * col;
				xj = x + j * col;
				for( k = 0; k < col; k++, xi++, xj++ ) px[0] += xi[0] * xj[0];
				pv[0] = px[0];
			}
		}
		Jacobid( xx, w, u, n, 30 );
		for( i = 0; i < n; i++ ) w[i] = ( w[i] > e ) ? sqrt( w[i] ) : 0.0;
		// sorting
		for( i = 0; i < n - 1; i++ ){
			// pick max wi
			for( j = i + 1, k = i; j < n; j++ )	if( w[j] > w[k] ) k = j;
			if( i == k ) continue;			// don't need to swap
			// swap w
			temp = w[i];
			w[i] = w[k];
			w[k] = temp;
			// swap u
			pv = u + i;
			pu = u + k;
			for( j = 0; j < n; j++, pv += n, pu += n ) {
				temp = pu[0];
				pu[0] = pv[0];
				pv[0] = temp;
			}
		}
		// vT = uT * x
		pv = v;
		for( i = 0; i < row; i++ ){				// XT * ui
			for( j = 0; j < col; j++, pv++ ){	// u1 ... un
				pu = u + i;
				px = x + j;
				pv[0] = 0.0f;
				for( k = 0; k < row; k++, pu += n, px += col ) {
					pv[0] += pu[0] * px[0];
				}
				if( fabs( w[i] ) > e )	pv[0] /= w[i];			// normalize by eigenvalue
			}
		}
	}
	else{		// col < row
		// XTX
		for( i = 0; i < col; i++ ){
			px = pv = xx + i * n + i;
			for( j = i; j < col; j++, px++, pv += n ){
				xi = x + i;
				xj = x + j;
				for( k = 0; k < row; k++, xi += n, xj += n ) px[0] += xi[0] * xj[0];
				pv[0] = px[0];
			}
		}
		Jacobid( xx, w, v, n, 30 );
		for( i = 0; i < n; i++ ) w[i] = ( w[i] > e ) ? sqrt( w[i] ) : 0.0;
		// sorting
		for( i = 0; i < n - 1; i++ ){
			// pick max wi
			for( j = i + 1, k = i; j < n; j++ )	if( w[j] > w[k] ) k = j;
			if( i == k ) continue;			// don't need to swap
			// swap w
			temp = w[i];
			w[i] = w[k];
			w[k] = temp;
			// swap v
			pv = v + i;
			pu = v + k;
			for( j = 0; j < n; j++, pv += n, pu += n ) {
				temp = pu[0];
				pu[0] = pv[0];
				pv[0] = temp;
			}
		}/**/
		// transpose v to vT
		for( i = 0; i < n; i++ ){
			pu = v + i * n + i + 1;
			pv = v + ( i + 1 ) * n + i;
			for( j = i + 1; j < n; j++, pu++, pv += n ){
				temp = pu[0];
				pu[0] = pv[0];
				pv[0] = temp;
			}
		}/**/
		// u = XvT
		pu = u;
		for( i = 0; i < row; i++ ){				// X * vi
			for( j = 0; j < col; j++, pu++ ){	// v1 ... vn
				pv = v + j * col;
				px = x + i * col;
				pu[0] = 0.0f;
				for( k = 0; k < col; k++, pv++, px++ ) pu[0] += pv[0] * px[0];
				if( fabs( w[j] ) > e )	pu[0] /= w[j];			// normalize by eigenvalue
			}
		}
	}
}


//
// Jacobi Algorithm
//

void Jacobif( float *A, float *eigen, float *v, int n, int time ){
	int i, j, p, q;
	float prev = 0.0f, sum, tangent, cosine, sine, tau, ptemp, qtemp,
		*pp, *pq, *qq, *a, *b;

	if(n<=1) return;

	//initialize A, V, eigen
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			if(i==j) v[i * n + j] = 1.0f;
			else v[i * n + j] = 0.0f;
		}
		eigen[i] = A[i * n + i];
	}
	//jacobi
	for( i = 0; i < time; i++, prev = sum ){
		sum = 0.0;
		for( p = 0; p < n-1; p++ ){
			a = A + p * n + p + 1;
			for( q = p + 1; q < n; q++, a++ ){
				sum += ( float )fabs( *a );
			}
		}
		if( sum == 0.0 || ( sum >= prev && i > 1 ) ) break;
		for( p = 0; p < n - 1; p++ ){
			pp = A + p * n + p;
			for( q = p + 1; q < n; q++ ){
				qq = A + q * n + q;
				pq = A + p * n + q;
				if( *pq == 0.0f && *qq - *pp == 0.0f ) tangent = ( float )tan( atan( 2.0f ) / 2.0f );
				else if( *qq - *pp == 0.0f ) tangent = 1.0f;
				else tangent = ( float )tan( atan( 2.0f * ( *pq ) / ( *qq - *pp ) ) / 2.0f );
				cosine = 1.0f / ( float )sqrt( tangent * tangent + 1.0f );
				sine = tangent * cosine;
				tau = sine / ( 1.0f + cosine );
				
				a = A + p;	b = A + q;
				for( j = 0; j < p; j++, a += n, b += n ){
					ptemp = *a; qtemp = *b;
					*a = ptemp - sine * ( qtemp + tau * ptemp );
					*b = qtemp + sine * ( ptemp - tau * qtemp );
				}
				*pp = *pp - tangent * ( *pq );
				a = A + p * n + p + 1;	b = A + ( p + 1 ) * n + q;
				for( j = p + 1; j < q; j++, a++, b += n ){
					ptemp = *a; qtemp = *b;
					*a = ptemp - sine * ( qtemp + tau * ptemp );
					*b = qtemp + sine * ( ptemp - tau * qtemp );
				}
				*qq = *qq + tangent * ( *pq );
				a = A + p * n + q + 1;	b = A + q * n + q + 1;
				for( j = q + 1; j < n; j++, a++, b++ ){
					ptemp = *a; qtemp = *b;
					*a = ptemp - sine * ( qtemp + tau * ptemp );
					*b = qtemp + sine * ( ptemp - tau * qtemp );
				}
				*pq = 0.0;
				for(j=0; j<n; j++){
					ptemp=v[j * n + p]; qtemp=v[j * n + q];
					v[j * n + p]=cosine*ptemp-sine*qtemp;
					v[j * n + q]=sine*ptemp+cosine*qtemp;
				}
			}
			a = A;
			for( j = 0; j < n; j++, a += n + 1 ) eigen[j] = *a;	// ( j, j )
		}
	}
}

void Jacobid( double *A, double *eigen, double *v, int n, int time ){
	int i, j, p, q;
	double prev = 0.0, sum, tangent, cosine, sine, tau, ptemp, qtemp,
		*pp, *pq, *qq, *a, *b;

	if(n<=1) return;

	//initialize A, V, eigen
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			if(i==j) v[i * n + j] = 1.0;
			else v[i * n + j] = 0.0;
		}
		eigen[i] = A[i * n + i];
	}
	//jacobi
	for( i = 0; i < time; i++, prev = sum ){
		sum = 0.0;
		for( p = 0; p < n-1; p++ ){
			a = A + p * n + p + 1;
			for( q = p + 1; q < n; q++, a++ ){
				sum += fabs( *a );
			}
		}
		if( sum == 0.0 || ( sum >= prev && i > 1 ) ) break;
		for( p = 0; p < n - 1; p++ ){
			pp = A + p * n + p;
			for( q = p + 1; q < n; q++ ){
				qq = A + q * n + q;
				pq = A + p * n + q;
				if( *pq == 0.0f && *qq - *pp == 0.0 ) tangent = tan( atan( 2.0 ) / 2.0 );
				else if( *qq - *pp == 0.0 ) tangent = 1.0;
				else tangent = tan( atan( 2.0 * ( *pq ) / ( *qq - *pp ) ) / 2.0 );
				cosine = 1.0 / sqrt( tangent * tangent + 1.0 );
				sine = tangent * cosine;
				tau = sine / ( 1.0 + cosine );
				
				a = A + p;	b = A + q;
				for( j = 0; j < p; j++, a += n, b += n ){
					ptemp = *a; qtemp = *b;
					*a = ptemp - sine * ( qtemp + tau * ptemp );
					*b = qtemp + sine * ( ptemp - tau * qtemp );
				}
				*pp = *pp - tangent * ( *pq );
				a = A + p * n + p + 1;	b = A + ( p + 1 ) * n + q;
				for( j = p + 1; j < q; j++, a++, b += n ){
					ptemp = *a; qtemp = *b;
					*a = ptemp - sine * ( qtemp + tau * ptemp );
					*b = qtemp + sine * ( ptemp - tau * qtemp );
				}
				*qq = *qq + tangent * ( *pq );
				a = A + p * n + q + 1;	b = A + q * n + q + 1;
				for( j = q + 1; j < n; j++, a++, b++ ){
					ptemp = *a; qtemp = *b;
					*a = ptemp - sine * ( qtemp + tau * ptemp );
					*b = qtemp + sine * ( ptemp - tau * qtemp );
				}
				*pq = 0.0;
				for(j=0; j<n; j++){
					ptemp=v[j * n + p]; qtemp=v[j * n + q];
					v[j * n + p]=cosine*ptemp-sine*qtemp;
					v[j * n + q]=sine*ptemp+cosine*qtemp;
				}
			}
			a = A;
			for( j = 0; j < n; j++, a += n + 1 ) eigen[j] = *a;	// ( j, j )
		}
	}
}

//
// Power Iteration
//

float PowerIterf( Matf A, Vecf u, int n, int time, Errf e, Vecf v ){
	float *a, max = 0.0f, prev = FLT_MAX;
	int i, j, k;

	// initial guess
	for( i = 0; i < n; i++ ) u[i] = ( float )rand() / ( float )RAND_MAX;

	for( k = 0; k < time; k++ ){
		if( ( float )fabs( max - prev ) < e ) break;
		// calculate v
		a = A;
		memset( v, 0x00, n * sizeof( float ) );
		for( i = 0; i < n; i++ )
			for( j = 0; j < n; j++, a++ ) v[i] += u[j] * a[0];
		// calculate u
		prev = max;
		max = 0.0f;
		for( i = 0; i < n; i++ ) {
			if( ( float )fabs( v[i] ) > max ) max = ( float )fabs( v[i] );
		}
		if( max > e ) for( i = 0; i < n; i++ ) u[i] = v[i] / max;
	}
	return max;
}

double PowerIterd( Matd A, Vecd u, int n, int time, Errd e, Vecd v ){
	double *a, max = 0.0, prev = DBL_MAX;
	int i, j, k;

	// initial guess
	for( i = 0; i < n; i++ ) u[i] = ( double )rand() / ( double )RAND_MAX;

	for( k = 0; k < time; k++ ){
		if( fabs( max - prev ) < e ) break;
		// calculate v
		a = A;
		memset( v, 0x00, n * sizeof( double ) );
		for( i = 0; i < n; i++ )
			for( j = 0; j < n; j++, a++ ) v[i] += u[j] * a[0];
		// calculate u
		prev = max;
		max = 0.0f;
		for( i = 0; i < n; i++ ) {
			if( fabs( v[i] ) > max ) max = fabs( v[i] );
		}
		if( max > e ) for( i = 0; i < n; i++ ) u[i] = v[i] / max;
	}
	return max;
}

//
// Inverse Iteration
//

float InvIterf( float *A, float *b0, float shift, int n, int time, float err, float *b1, int *ind ){
	float sum, g, prev;
	int i, k;

	// initialize
	for( i = 0; i < n; i++ ) A[ i * n + i ] -= shift;
	LUfactf( A, ind, n, b1 );

	// keep guess when grown factor too small
	do{
		for( i = 0; i < n; i++ ) b0[i] = b1[i] = ( float )rand() / ( float )RAND_MAX;
		// normalize for both
		sum = 0.0f;
		for( i = 0; i < n; i++ ) sum += SQR( b0[i] );
		sum = ( float )sqrt( sum );
		for( i = 0; i < n; i++ ) b0[i] = b1[i] = b0[i] / sum;
		// solve ( A- aI )b1 = b0
		SolveLUf( A, ind, b1, n );
		// normalized b1
		g = 0.0f;					// grown factor: | b0 |
		for( i = 0; i < n; i++ ) g += SQR( b1[i] );
		g = ( float )sqrt( g );
	} while( g < err );	
	
	// normalize
	prev = 0.0f;							// | b1 - b0 |
	for( i = 0; i < n; i++ ) {	
		b1[i] /= g;						// nomralize
		prev += SQR( b0[i] - b1[i] );	
		b0[i] = b1[i];					// updata b0
	}

	for( k = 0; k < time; k++ ){
		// solve ( A- aI )y = b
		SolveLUf( A, ind, b1, n );
		// normalized
		sum = 0.0f;
		for( i = 0; i < n; i++ ) sum += SQR( b1[i] );
		sum = ( float )sqrt( sum );
		g = 0.0f;					// | b1 - b0 |
		for( i = 0; i < n; i++ ) {
			b1[i] /= sum;			// normalize 
			g += SQR( b0[i] - b1[i] );	// | b1 - b0 |
		}
		// break contition
		if( g < err ) break;
		// updata shift value when shift of | b1 - b0 | < error bound
		if( prev - g < err ){
			prev = g;
			g = 0.0f;
			for( i = 0; i < n; i++ ) g += b0[i] * b1[i];
			g = 1.0f / ( g * sum );	// unnormalized y 
			LUcombf( A, ind, n );
			for( i = 0; i < n; i++ ) A[ i * n + i ] += shift;
			shift += g;
			for( i = 0; i < n; i++ ) A[ i * n + i ] -= shift;
			LUfactf( A, ind, n, b0 );
		}
		else prev = g;
		// updata b0
		for( i = 0; i < n; i++ ) b0[i] = b1[i];
	}
	return shift;
}

double InvIterd( double *A, double *b0, double shift, int n, int time, double err, double *b1, int *ind ){
	double sum, g, prev;
	int i, k;

	// initialize
	for( i = 0; i < n; i++ ) A[ i * n + i ] -= shift;
	LUfactd( A, ind, n, b1 );

	// keep guess when grown factor too small
	do{
		for( i = 0; i < n; i++ ) b0[i] = b1[i] = ( double )rand() / ( double )RAND_MAX;
		// normalize for both
		sum = 0.0;
		for( i = 0; i < n; i++ ) sum += SQR( b0[i] );
		sum = sqrt( sum );
		for( i = 0; i < n; i++ ) b0[i] = b1[i] = b0[i] / sum;
		// solve ( A- aI )b1 = b0
		SolveLUd( A, ind, b1, n );
		// normalized b1
		g = 0.0;					// grown factor: | b0 |
		for( i = 0; i < n; i++ ) g += SQR( b1[i] );
		g = sqrt( g );
	} while( g < err );	
	
	// normalize
	prev = 0.0f;							// | b1 - b0 |
	for( i = 0; i < n; i++ ) {	
		b1[i] /= g;						// nomralize
		prev += SQR( b0[i] - b1[i] );	
		b0[i] = b1[i];					// updata b0
	}

	for( k = 0; k < time; k++ ){
		// solve ( A- aI )y = b
		SolveLUd( A, ind, b1, n );
		// normalized
		sum = 0.0;
		for( i = 0; i < n; i++ ) sum += SQR( b1[i] );
		sum = sqrt( sum );
		g = 0.0;					// | b1 - b0 |
		for( i = 0; i < n; i++ ) {
			b1[i] /= sum;			// normalize 
			g += SQR( b0[i] - b1[i] );	// | b1 - b0 |
		}
		// break contition
		if( g < err ) break;
		// updata shift value when shift of | b1 - b0 | < error bound
		if( prev - g < err ){
			prev = g;
			g = 0.0;
			for( i = 0; i < n; i++ ) g += b0[i] * b1[i];
			g = 1.0 / ( g * sum );	// unnormalized y 
			LUcombd( A, ind, n );
			for( i = 0; i < n; i++ ) A[ i * n + i ] += shift;
			shift += g;
			for( i = 0; i < n; i++ ) A[ i * n + i ] -= shift;
			LUfactd( A, ind, n, b0 );
		}
		else prev = g;
		// updata b0
		for( i = 0; i < n; i++ ) b0[i] = b1[i];
	}
	return shift;
}

//
// LU Factorization
//

int LUfactf( float *A, int *index, int n, float *v ){
	float sum, big, *a, *b;
	int i, j, k, row;

	// calculate largest value per row
	a = A;
	for( i = 0; i < n; i++ ){
		big = 0.0;
		for( j = 0; j < n; j++, a++ ){
			if( ( sum = ( float )fabs( *a ) ) > big ) big = sum;
		}
		if( big == 0.0f ) return 0;
		v[i] = 1.0f / big;
	}
	for( j = 0; j < n; j++ ){
		// calculate U when i < j
		for( i = 0; i < j; i++ ){
			sum = A[ i * n + j ];
			for( k = 0; k < i; k++ ) sum -= A[ i * n + k ] * A[ k * n + j ];
			A[ i * n + j ] = sum;
		}
		big = 0.0f;
		// calculate L wneh i >= j;
		for( i = j; i < n; i++ ){
			sum = A[ i * n + j ];
			for( k = 0; k < j; k++ ) sum -= A[ i * n + k ] * A[ k * n + j ];
			A[ i * n + j ] = sum;
			sum = v[i] * ( float )fabs( sum ) ;
			if( sum >= big ){
				big = sum;	row = i;
			}
		}
		if( j != row ){		// need to change pivot
			a = A + row * n;
			b = A + j * n;
			for( k = 0; k < n; k++, a++, b++ ){	// row interchange
				sum = *a;
				*a = *b;
				*b = sum;
			}
			v[ row ] = v[j];
		}
		index[j] = row;				// row index
		a = A + j * n + j;
		if( *a == 0 ) *a = FLT_MIN;	// sigular case
		if( j < n - 1 ){			// devided by pivot
			sum = 1.0f / *a;
			a += n;
			for( i = j + 1; i < n; i++, a += n ) *a *= sum;
		}
	}
	return 1;
}

int LUfactd( double *A, int *index, int n, double *v ){
	double sum, big, *a, *b;
	int i, j, k, row;

	// calculate largest value per row
	a = A;
	for( i = 0; i < n; i++ ){
		big = 0.0;
		for( j = 0; j < n; j++, a++ ){
			if( ( sum = fabs( *a ) ) > big ) big = sum;
		}
		if( big == 0.0 ) return 0;
		v[i] = 1.0 / big;
	}
	for( j = 0; j < n; j++ ){
		// calculate U when i < j
		for( i = 0; i < j; i++ ){
			sum = A[ i * n + j ];
			for( k = 0; k < i; k++ ) sum -= A[ i * n + k ] * A[ k * n + j ];
			A[ i * n + j ] = sum;
		}
		big = 0.0f;
		// calculate L wneh i >= j;
		for( i = j; i < n; i++ ){
			sum = A[ i * n + j ];
			for( k = 0; k < j; k++ ) sum -= A[ i * n + k ] * A[ k * n + j ];
			A[ i * n + j ] = sum;
			sum = v[i] * fabs( sum ) ;
			if( sum >= big ){
				big = sum;	row = i;
			}
		}
		if( j != row ){		// need to change pivot
			a = A + row * n;
			b = A + j * n;
			for( k = 0; k < n; k++, a++, b++ ){	// row interchange
				sum = *a;
				*a = *b;
				*b = sum;
			}
			v[ row ] = v[j];
		}
		index[j] = row;				// row index
		a = A + j * n + j;
		if( *a == 0 ) *a = FLT_MIN;	// sigular case
		if( j < n - 1 ){			// devided by pivot
			sum = 1.0 / *a;
			a += n;
			for( i = j + 1; i < n; i++, a += n ) *a *= sum;
		}
	}
	return 1;
}

void SolveLUf( float *A, int *ind, float *b, int n ){
	float sum;
	int i, j, ip;

	for( i = 0; i < n; i++ ){		// solve Ly = b part
		ip = ind[i];
		sum = b[ip];
		b[ip] = b[i];
		if( i ){
			for( j = 0; j < i; j++ ) sum -= A[ i * n + j ] * b[j];
		}
		b[i] = sum;
	}
	for( i = n - 1; i >= 0; i-- ){	// Solve Ux = y part
		sum = b[i];
		for( j = i + 1; j < n; j++ ) sum -= A[ i * n + j ] * b[j];
		b[i] = sum / A[ i * n + i ];
	}
}

void SolveLUd( double *A, int *ind, double *b, int n ){
	double sum;
	int i, j, ip;

	for( i = 0; i < n; i++ ){		// solve Ly = b part
		ip = ind[i];
		sum = b[ip];
		b[ip] = b[i];
		if( i ){
			for( j = 0; j < i; j++ ) sum -= A[ i * n + j ] * b[j];
		}
		b[i] = sum;
	}
	for( i = n - 1; i >= 0; i-- ){	// Solve Ux = y part
		sum = b[i];
		for( j = i + 1; j < n; j++ ) sum -= A[ i * n + j ] * b[j];
		b[i] = sum / A[ i * n + i ];
	}
}

void LUcombf( float *A, int *ind, int n ){
	float sum, *a, *b;
	int i, j, k;

	for( j = n - 1; j >= 0; j-- ){
		for( i = n - 1; i > j; i-- ){		// case i >= j 
			sum = 0.0f;
			for( k = 0; k <= j; k++ ) sum += A[ i * n + k ] * A[ k * n + j ];
			A[ i * n + j ] = sum;
		}
		for( i = j; i >= 0; i-- ){			// acse i < j
			sum = 0.0;
			for( k = 0; k < i; k++ ) sum += A[ i * n + k ] * A[ k * n + j ];
			A[ i * n + j ] += sum;
		}
	}
	for( i = 0; i < n; i++ ){			// row interchange
		j = ind[i];
		if( j == i ) continue;
		// change!
		a = A + i * n;
		b = A + j * n;
		for( k = 0; k < n; k++, a++, b++ ){
			sum = *a;
			*a = *b;
			*b = sum;
		}
	}
}

void LUcombd( double *A, int *ind, int n ){
	double sum, *a, *b;
	int i, j, k;

	for( j = n - 1; j >= 0; j-- ){
		for( i = n - 1; i > j; i-- ){		// case i >= j 
			sum = 0.0;
			for( k = 0; k <= j; k++ ) sum += A[ i * n + k ] * A[ k * n + j ];
			A[ i * n + j ] = sum;
		}
		for( i = j; i >= 0; i-- ){			// acse i < j
			sum = 0.0;
			for( k = 0; k < i; k++ ) sum += A[ i * n + k ] * A[ k * n + j ];
			A[ i * n + j ] += sum;
		}
	}
	for( i = 0; i < n; i++ ){			// row interchange
		j = ind[i];
		if( j == i ) continue;
		// change!
		a = A + i * n;
		b = A + j * n;
		for( k = 0; k < n; k++, a++, b++ ){
			sum = *a;
			*a = *b;
			*b = sum;
		}
	}
}