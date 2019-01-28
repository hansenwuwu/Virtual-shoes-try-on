
// feature selection and generation

#include "pr.h"
#include <numerical\linear.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <stdio.h>

// 
// macro
//

#define SQR(x)		( (x) * (x) )
#define min( a, b )	( ( (a) < (b) ) ? (a) : (b) )
#define BIG			10000000000

//
// functions declaration
//

void kNNMean( double *u, double *x, int xi, int n, int d, int k, void *mem );
void kNNMean( double *u, double *x, double *xi, int n, int d, int k, void *mem );

//
// Feature Selection
//

// Scatter Matrix

void ScatterB( double *s, double *x, int m, int n, int d ){
	double *px, *u, *pu, *u0, *ps1, *ps2, p;
	int i, j, k;
	void *mem;

	// memory pool
	mem = malloc( d * ( m + 1 ) * sizeof( double ) );
	u0 = ( double* )mem;	// d
	u = u0 + d;				// d * m

	// initialize
	for( i = 0, j = d * d; i < j; i++ ) s[i] = 0.0;
	for( i = 0, j = d * m + d; i < j; i++ ) u0[i] = 0.0;

	// compute mean each class
	px = x;
	pu = u;
	p = 1.0 / ( double )n;
	for( i = 0; i < m; i++, pu += d ){
		for( j = 0; j < n; j++ ){
			for( k = 0; k < d; k++, px++ ) pu[k] += px[0];
		}
		for( j = 0; j < d; j++ ) pu[j] *= p;	// normalize
	}
	// total mean
	pu = u;
	p = 1.0 / ( double )m;
	for( i = 0; i < m; i++, pu += d ){
		for( j = 0; j < d; j++ ) u0[j] += pu[j];
	}
	for( i = 0; i < d; i++ ) u0[i] *= p;

	// Sb = E[ ( ui - u0 )( ui - u0 )T ]
	pu = u;
	for( k = 0; k < m; k++, pu += d ){
		for( i = 0; i < d; i++ ){
			ps1 = ps2 = s + i * d + i;
			for( j = i; j < d; j++, ps1++, ps2 += d ){
				ps1[0] += ( pu[i] - u0[i] ) * ( pu[j] - u0[j] );
				ps2[0] = ps1[0];
			}
		}
	}
	// normalize
	// p = 1.0 / ( double )m;
	for( i = 0, j = d * d; i < j; i++ ) s[i] *= p;
	// free memory pool
	free( mem );
}

void ScatterW( double *s, double *x, int m, int n, int d ){
	double *px, *u, *pu, *ps1, *ps2, p;
	int i, j, k, t;
	void *mem;

	// memory pool
	mem = malloc( d * m * sizeof( double ) );
	u = ( double* )mem;							// d * m

	// initialize
	for( i = 0, j = d * d; i < j; i++ ) s[i] = 0.0;
	for( i = 0, j = d * m; i < j; i++ ) u[i] = 0.0;

	// compute mean each class
	px = x;
	pu = u;
	p = 1.0 / ( double )n;
	for( i = 0; i < m; i++, pu += d ){
		for( j = 0; j < n; j++ ){
			for( k = 0; k < d; k++, px++ ) pu[k] += px[0];
		}
		for( j = 0; j < d; j++ ) pu[j] *= p;	// normalize
	}
	// Sw = E[ ( x - ui )( x - ui )T ]
	px = x;
	pu = u;
	t = m * n;
	for( k = 0; k < t; k++, px += d ){		// for every data
		if( k % n == 0 && k > 0 ) pu += d;	// mean of class
		for( i = 0; i < d; i++ ){
			ps1 = ps2 = s + i * d + i;
			for( j = i; j < d; j++, ps1++, ps2 += d ){
				ps1[0] += ( px[i] - pu[i] ) * ( px[j] - pu[j] );
				ps2[0] = ps1[0];
			}
		}
	}
	// normalize
	p = 1.0 / ( double )t;
	for( i = 0, j = d * d; i < j; i++ ) s[i] *= p;
	// free memory pool
	free( mem );
}

void ScatterM( double *s, double *x, int m, int n, int d ){
	double *px, *u0, *ps1, *ps2, p;
	int i, j, k;
	void *mem;

	// memory pool
	mem = malloc( d * sizeof( double ) );
	u0 = ( double* )mem;	// d

	// initialize
	for( i = 0, j = d * d; i < j; i++ ) s[i] = 0.0;
	for( i = 0, j = d; i < j; i++ ) u0[i] = 0.0;

	// compute mean each class
	px = x;
	p = 1.0 / ( double )( n * m );
	for( i = 0; i < m; i++ ){
		for( j = 0; j < n; j++ ){
			for( k = 0; k < d; k++, px++ ) u0[k] += px[0];
		}
	}
	for( j = 0; j < d; j++ ) u0[j] *= p;	// normalize

	// Sm = E[ ( x - u0 )( x - u0 )T ]
	px = x;
	for( k = m * n; k > 0; k--, px += d ){
		for( i = 0; i < d; i++ ){
			ps1 = ps2 = s + i * d + i;
			for( j = i; j < d; j++, ps1++, ps2 += d ){
				ps1[0] += ( px[i] - u0[i] ) * ( px[j] - u0[j] );
				ps2[0] = ps1[0];
			}
		}
	}
	// normalize
	// p = 1.0 / ( double )( n * m );
	for( i = 0, j = d * d; i < j; i++ ) s[i] *= p;
	// free memory pool
	free( mem );
}

void CompSbSm( double *sb, double *sm, int d, int m, va_list va ){
	double **x, *px, *u, *pu, *sij, *sji, *u0, p;
	int *n, i, j, k, l, total;
	void *mem;

	// memory pool
	mem = malloc( m * ( sizeof( double* ) + sizeof( int ) ) + ( m + 1 ) * d * sizeof( double ) );
	x = ( double** )mem;			// m
	n = ( int* )( x + m );			// m
	u0 = ( double* )( n + m );		// d
	u = u0 + d;						// m * d

	// initialize
	for( i = 0; i < m; i++ ) x[i] = va_arg( va, double * );
	for( i = 0; i < m; i++ ) n[i] = va_arg( va, int );	
	// zero
	for( i = 0, j = d * d; i < j; i++ ) sb[i] = sm[i] = 0.0;
	for( i = 0; i < d; i++ ) u0[i] = 0.0;
	for( i = 0, j = m * d; i < j; i++ ) u[i] = 0.0;

	// compute u and u0
	pu = u;
	for( i = total = 0; i < m; i++, pu += d ){	// for class i
		px = x[i];
		p = 1.0 / ( double )n[i];
		for( j = 0; j < n[i]; j++, px += d ){	// data in class i
			for( k = 0; k < d; k++ ){			// d dimension
				pu[k] += px[k];
				u0[k] += px[k];
			}
		}
		for( k = 0; k < d; k++ ) pu[k] *= p;
		total += j;
	}
	p = 1.0 / ( double )total;
	for( i = 0; i < d; i++ ) u0[i] *= p;

	// compute Sb = Sum Pk * ( uk - u0 )( uk - u0 )T
	pu = u;
	for( k = 0; k < m; k++, pu += d ){				// for each class k
		p = ( double )n[k] / ( double )total;
		for( i = 0; i < d; i++ ) {					// d * d scatter matrix Sb
			sij = sji = sb + i * d + i;
			for( j = i; j < d; j++, sij++, sji += d ){
				sij[0] += p * ( pu[i] - u0[i] ) * ( pu[j] - u0[j] );
				sji[0] = sij[0];
			}
		}
	}
	// compute Sm = E[ ( x - u0 )( x - u0 )T ]
	for( k = 0; k < m; k++ ){
		px = x[k];
		for( l = 0; l < n[k]; l++, px += d ){
			for( i = 0; i < d; i++ ) {				// d * d scatter matrix Sm
				sij = sm + i * d + i;
				for( j = i; j < d; j++, sij++ )
					sij[0] += ( px[i] - u0[i] ) * ( px[j] - u0[j] );
			}
		}
	}
	// normalize Sm
	p = 1.0 / ( double )total;
	for( i = 0; i < d; i++ ) {				// d * d scatter matrix Sm
		sij = sji = sm + i * d + i;
		for( j = i; j < d; j++, sij++, sji += d )
			sji[0] = sij[0] = sij[0] * p;
	}
	free( mem );
}

// Nonparametric Scatter Matrix

void kNNMean( double *u, double *x, int xi, int n, int d, int k, void *mem ){
	double *nb, *pn1, *pn2, *px, *pxi, *dist, p;
	int i, j, b, a;

	// memory pool
	dist = ( double* )mem;	// knn
	nb = dist + k;			// knn * d

	// initialize
	for( i = 0; i < k; i++ ) dist[i] = DBL_MAX;	
	pxi = x + xi * d;

	px = x;
	for( i = 0; i < n; i++, px += d ){				// n data in class
		if( xi == i ) continue;
		// compute distance
		p = 0.0;
		for( j = 0; j < d; j++ ) p += SQR( pxi[j] - px[j] );
		// update knn
		for( j = 0; j < k; j++ ){
			if( p < dist[j] ){
				// offset
				pn1 = nb + ( k - 2 ) * d;
				pn2 = pn1 + d;
				for( a = k - 1; a > j; a--, pn1 -= d, pn2 -= d ){
					for( b = 0; b < d; b++ ) pn2[b] = pn1[b];	// copy data to next
					dist[a] = dist[ a - 1 ];					// copy distance to next
				}
				// save data
				pn1 = nb + j * d;
				for( a = 0; a < d; a++ ) pn1[a] = px[a];		// data
				dist[j] = p;									// distance
				break;
			}
		}
	}
	// compute mean of knn
	p = 1.0 / ( double )k;
	for( i = 0; i < d; i++ ){
		u[i] = 0.0;
		pn1 = nb + i;
		for( j = 0; j < k; j++, pn1 += d ) u[i] += pn1[0];
		u[i] *= p;
	}
}

void kNNMean( double *u, double *x, double *xi, int n, int d, int k, void *mem ){
	double *nb, *pn1, *pn2, *px, *dist, p;
	int i, j, b, a;

	// memory pool
	dist = ( double* )mem;	// knn
	nb = dist + k;			// knn * d

	// initialize
	for( i = 0; i < k; i++ ) dist[i] = DBL_MAX;	

	px = x;
	for( i = 0; i < n; i++, px += d ){				// n data in class
		// compute distance
		p = 0.0;
		for( j = 0; j < d; j++ ) p += SQR( xi[j] - px[j] );
		// update knn
		for( j = 0; j < k; j++ ){
			if( p < dist[j] ){
				// offset
				pn1 = nb + ( k - 2 ) * d;
				pn2 = pn1 + d;
				for( a = k - 1; a > j; a--, pn1 -= d, pn2 -= d ){
					for( b = 0; b < d; b++ ) pn2[b] = pn1[b];	// copy data to next
					dist[a] = dist[ a - 1 ];					// copy distance to next
				}
				// save data
				pn1 = nb + j * d;
				for( a = 0; a < d; a++ ) pn1[a] = px[a];		// data
				dist[j] = p;									// distance
				break;
			}
		}
	}
	// compute mean of knn
	p = 1.0 / ( double )k;
	for( i = 0; i < d; i++ ){
		u[i] = 0.0;
		pn1 = nb + i;
		for( j = 0; j < k; j++, pn1 += d ) u[i] += pn1[0];
		u[i] *= p;
	}
}

void NPScatterB( double *s, double *x, int m, int n, int d, int knn, double param ){
	double *px, *u1, *u2, *ps1, *ps2, *nb, *dist, p, d1, d2;
	int i, j, k, ni, a, b;
	void *mem;

	// memory pool
	mem = malloc( ( d * ( 2 + knn ) + knn ) * sizeof( double ) );
	u1 = ( double* )mem;	// d
	u2 = u1 + d;			// d
	nb = u2 + d;			// d * knn
	dist = nb + d * knn;	// knn

	// initialize
	for( i = 0, j = d * d; i < j; i++ ) s[i] = 0.0;
	param *= 0.5;

	// Sb' = param * E[ ( x - xE )( x - xE )T ]
	px = x;
	for( i = 0; i < m; i++ ){
		// data in class i
		for( ni = 0; ni < n; ni++, px += d ){
			// kNN mean in class i
			kNNMean( u1, x + i * n * d, ni, n, d, knn, nb ); 
			// compute distance from xi to knn mean in class i
			d1 = 0.0;
			for( k = 0; k < d; k++ ) d1 += SQR( px[k] - u1[k] );
			d1 = pow( d1, param );
			// class j
			for( j = i + 1; j < m; j++ ){
				// kNN mean in class j
				kNNMean( u2, x + j * n * d, px, n, d, knn, nb );
				// compute distance from xi to knn mean in class j
				d2 = 0.0;
				for( k = 0; k < d; k++ ) d2 += SQR( px[k] - u2[k] );
				d2 = pow( d2, param );
				// weight
				p = min( d1, d2 ) / ( d1 + d2 );
				// Scatter matrix
				for( a = 0; a < d; a++ ){
					ps1 = ps2 = s + a * d + a;
					for( b = a; b < d; b++, ps1++, ps2 += d ){
						ps1[0] += p * ( px[a] - u2[a] ) * ( px[b] - u2[b] );
						ps2[0] = ps1[0];
					}
				}
			}
		}
	}
	// free memory pool
	free( mem );
}

// Feature Criteria

double J1( double *x, int m, int n, int d ){
	double *sm, *sw, ans, p;
	int i;
	void *mem;

	// memory pool
	mem = malloc( 2 * d * d * sizeof( double ) );
	sm = ( double* )mem;
	sw = sm + d * d;

	// compute J1
	ScatterW( sw, x, m, n, d );
	ScatterM( sm, x, m, n, d );
	ans = p = 0.0;
	for( i = 0; i < d * d; i += d + 1 ){
		ans += sm[i];
		p += sw[i];
	}
	free( mem );
	return ans / p;
}

double J1( int d, int m, ... ){
	double *sb, *sm, *psb, *psm, tsm, tsw;
	int i;
	va_list va;
	void *mem;

	// memory pool
	mem = malloc( 2 * d * d * sizeof( double ) );
	sb = ( double* )mem;			// d * d
	sm = sb + d * d;				// d * d

	// initialize
	va_start( va, m );
	CompSbSm( sb, sm, d, m, va );
	va_end( va );

	psb = sb;	psm = sm;
	tsm = tsw = 0.0;
	for( i = 0; i < d; i++, psb += d + 1, psm += d + 1 ){
		tsm += psm[0];
		tsw += psm[0] - psb[0];		// Sw = Sm - Sb
	}
	free( mem );	
	return tsm / tsw;
}

double J2( double *x, int m, int n, int d ){
	double *sm, *sw, *s, *psm, *psw, *ps, ans;
	int i, j, k;
	void *mem;

	// memory pool
	mem = malloc( 3 * d * d * sizeof( double ) );
	sm = ( double* )mem;
	sw = sm + d * d;
	s = sw + d * d;

	// compute J2 = | Sw-1 * Sm |
	ScatterW( sw, x, m, n, d );
	inverse( sw, d, s );
	ScatterM( sm, x, m, n, d );

	ps = s;
	for( i = 0; i < d; i++ ){
		for( j = 0; j < d; j++, ps++ ){
			psw = sw + i * d;
			psm = sm + j;
			ps[0] = 0.0;
			for( k = 0; k < d; k++, psw++, psm += d ) 
				ps[0] += psw[0] * psm[0];
		}
	}
	ans = det( s, d, sw );
	free( mem );
	return ans;
}

double J2( int d, int m, ... ){
	double *sw, *sm, *s, *psw, *psm, *ps, ans;
	int i, j, k;
	va_list va;
	void *mem;

	// memory pool
	mem = malloc( 3 * d * d * sizeof( double ) );
	sw = ( double* )mem;			// d * d
	sm = sw + d * d;				// d * d
	s = sm + d * d;					// d * d

	// initialize
	va_start( va, m );
	CompSbSm( sw, sm, d, m, va );	// now Sw := Sb
	va_end( va );

	// compute Sw
	psw = sw;	psm = sm;
	for( i = 0, j = d * d; i < j; i++, psw++, psm++ )
		psw[0] = psm[0] - psw[0];		// Sw = Sm - Sb
	// compute Sw^-1
	inverse( sw, d, s );

	// Sw^-1 * Sm
	ps = s;
	for( i = 0; i < d; i++ ){
		for( j = 0; j < d; j++, ps++ ){
			ps[0] = 0.0;
			psw = sw + i * d;
			psm = sm + j;
			for( k = 0; k < d; k++, psw++, psm += d ) ps[0] += psw[0] * psm[0];
		}
	}
	ans = det( s, d, sw );
	free( mem );
	return ans;
}

double J3( double *x, int m, int n, int d ){
	double *sb, *sw, *s, *psb, *psw, *ps, ans;
	int i, j, k;
	void *mem;

	// memory pool
	mem = malloc( 3 * d * d * sizeof( double ) );
	sb = ( double* )mem;
	sw = sb + d * d;
	s = sw + d * d;

	// compute J3 = trace( Sw-1 * Sb )
	ScatterW( sw, x, m, n, d );
	inverse( sw, d, s );
	ScatterB( sb, x, m, n, d );

	ps = s;
	for( i = 0; i < d; i++ ){
		for( j = 0; j < d; j++, ps++ ){
			psw = sw + i * d;
			psb = sb + j;
			ps[0] = 0.0;
			for( k = 0; k < d; k++, psw++, psb += d ) 
				ps[0] += psw[0] * psb[0];
		}
	}
	// trace
	ans = 0.0;
	ps = s;
	for( i = 0; i < d; i++, ps += d + 1 ) ans += ps[0];

	free( mem );
	return ans;
}

double J3( int d, int m, ... ){
	double *sw, *sb, *s, *psw, *psb, *ps, tr = 0.0;
	int i, j, k;
	va_list va;
	void *mem;

	// memory pool
	mem = malloc( 3 * d * d * sizeof( double ) );
	sw = ( double* )mem;			// d * d
	sb = sw + d * d;				// d * d
	s = sb + d * d;					// d * d

	// initialize
	va_start( va, m );
	CompSbSm( sb, sw, d, m, va );	// now Sw := Sm
	va_end( va );

	// compute Sw
	psw = sw;	psb = sb;
	for( i = 0, j = d * d; i < j; i++, psw++, psb++ )
		psw[0] -= psb[0];			// Sw = Sm - Sb
	// compute Sw^-1
	inverse( sw, d, s );

	// Sw^-1 * Sb
	ps = s;
	for( i = 0; i < d; i++ ){
		for( j = 0; j < d; j++, ps++ ){
			ps[0] = 0.0;
			psw = sw + i * d;
			psb = sb + j;
			for( k = 0; k < d; k++, psw++, psb += d ) ps[0] += psw[0] * psb[0];
		}
	}
	// trace
	ps = s;
	for( i = 0; i < d; i++, ps += d + 1 ) tr += ps[0];
	free( mem );
	return tr;
}

double Fdr( int d, int m, ... ){
	double *u, *v, **x, *px, *pu, *pv, *pu2, *pv2, p, fdr, temp;
	int *n, i, j, k;
	va_list va;
	void *mem;

	// memory pool
	mem = malloc( m * ( 2 * d * sizeof( double ) + sizeof( double* ) + sizeof( int ) ) );
	u = ( double* )mem;				// m * d
	v = u + m * d;					// m * d
	x = ( double** )( v + m * d );	// m
	n = ( int* )( x + m );			// m
	
	// initialize
	va_start( va, m );
	for( i = 0; i < m; i++ ) x[i] = va_arg( va, double* );
	for( i = 0; i < m; i++ ) n[i] = va_arg( va, int );
	va_end( va );
	for( i = 0, j = m * d; i < j; i++ ) u[i] = v[i] = 0.0;

	// compute E[x] and V[x]
	pu = u;	pv = v;
	for( i = 0; i < m; i++, pu += d, pv += d ){	// for each class
		px = x[i];	
		for( j = 0; j < n[i]; j++ ){			// for each data
			for( k = 0; k < d; k++, px++ ){		// for each dimension
				pu[k] += px[0];
				pv[k] += px[0] * px[0];
			}
		}
		p = 1.0 / ( double )j;
		for( j = 0; j < d; j++ ){
			pu[j] *= p;							// E[x]
			pv[j] = pv[j] * p - SQR( pu[j] );	// E[x2] - E(x)2
		}
	}

	// FDR = sum( i, j ) ( ui - uj )2 / ( vi + vj )
	fdr = p = 0.0;
	for( i = 0; i < m; i++ ){			// for each class pair
		pu = u + i * d;
		pv = v + i * d;
		for( j = i; j < m; j++ ){
			pu2 = u + j * d;
			pv2 = v + j * d;
			for( k = 0; k < d; k++ ) {	// for each dimension
				temp = pu[k] - pu2[k];
				fdr += SQR( temp );
				p += pv[k] + pv2[k];
			}
		}
	}
	// free memory pool
	free( mem );
	return fdr / p;
}

//
// PCA
//

int pca( double *pc, double *var, double *x, int n, int d ){
	double *r, *u, *px, *pr, *eig, p, temp;
	int i, j, k;
	void *mem;

	// memory pool
	mem = calloc( ( d + 1 ) * d, sizeof( double ) );
	u = ( double* )mem;		// d
	r = u + d;				// d * d
	eig = u;
	if( mem == NULL ) return 1;

	// compute mean
	px = x;
	for( i = 0; i < n; i++ ){
		for( j = 0; j < d; j++, px++ ){
			u[j] += px[0];
		}
	}
	p = 1.0 / ( double )n;
	for( i = 0; i < d; i++ ) u[i] *= p;

	// compute correlation matrix R
	for( i = 0; i < d; i++ ){				// Rij
		pr = r + i * d + i;
		for( j = i; j < d; j++, pr++ ){			// i-d and j-d of all x
			px = x;
			for( k = 0; k < n; k++, px += d )
				// ( xi-u ) - ( xj-u )
				pr[0] += ( px[i] - u[i] ) * ( px[j] - u[j] );
		}
	}
	// complete R
	for( i = 0; i < d; i++ ){
		pr = px = r + i * d + i;
		for( j = i; j < d; j++, pr++, px += d ) px[0] = pr[0] = pr[0] * p;
	}
	// pca
	Jacobi( r, eig, pc, d, 30 );
	// sorting by eigenvalue from large to small
	for( i = 0; i < d; i++ ){
		for( j = i + 1, k = i; j < d; j++ ){// k is max index
			if( eig[j] > eig[k] ) k = j;
		}
		var[i] = eig[k];
		if( i == k ) continue;
		// swap eigenvalue
		temp = eig[k];
		eig[k] = eig[i];
		eig[i] = temp;
		px = pc + i;
		pr = pc + k;
		// swap column
		for( j = 0; j < d; j++, px += d, pr += d ){
			temp = px[0];
			px[0] = pr[0];
			pr[0] = temp;
		}
	}
	return 0;
}

int pca( double *pc, double *x, int n, int d, int *l, const double rv ){
	double *w, *px, *pr, p, temp;
	int i, j;
	void *mem;

	// memory pool
	mem = calloc( ( d + 1 ) * d, sizeof( double ) );
	w = ( double* )mem;		// d
	if( mem == NULL ) return 1;

	if( pca( pc, w, x, d, n ) ) {
		free( mem );
		return 1;
	}

	p = 0.0;
	for( i = 0; i < d; i++ ) p += w[i];
	// select top k pc with reasonable variance > rv
	p = 1.0 / p;
	temp = 0.0;
	for( i = 0; i < d; ) {
		temp += w[ i++ ];
		if( temp * p > rv ) break;
	}
	*l = i;
	// solution saving
	px = pc;
	for( i = 0; i < d; i++ ){
		pr = pc + i * d;
		for( j = 0; j < *l; j++, pr++, px++ ) px[0] = pr[0];
	}
	// free memory pool
	free( mem );
	return 0;
}

int pca( double *pc, double *x, int n, int d, int l ){
	double *a, *ppc, *pa;
	int i, j;

	// memory pool
	a = ( double* )malloc( d * d * sizeof( double ) );
	if( a == NULL ) return 1;
	if( pca( a, x, d, n, &i, 1.0 ) ) {
		free( a );
		return 1;
	}

	// copy solution
	ppc = pc;
	for( i = 0; i < d; i++ ){
		pa = a + i * d;
		for( j = 0; j < l; j++, pa++, ppc++ ) ppc[0] = pa[0];
	}
	free( a );
	return 0;
}

//
// PDA: Parametrix discriminant analysis 
//

int pda( double *a, double *x, int m, int n, int d ){
	double *sb, *sw, *s, *psb, *psw, *ps;
	int i, j, k;
	void *mem;

	// memory pool
	mem = malloc( 3 * d * d * sizeof( double ) );
	sb = ( double* )mem;		// d * d
	sw = sb + d * d;			// d * d
	s = sw + d * d;				// d * d
	if( mem == NULL ) return 1;

	// compute scatter matrix
	ScatterB( sb, x, m, n, d );
	ScatterW( sw, x, m, n, d );
	inverse( sw, d, s );

	// S = Sw^-1 * Sb
	ps = s;
	for( i = 0; i < d; i++ ){
		for( j = 0; j < d; j++, ps++ ){
			ps[0] = 0.0;
			psw = sw + i * d;
			psb = sb + j;
			for( k = 0; k < d; k++, psw++, psb += d )
				ps[0] += psw[0] * psb[0];
		}
	}
	// eigenvector of S
	Jacobi( s, sb, a, d, 50 );
	// free memory pool
	free( mem );
	return 0;
}

//
// NDA: Nonparametrix discriminant analysis 
//

int nda( double *a, double *x, int m, int n, int d, int knn, double param ){
	double *sb, *sw, *psb, *psw, *ps;
	int i, j, k;
	void *mem;

	// memory pool
	mem = malloc( 2 * d * d * sizeof( double ) );
	sb = ( double* )mem;		// d * d
	sw = sb + d * d;			// d * d
	if( mem == NULL ) return 1;

	// compute scatter matrix
	NPScatterB( sb, x, m, n, d, knn, param );
	ScatterW( sw, x, m, n, d );
	inverse( sw, d, a );

	// S = Sw^-1 * Sb
	ps = a;
	for( i = 0; i < d; i++ ){
		for( j = 0; j < d; j++, ps++ ){
			ps[0] = 0.0;
			psw = sw + i * d;
			psb = sb + j;
			for( k = 0; k < d; k++, psw++, psb += d )
				ps[0] += psw[0] * psb[0];
		}
	}
	// eigenvector of S
	Jacobi( a, sb, sw, d, 50 );
	for( i = 0, j = d * d; i < j; i++ ) a[i] = sw[i];
	// free memory pool
	free( mem );
	return 0;
}
