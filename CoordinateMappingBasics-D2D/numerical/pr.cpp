
#include "pr.h"
#include <numerical\linear.h>
#include <numerical\linear2.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <stdio.h>
#include <string.h>
using namespace distribution;

// 
// macro
//

#define SQR(x) ( (x) * (x) )
#define BIG		10000000000

//
// inline functions
//

inline void MultiMatrix( double *a, double *b, double *c, int n ){
	int i, j, k;
	for( i = 0; i < n; i++ ){
		for( j = 0; j < n; j++ ){
			c[i * n + j] = 0.0;
			for( k = 0; k < n; k++ ) c[i * n + j] += a[i * n + k] * b[k * n + j];
		}
	}
}

inline void MultiVector( double *a, double *b, double *c, int n ){
	int i, j;
	for( i = 0; i < n; i++ ){
		c[i] = 0.0;
		for( j = 0; j < n; j++ ) c[i] += a[i * n + j] * b[j];
	}
}

inline void Transpose( double *a, int n ){
	int i, j;
	double temp;
	for( i = 0; i < n; i++ ){
		for( j = i + 1; j < n; j++ ){
			temp = a[i * n + j];
			a[i * n + j] = a[j * n + i];
			a[j * n + i] = temp;
		}
	}
}

inline double Dot( double *a, double *b, int n ){
	double temp = 0.0;
	int i;
	for( i = 0; i < n; i++ ) temp += a[i] * b[i];
	return temp;
}

inline double logfunc( double x ){
	return 1.0 / ( 1.0 + exp( -x ) );
}

//
// Function declarations
//

void CompSbSm( double *sb, double *sm, int d, int m, va_list );

// 
// Structures
//

void NewLClassifier( LClassifier *w, int d ){
	w->w = ( double* )malloc( ( d + 1 ) * sizeof( double ) );
	w->d = d;
}

void DeleteLClassifier( LClassifier w ){
	free( w.w );
}

void NewLClassifierf( LClassifierf *w, int d ){
	w->w = ( float* )malloc( ( d + 1 ) * sizeof( float ) );
	w->d = d;
}

void DeleteLClassifierf( LClassifierf w ){
	free( w.w );
}

void NewNeuralNet( NeuralNet *net, int d, int k ){
	net->w = ( double* )malloc( ( ( d + 3 ) * k + 1 ) * sizeof( double ) );
	net->d = d;
	net->k = k;
}

void DeleteNeuralNet( NeuralNet net ){
	free( net.w );
}

void NewAdaBoost( AdaBoost *ab, int d, int k ){
	ab->w = ( double* )malloc( ( d + 2 ) * k * sizeof( double ) ); 
	ab->a = ab->w + ( d + 1 ) * k;
	ab->d = d;
	ab->k = k;
}

void DeleteAdaBoost( AdaBoost ab ){
	free( ab.w );
}

//
// Classification
//

void Bayesian( Normal_ND &Ni, Normal_ND &Nj, double Pi, double Pj, double x[], double w[] ){
	int i, j, d = Ni.GetDimension();
	double temp, *var, *ui, *uj;
	void *mem;
	
	// memory pool
	mem = malloc( ( 2 * d + 2 * d * d ) * sizeof( double ) );
	var = ( double* )mem + d * d;
	ui = var + d * d;
	uj = ui + d;

	// initialize
	Ni.Mean( ui );
	Nj.Mean( uj );
	Ni.Var( var );
	inverse( var, d, mem );

	for( i = 0; i < d; i++ ) x[i] = ui[i] - uj[i];		// ( ui - uj )
	for( i = 0; i < d; i++ ){
		w[i] = 0.0;
		// dot w
		for( j = 0; j < d; j++ ) w[i] += var[i * d + j] * x[j];		
	}

	// -1 / [ (ui - uj)T * V^-1 * (ui - uj) ] * ln( Pi / Pj ) * (ui - uj) 
	temp = 0.0;
	for( i = 0; i < d; i++ ) temp += x[i] * w[i];	// dot
	temp = -( log( Pi / Pj ) / temp );
	for( i = 0; i < d; i++ ) x[i] *= temp;
	// generate x
	for( i = 0; i < d; i++ ) x[i] = 0.5 * ( ui[i] + uj[i] ) + x[i];
	
	free( mem );		// free memory pool
}

int Euclidean( Normal_ND N[], double x[], int c ){
	double *u, temp, dist, min =DBL_MAX;
	int i, j, ans = -1, d = N[0].GetDimension();

	// initialize
	u = ( double* )malloc( d * sizeof( double ) );

	for( i = 0; i < c; i++ ){
		// initialize
		N[i].Mean( u );
		dist = 0.0;
		// calculate euclidean distance
		for( j = 0; j < d; j++ ) {
			temp = u[j] - x[j];
			dist += SQR( temp );
		}
		if( dist < min ){
			min = dist;
			ans = i;
		}
	}

	free( u );
	return ans;
}

int Mahalanobis( Normal_ND N[], double x[], int c ){
	double *u, *var, *temp, dist, min =DBL_MAX;
	void *mem;
	int i, j, ans = -1, d = N[0].GetDimension();

	// initialize
	mem = malloc( ( 2 * d + 2 * d * d ) * sizeof( double ) );
	u = ( double* )mem + d * d;
	temp = u + d;
	var = temp + d;

	for( i = 0; i < c; i++ ){
		// initialize
		N[i].Mean( u );
		N[i].Var( var );
		inverse( var, d, mem );
		// calculate mahalanobis distance
		dist = 0.0;
		for( j = 0; j < d; j++ ) u[j] = x[j] - u[j];	// ( x - u )
		MultiVector( var, u, temp, d );
		for( j = 0; j < d; j++ ) dist += u[j] * temp[j];
		if( dist < min ){
			min = dist;
			ans = i;
		}
	}

	free( mem );
	return ans;
}

//
// Linear Classifier
//

//
// perceptron
//

void Perceptronf( LClassifierf *c, FeatureArrayf x, ClassArrayf y, ArrayNum n, LearningRatef lr, int time ){
	FeatureDim d = ( *c ).d;
	float temp, *w = ( *c ).w;
	int i, j, k;

	for( k = 0; k < time; k++ ){
		for( i = 0; i < n; i++ ){
			temp = y[i] * lr( k );
			for( j = 0; j < d; j++ ) w[j] += temp * x[i][j];
			w[j] += temp;
		}
	}
}

void Perceptronf( LClassifierf *c, TrueFeatureSetf x1, FalseFeatureSetf x2, LearningRatef lr, int time ){
	FeatureDim d = ( *c ).d;
	float temp, *w = ( *c ).w;
	int i, j, k;

	for( k = 0; k < time; k++ ){
		temp = lr( k );
		for( i = 0; i < x1.n; i++ ){
			for( j = 0; j < d; j++ ) w[j] += temp * x1.x[i][j];
			w[j] += temp;
		}
		for( i = 0; i < x2.n; i++ ){
			for( j = 0; j < d; j++ ) w[j] -= temp * x2.x[i][j];
			w[j] -= temp;
		}
	}
}

void Perceptron( double x[], double y[], LClassifier *c, int nx, double p, int time ){
	double *px, *w, temp;
	int i, j, k, d = c->d;

	w = c->w;
	for( i = 0; i <= d; i++ ) w[i] = 0.0;		// initialize
	for( k = 0; k < time; k++ ){
		for( i = 0; i < nx; i++ ){
			px = x + i * d;
			temp = y[i] * p;
			for( j = 0; j < d; j++ ) w[j] += temp * px[j];
			w[j] -= temp;
		}
	}
	// normalize
	temp = 0.0;
	for( i = 0; i <= d; i++ ) temp += w[i] * w[i];
	temp = sqrt( temp );
	for( i = 0; i <= d; i++ ) w[i] /= temp;
}

//
// LMS Algorithm
//

void LMS_Classifyf( LClassifierf *c, FeatureArrayf x, ClassArrayf y, ArrayNum n, LearningRatef lr, int time ){
	FeatureDim d = ( *c ).d;
	float temp, p, *w = ( *c ).w;
	int i, j, k;

	for( k = 0; k < time; k++ ){
		p = lr( k );
		for( i = 0; i < n; i++ ){
			// y(x)
			temp = 0.0f;
			for( j = 0; j < d; j++ ) temp += x[i][j] * w[j];
			temp += w[j];
			// w += lr * x * ( y - y(x) )
			temp = p * ( y[i] - temp );
			for( j = 0; j < d; j++ ) w[j] += temp * x[i][j];
			w[j] += temp;
		}
	}
}

void LMS_Classifyf( LClassifierf *c, TrueFeatureSetf x1, FalseFeatureSetf x2, LearningRatef lr, int time ){
	FeatureDim d = ( *c ).d;
	float temp, p, *w = ( *c ).w;
	int i, j, k;

	for( k = 0; k < time; k++ ){
		p = lr( k );
		// for data set 1
		for( i = 0; i < x1.n; i++ ){
			// y(x)
			temp = 0.0f;
			for( j = 0; j < d; j++ ) temp += x1.x[i][j] * w[j];
			temp += w[j];
			// w += lr * x * ( y - y(x) )
			temp = p * ( 1.0f - temp );
			for( j = 0; j < d; j++ ) w[j] += temp * x1.x[i][j];
			w[j] += temp;
		}
		// for data set 2
		for( i = 0; i < x2.n; i++ ){
			// y(x)
			temp = 0.0f;
			for( j = 0; j < d; j++ ) temp += x2.x[i][j] * w[j];
			temp += w[j];
			// w += lr * x * ( y - y(x) )
			temp = p * ( -1.0f - temp );
			for( j = 0; j < d; j++ ) w[j] += temp * x2.x[i][j];
			w[j] += temp;
		}
	}
}

void LMS_Train( double x[], double y[], LClassifier *c, int nx, double p ){
	double *px, *w, temp;
	int i, j, d = c->d;

	w = c->w;
	for( i = 0; i <= d; i++ ) w[i] = 1.0;		// initialize
	for( i = 0; i < nx; i++ ){
		px = x + i * d;
		// xkT * u( k - 1 )
		temp = 0.0;
		for( j = 0; j < d; j++ ) temp += px[j] * w[j];
		temp += w[j];
		// p( y - xkT * u( k - 1 ) )
		temp = p * ( y[i] - temp );
		for( j = 0; j < d; j++ ) w[j] += temp * px[j];
		w[j] += temp;
		// normalize
		temp = 0.0;
		for( j = 0; j <= d; j++ ) temp += w[j] * w[j];
		temp = sqrt( temp );
		for( j = 0; j <= d; j++ ) w[j] /= temp;
	}
}

//
// Least Square Classifier
//

void LS_Classifyf( LClassifierf *c, FeatureArrayf x, ClassArrayf y, ArrayNum n ){
	Matf xTx, inv;
	Vecf w = ( *c ).w;
	FeatureDim d = ( *c ).d, d1 = d + 1;
	float *px;
	void *mem;
	int i, j, k;

	// zerolize w
	memset( w, 0x00, d1 * sizeof( float ) );

	// memory pool
	mem = malloc( SQR( d1 ) * 2 * sizeof( float ) );
	if( mem == NULL ) return;
	xTx = ( float * )mem;
	inv = xTx + SQR( d1 );

	// compute xTx
	px = xTx;
	memset( xTx, 0x00, SQR( d1 ) * sizeof( float ) );
	for( i = 0; i < d; i++ ){
		for( j = 0; j < d; j++, px++ ){
			for( k = 0; k < n; k++ ){
				px[0] += x[k][i] * x[k][j];
			}
		}
		for( k = 0; k < n; k++ ){
			px[0] += x[k][i];
		}
		px++;
	}
	// last row
	for( i = 0; i < d; i++, px++ ){
		for( j = 0; j < n; j++ ){
			px[0] += x[j][i];
		}
	}
	px[0] = ( float )n;

	// inv( xTx )
	invf( inv, xTx, d1 );
	// xT * y
	memset( xTx, 0x00, d1 * sizeof( float ) );
	for( i = 0; i < d; i++ ){
		for( j = 0; j < n; j++ ){
			xTx[i] += x[j][i] * y[j];
		}
	}
	for( j = 0; j < n; j++ ) xTx[i] += y[j]; 
	// inv( xTx ) * ( xT * y )
	px = inv;
	for( i = 0; i <= d; i++ ){
		for( j = 0; j <= d; j++, px++ ){
			w[i] += px[0] * xTx[j];
		}
	}
	free( mem );
}

void LS_Classifyf( LClassifierf *c, TrueFeatureSetf x1, FalseFeatureSetf x2 ){
	Matf xTx, inv;
	Vecf w = ( *c ).w;
	FeatureDim d = ( *c ).d, d1 = d + 1;
	float *px;
	void *mem;
	int i, j, k;

	// zerolize w
	memset( w, 0x00, d1 * sizeof( float ) );

	// memory pool
	mem = malloc( SQR( d1 ) * 2 * sizeof( float ) );
	if( mem == NULL ) return;
	xTx = ( float * )mem;
	inv = xTx + SQR( d1 );

	// compute xTx
	px = xTx;
	memset( xTx, 0x00, SQR( d1 ) * sizeof( float ) );
	for( i = 0; i < d; i++ ){
		for( j = 0; j < d; j++, px++ ){
			for( k = 0; k < x1.n; k++ ) px[0] += x1.x[k][i] * x1.x[k][j];
			for( k = 0; k < x2.n; k++ ) px[0] += x2.x[k][i] * x2.x[k][j];
		}
		for( k = 0; k < x1.n; k++ ) px[0] += x1.x[k][i];
		for( k = 0; k < x2.n; k++ ) px[0] += x2.x[k][i];
		px++;
	}
	// last row
	for( i = 0; i < d; i++, px++ ){
		for( j = 0; j < x1.n; j++ ) px[0] += x1.x[j][i];
		for( j = 0; j < x2.n; j++ ) px[0] += x2.x[j][i];
	}
	px[0] = ( float )( x1.n + x2.n );

	// inv( xTx )
	invf( inv, xTx, d1 );
	// xT * y
	memset( xTx, 0x00, d1 * sizeof( float ) );
	for( i = 0; i < d; i++ ){
		for( j = 0; j < x1.n; j++ ) xTx[i] += x1.x[j][i];
		for( j = 0; j < x2.n; j++ ) xTx[i] -= x2.x[j][i];
	}
	xTx[i] = ( float )( x1.n - x2.n );
	// inv( xTx ) * ( xT * y )
	px = inv;
	for( i = 0; i <= d; i++ ){
		for( j = 0; j <= d; j++, px++ ){
			w[i] += px[0] * xTx[j];
		}
	}
	free( mem );
}

void LS_Train( double x[], double y[], LClassifier *c, int nx ){
	double *X, *px, *pX, *w = c->w;
	int i, j, d = c->d;

	X = ( double* )malloc( nx * ( d + 1 ) * sizeof( double ) );

	px = x;	pX = X;
	for( i = 0; i < nx; i++ ){
		for( j = 0; j < d; j++, px++, pX++ ) *pX = *px;
		*( pX++ ) = 1.0;
	}

	LS_method( X, y, w, nx, d + 1 );

	free( X );
}

//
// Classify Function
//

float Classifyf( Featuref x, LClassifierf c ){
	float g = 0.0f;
	int i, d = c.d;
	for( i = 0; i < d; i++ ) g += x[i] * c.w[i];
	return g + c.w[d];
}

double Classify( double x[], LClassifier c ){
	double g = 0.0;
	int i, d = c.d;
	for( i = 0; i < d; i++ ) g += x[i] * c.w[i];
	return g + c.w[d];
}

//
// Nonlinear Classifier
//

// Neural Network, 2LP model

double BackProp( double x[], double y[], NeuralNet *net, int nx, 
			  double lr, double ri, double rd, double c, int time ){
	int i, j, t, it, d = net->d, k = net->k;
	double *pw, *px, *py, *a, *v, *pv, *pa, *yy, *pyy, *w = net->w,
		e, e_prev, temp;
	void *mem;

	// w has ( d + 1 ) * k + ( k + 1 ) elements
	for( i = ( d + 1 ) * k + k; i >= 0; i-- ) w[i] = ( double )rand() / RAND_MAX; 
	mem = malloc( ( 3 * nx * ( k + 1 ) ) * sizeof( double ) );
	v = ( double* )mem;			// ( k + 1 ) * n
	a = v + ( k + 1 ) * nx;		// ( k + 1 ) * n
	yy = a + ( k + 1 ) * nx;	// ( k + 1 ) * n
	e = 0.0;

	for( it = 0; it < time; it++ ){
		e_prev = e;
		e = 0.0;
		px = x;	py = y;	pv = v;	pa = a;	pyy = yy;
		for( t = 0; t < nx; t++, px += d, py++, pv += k + 1, pa += k + 1, pyy += k + 1 ){
			// forward computation
			// input layer
			pw = w;
			for( i = 0; i < k; i++, pw += d + 1 ){
				pv[i] = 0.0;
				for( j = 0; j < d; j++ ) pv[i] += px[j] * pw[j];
				pv[i] += pw[d];
				pyy[i] = logfunc( pv[i] );
			}
			// first
			pv[k] = 0.0;
			for( i = 0; i < k; i++ ) pv[k] += pyy[i] * pw[i];
			pv[k] += pw[k];
			pyy[k] = logfunc( pv[k] );
			// least square energy function
			if( py[0] > 0.0 && pyy[k] > 0.0 ) e += pyy[k] * log( py[0] / pyy[k] );
			if( py[0] < 1.0 && pyy[k] < 1.0 ) e += ( 1.0 - pyy[k] ) * log( ( 1.0 - py[0] ) / ( 1.0 - pyy[k] ) );
			
			// backword computation
			pa[k] = ( pyy[k] - py[0] ) * pyy[k] * ( 1.0 - pyy[k] );
			for( i = 0; i < k; i++ ) {
				pa[i] = pyy[i] * ( 1.0 - pyy[i] ) * ( pw[i] * pa[k] ); 
			}
		}
		if( e < 0.0000001 ) break;	// stop critiria
		// update learning rate
		temp = e / e_prev;
		if( temp < 1.0 ) lr *= ri;
		else if( temp > c ) lr *= rd;
		// update 
		px = x;	pa = a;	pv = v;	pyy = yy;
		for( t = 0; t < nx; t++, px += d, pv += k + 1, pa += k + 1, pyy += k + 1 ){
			// input layer
			pw = w;
			for( i = 0; i < k; i++, pw += d + 1 ){
				for( j = 0; j < d; j++ ) pw[j] -= lr * pa[i] * px[j];
				pw[d] -= lr * pa[i];
			}
			// update first layer
			for( i = 0; i < k; i++ ) pw[i] -= lr * pa[k] * pyy[i];
			pw[k] -= lr * pa[k];
		}
	}
	free( mem );
	return e;
}

double Classify( double x[], NeuralNet net ){
	int i, j, d = net.d, k = net.k;
	double ans = 0.0, *y, *pw, *w = net.w;

	y = w + ( d + 2 ) * k + 1; 

	// input layer
	pw = w;
	for( i = 0; i < k; i++, pw += d + 1 ){
		y[i] = 0.0;
		for( j = 0; j < d; j++ ) y[i] += x[j] * pw[j];
		y[i] += pw[d];
		y[i] = logfunc( y[i] );
	}
	// first layer
	for( i = 0; i < k; i++ ) ans += y[i] * pw[i];
	ans = logfunc( ans + pw[k] );

	return ans;
}

// AdaBoost

void Boosting( double x[], double y[], AdaBoost *ab, int n ){
	int i, j, s, d = ab->d, k = ab->k;
	double *px, *w, *pw, *g, *a = ab->a, p;
	void *mem;
	LClassifier h;

	// memory pool
	mem = malloc( 2 * n * sizeof( double ) );
	w = ( double* )mem;								// n
	g = w + n;										// n

	// initialize
	p = 1.0 / ( double )n;
	for( i = 0; i < n; i++ ) w[i] = p;
	
	h.d = d;
	h.w = pw = ab->w;
	for( i = 0; i < k; i++, pw += d + 1 ){
		/**/
		// LMS weighted training
		for( j = 0; j <= d; j++ ) pw[j] = ( double )rand() / RAND_MAX;		// initialize
		px = x;
		for( j = 0; j < n; j++, px += d ){
			// skip the data that is correctly classified 
			if( i > 0 && y[j] * g[j] > 0.0 ) continue;	
			// xkT * u( k - 1 )
			p = 0.0;
			for( s = 0; s < d; s++ ) p += px[s] * pw[s];
			p += pw[s];
			// wi( y - xkT * u( k - 1 ) )
			p = w[j] * ( y[j] - p );
			for( s = 0; s < d; s++ ) pw[s] += p * px[s];
			pw[s] += p;
			// normalize
			p = 0.0;
			for( s = 0; s <= d; s++ ) p += pw[s] * pw[s];
			p = sqrt( p );
			for( s = 0; s <= d; s++ ) pw[s] /= p;
		}
		/**/
		h.w = pw;
		// calculate p = sum( w(i) ) for misclassified data
		px = x;
		p = 0.0;
		for( j = 0; j < n; j++, px += d ){
			g[j] = Classify( px, h );
			if( y[j] * g[j] < 0.0 ) p += w[j];
		}
		if( p == 1.0 ) a[i] = 0.0;
		else if( p > 0.0 ) a[i] = 0.5 * log( ( 1.0 - p ) / p );			// get a
		else {
			// single linear classifier
			for( j = 0; j <= d; j++ ) ab->w[j] = h.w[j];
			ab->k = 1;
			a[0] = 1.0;
			break;
		}
		// update w and training data
		p = 0.0;
		for( j = 0; j < n; j++ ){
			w[j] *= ( y[j] * g[j] > 0.0 ) ? exp( -a[i] ) : exp( a[i] );
			p += w[j];
		}
		// normalize w
		for( j = 0; j < n; j++ ) w[j] /= p;
	}
	// free memory
	free( w );
}

double Classify( double x[], AdaBoost ab ){
	int i, j, d = ab.d, k = ab.k;
	double *w = ab.w, *a = ab.a, ans, temp;

	ans = 0.0;
	for( i = 0; i < k; i++, w += d + 1 ){
		temp = 0.0;
		for( j = 0; j < d; j++ ) temp += x[j] * w[j];
//		ans += a[i] * ( temp + w[j] );
		ans += ( temp + w[j] > 0.0 ) ? a[i] : -a[i];
	}
	return ans;
}

//
// Estimation
//

// Expectation-Maximization

void EM_Estimate( double x[], int nx, int d, int n, Mixture *mix, int time ){
	double *p, *uj, *vj, *xk, *px, *pjk, temp, sum;
	void *mem;
	int i, j, k, a, b, it;

	// memory pool
	mem = malloc( ( d * ( 2 + d ) + n * ( nx + 1 ) ) * sizeof( double ) );
	uj = ( double* )mem;	// d
	p = uj + d;				// n
	vj = p + n;				// d * d
	pjk = vj + d * d;		// nx * n
	xk = pjk + nx * n;		// d

	// initialize
	for( i = 0; i < n; i++ ) p[i] = mix->Prob( i );

	for( it = 0; it < time; it++ ){
		// Update Pjk
		for( k = 0; k < nx; k++ ){			// nx
			sum = 0.0;
			// j-th model
			for( j = 0; j < n; j++ ){
				for( i = 0; i < d; i++ ) xk[i] = x[ k * d + i ];
				sum += mix->f( j, xk ) * p[j];
			}
			for( j = 0; j < n; j++ ){
				for( i = 0; i < d; i++ ) xk[i] = x[ k * d + i ];
				pjk[ j * nx + k ] = mix->f( j, xk ) * p[j] / sum;
			}
		}		

		// for each mixed gaussian models	
		for( j = 0; j < n; j++ ){
			// initialize
			for( i = 0; i < d; i++ ) uj[i] = 0.0;
			k = d * d;
			for( i = 0; i < k; i++ ) vj[i] = 0.0;
			sum = 0.0;
			for( k = 0; k < nx; k++ ) sum += pjk[ j * nx + k ];

			// mean
			for( k = 0; k < nx; k++ ){
				px = x + k * d;
				temp = pjk[j * nx + k] / sum;
				for( i = 0; i < d; i++ ) uj[i] += px[i] * temp;
			}
			mix->SetMean( j, uj );

			// variance
			for( k = 0; k < nx; k++ ){
				px = x + k * d;
				temp = pjk[j * nx + k] / sum;
				for( a = 0; a < d; a++ ){
					for( b = a; b < d; b++ ){
						vj[a * d + b] += temp * ( px[a] - uj[a] ) * ( px[b] - uj[b] );
						vj[b * d + a] = vj[a * d + b];
					}
				}
			}
			mix->SetVar( j, vj );
			
			// probability
			p[j] = sum / ( double )nx;
		}
		mix->SetProb( p );
	}

	free( mem );
}

// Parsen Window

double ParzenWin( double xi[], double x[], int nx, int d, double h ){
	Normal_ND N( d );
	double p = 0.0, *y;
	int i, j;

	y = ( double* )calloc( d, sizeof( double ) );

	for( i = 0; i < nx; i++ ){
		for( j = 0; j < d; j++ ) y[j] = ( x[j] - xi[i * d + j] ) / h;		// x - xi
		p += N.f( y );
	}

	p /= ( double )nx;
	for( i = 0; i < d; i++ ) p /= h;

	free( y );
	return p;
}
