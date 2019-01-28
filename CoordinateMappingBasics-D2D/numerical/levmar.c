#include "linear2.h"
#include "nonlinear.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

#include <stdio.h>

#define SMALL		0.0001
#define SQR( x )	( (x) * (x) )

//
// globle vatiables
//

__declspec( thread ) static LevMarFunc2f	_levmar2f;
__declspec( thread ) static LevMarFunc2d	_levmar2d;
__declspec( thread ) static Mat2f			_yf;
__declspec( thread ) static Mat2d			_yd;
__declspec( thread ) static Dim			_row, _col;

//
// residual function
//

static Errf Resf( LevMarFuncf, Vecf_in, Vecf_out, Dim_in, Dim_out );
static Errd Resd( LevMarFuncd, Vecd_in, Vecd_out, Dim_in, Dim_out );

//
// local function warper
//

static void GoalFunc2f( Vecf_in, Vecf_out );
static void GoalFunc2d( Vecd_in, Vecd_out );

//
// y = f( x ), x, y are vector
// y are known, solve x
//

Errf LevMarf( LevMarFuncf goal, Vecf_unknown x, Vecf_value y, Dim_in nx, Dim_out ny, 
			 IterTime time, Paramf tau, Errf e1, Errf e2, Errf e3 ){
	float u = ( float )fabs( tau ), v = 2.0f;
	float q, temp;
	float *pj, *p, *pn;
	float 
		*x_new,			// nx
		*dx,			// nx
		*g,				// nx
		*e,				// ny
		*e_new,			// ny
		*y_head,		// ny
		*y2,			// ny
		*jacobian,		// nx x ny
		*n,				// nx^2
		*ninv;			// nx^2
	void *mem;
	int i, j, k, it;

	mem = malloc( ( 3 * nx + 4 * ny + 2 * SQR( nx ) + nx * ny ) * sizeof( float ) );
	if( mem == NULL ) return FLT_MAX;
	x_new = ( float * )mem;
	dx = x_new + nx;
	g = dx + nx;
	e = g + nx;
	e_new = e + ny;
	y_head = e_new + ny;
	y2 = y_head + ny;
	jacobian = y2 + ny;
	n = jacobian + nx * ny;
	ninv = n + SQR( nx );

	// estimate by initial guess
	goal( x, y_head );
	// compute error
	for( i = 0; i < ny; i++ ) e[i] = y[i] - y_head[i];
	// copy x solution
	for( i = 0; i < nx; i++ ) x_new[i] = x[i];
	// compute jacobian matrix
	temp = 1.0f / ( float )SMALL;
	for( j = 0; j < nx; j++, pj++ ){
		pj = jacobian + j;
		x_new[j] += ( float )SMALL;
		if( j != 0 ) x_new[ j - 1 ] = x[ j - 1 ];
		goal( x_new, y2 );
		for( i = 0; i < ny; i++, pj += nx ) pj[0] = ( y2[i] - y_head[i] ) * temp;
	}
	// generate N = JT * J matrix
	memset( n, 0x00, SQR( nx ) * sizeof( float ) );
	for( i = 0; i < nx; i++ ){
		pn = n + i * nx + i;
		for( j = i; j < nx; j++, pn++ ){
			pj = jacobian + i;
			p = jacobian + j;
			for( k = 0; k < ny; k++, pj += nx, p += nx ) {
				pn[0] += pj[0] * p[0];
			}
			n[ j * nx + i ] = pn[0];
		}
	}
	// set parameter u and v
	pn = n;
	for( i = 0, u = n[0]; i < nx; i++, pn += nx + 1 ) if( pn[0] > u ) u = pn[0]; 
	u *= tau; v = 2.0f;

	// g = JT * e
	temp = 0.0f;
	memset( g, 0x00, nx * sizeof( float ) );
	for( i = 0; i < nx; i++ ){
		pj = jacobian + i;
		for( j = 0; j < ny; j++, pj += nx ) g[i] += pj[0] * e[j];
		temp += SQR( g[i] );
	}

	// exit condition 1
	if( temp < SQR( e1 ) ){
		free( mem );
		return Resf( goal, x, y, nx, ny );
	}

	// start iteration
	for( it = 0; it < time; it++ ){
		while(1){
			// generate ( N + uI );
			pn = n;
			for( i = 0; i < nx; i++, pn += nx + 1 ) pn[0] += u;
			// solve ( N + uI )dx = g;
			temp = 0.0;
			invf( ninv, n, nx );
			memset( dx, 0x00, nx * sizeof( float ) );
			pn = ninv;
			for( i = 0; i < nx; i++ ){
				for( j = 0; j < nx; j++, pn++ ) dx[i] += pn[0] * g[j];
				// compute for stop condition 2
				temp += SQR( dx[i] ) - SQR( e2 ) * SQR( x[i] );
			}
			// stop condition 2
			if( temp < 0.0f ){
				free( mem );
				return Resf( goal, x, y, nx, ny );
			}
			// compute new x
			for( i = 0; i < nx; i++ ) x_new[i] = x[i] + dx[i];
			// and estimate new y_head
			q = 0.0f;					// error trend
			goal( x_new, y_head );
			for( i = 0; i < ny; i++ ){	// compute error
				e_new[i] = y[i] - y_head[i];
				q += SQR( e[i] ) - SQR( e_new[i] );	// compute error trend
			}

			// normalize error trend by devite dxT( udx + g )
			temp = 0.0;
			for( i = 0; i < nx; i++ ) temp += dx[i] * ( u * dx[i] + g[i] );
			q /= temp;

			// break this iteration if error decrease
			if( q > 0.0f ){
				memcpy( x, x_new, nx * sizeof( float ) );	// update x
				// compute jacobian matrix
				temp = 1.0f / ( float )SMALL;
				for( j = 0; j < nx; j++, pj++ ){
					pj = jacobian + j;
					x_new[j] += ( float )SMALL;
					if( j != 0 ) x_new[ j - 1 ] = x[ j - 1 ];
					goal( x_new, y2 );
					for( i = 0; i < ny; i++, pj += nx ) pj[0] = ( y2[i] - y_head[i] ) * temp;
				}
				// generate N = JT * J matrix
				memset( n, 0x00, SQR( nx ) * sizeof( float ) );
				for( i = 0; i < nx; i++ ){
					pn = n + i * nx + i;
					for( j = i; j < nx; j++, pn++ ){
						pj = jacobian + i;
						p = jacobian + j;
						for( k = 0; k < ny; k++, pj += nx, p += nx ) {
							pn[0] += pj[0] * p[0];
						}
						n[ j * nx + i ] = pn[0];
					}
				}
				// update error
				memcpy( e, e_new, ny * sizeof( float ) );
				temp = 0.0f;
				for( i = 0; i < ny; i++ ) temp += SQR( e[i] );
				// stop condition 3
				if( temp < e3 ){
					free( mem );
					return Resf( goal, x, y, nx, ny );
				}
				// compte g
				temp = 0.0f;
				memset( g, 0x00, nx * sizeof( float ) );
				for( i = 0; i < nx; i++ ){
					pj = jacobian + i;
					for( j = 0; j < ny; j++, pj += nx ) g[i] += pj[0] * e[j];
					temp += SQR( g[i] );
				}
				// stop condition 1
				if( temp < SQR( e1 ) ){
					free( mem );
					return Resf( goal, x, y, nx, ny );
				}

				// update u, v
				temp = 2.0f * q - 1.0f;
				v = 1.0f - temp * SQR( temp );
				u *= ( 0.33333f > v ? 0.33333f : v );
				v = 2.0f;
				break;
			}
			else u *= v; v *= 2;
		}
	}
	free( mem );
	return Resf( goal, x, y, nx, ny );
}

Errd LevMard( LevMarFuncd goal, Vecd_unknown x, Vecd_value y, Dim_in nx, Dim_out ny, 
			 IterTime time, Paramd tau, Errd e1, Errd e2, Errd e3 ){
	double u = fabs( tau ), v = 2.0;
	double q, temp;
	double *pj, *p, *pn;
	double 
		*x_new,			// nx
		*dx,			// nx
		*g,				// nx
		*e,				// ny
		*e_new,			// ny
		*y_head,		// ny
		*y2,			// ny
		*jacobian,		// nx x ny
		*n,				// nx^2
		*ninv;			// nx^2
	void *mem;
	int i, j, k, it;

	mem = malloc( ( 3 * nx + 4 * ny + 2 * SQR( nx ) + nx * ny ) * sizeof( double ) );
	if( mem == NULL ) return DBL_MAX;
	x_new = ( double * )mem;
	dx = x_new + nx;
	g = dx + nx;
	e = g + nx;
	e_new = e + ny;
	y_head = e_new + ny;
	y2 = y_head + ny;
	jacobian = y2 + ny;
	n = jacobian + nx * ny;
	ninv = n + SQR( nx );

	// estimate by initial guess
	goal( x, y_head );
	// compute error
	for( i = 0; i < ny; i++ ) e[i] = y[i] - y_head[i];
	// copy x solution
	for( i = 0; i < nx; i++ ) x_new[i] = x[i];
	// compute jacobian matrix
	temp = 1.0 / SMALL;
	for( j = 0; j < nx; j++, pj++ ){
		pj = jacobian + j;
		x_new[j] += SMALL;
		if( j != 0 ) x_new[ j - 1 ] = x[ j - 1 ];
		goal( x_new, y2 );
		for( i = 0; i < ny; i++, pj += nx ) pj[0] = ( y2[i] - y_head[i] ) * temp;
	}
	// generate N = JT * J matrix
	memset( n, 0x00, SQR( nx ) * sizeof( double ) );
	for( i = 0; i < nx; i++ ){
		pn = n + i * nx + i;
		for( j = i; j < nx; j++, pn++ ){
			pj = jacobian + i;
			p = jacobian + j;
			for( k = 0; k < ny; k++, pj += nx, p += nx ) {
				pn[0] += pj[0] * p[0];
			}
			n[ j * nx + i ] = pn[0];
		}
	}
	// set parameter u and v
	pn = n;
	for( i = 0, u = n[0]; i < nx; i++, pn += nx + 1 ) if( pn[0] > u ) u = pn[0]; 
	u *= tau; v = 2.0;

	// g = JT * e
	temp = 0.0f;
	memset( g, 0x00, nx * sizeof( double ) );
	for( i = 0; i < nx; i++ ){
		pj = jacobian + i;
		for( j = 0; j < ny; j++, pj += nx ) g[i] += pj[0] * e[j];
		temp += SQR( g[i] );
	}

	// exit condition 1
	if( temp < SQR( e1 ) ){
		free( mem );
		return Resd( goal, x, y, nx, ny );
	}

	// start iteration
	for( it = 0; it < time; it++ ){
		while(1){
			// generate ( N + uI );
			pn = n;
			for( i = 0; i < nx; i++, pn += nx + 1 ) pn[0] += u;
			// solve ( N + uI )dx = g;
			temp = 0.0;
			invd( ninv, n, nx );
			memset( dx, 0x00, nx * sizeof( double ) );
			pn = ninv;
			for( i = 0; i < nx; i++ ){
				for( j = 0; j < nx; j++, pn++ ) dx[i] += pn[0] * g[j];
				// compute for stop condition 2
				temp += SQR( dx[i] ) - SQR( e2 ) * SQR( x[i] );
			}
			// stop condition 2
			if( temp < 0.0f ){
				free( mem );
				return Resd( goal, x, y, nx, ny );
			}
			// compute new x
			for( i = 0; i < nx; i++ ) x_new[i] = x[i] + dx[i];
			// and estimate new y_head
			q = 0.0;					// error trend
			goal( x_new, y_head );
			for( i = 0; i < ny; i++ ){	// compute error
				e_new[i] = y[i] - y_head[i];
				q += SQR( e[i] ) - SQR( e_new[i] );	// compute error trend
			}

			// normalize error trend by devite dxT( udx + g )
			temp = 0.0;
			for( i = 0; i < nx; i++ ) temp += dx[i] * ( u * dx[i] + g[i] );
			q /= temp;

			// break this iteration if error decrease
			if( q > 0.0 ){
				memcpy( x, x_new, nx * sizeof( double ) );	// update x
				// compute jacobian matrix
				temp = 1.0f / SMALL;
				for( j = 0; j < nx; j++, pj++ ){
					pj = jacobian + j;
					x_new[j] += SMALL;
					if( j != 0 ) x_new[ j - 1 ] = x[ j - 1 ];
					goal( x_new, y2 );
					for( i = 0; i < ny; i++, pj += nx ) pj[0] = ( y2[i] - y_head[i] ) * temp;
				}
				// generate N = JT * J matrix
				memset( n, 0x00, SQR( nx ) * sizeof( double ) );
				for( i = 0; i < nx; i++ ){
					pn = n + i * nx + i;
					for( j = i; j < nx; j++, pn++ ){
						pj = jacobian + i;
						p = jacobian + j;
						for( k = 0; k < ny; k++, pj += nx, p += nx ) {
							pn[0] += pj[0] * p[0];
						}
						n[ j * nx + i ] = pn[0];
					}
				}
				// update error
				memcpy( e, e_new, ny * sizeof( double ) );
				temp = 0.0f;
				for( i = 0; i < ny; i++ ) temp += SQR( e[i] );
				// stop condition 3
				if( temp < e3 ){
					free( mem );
					return Resd( goal, x, y, nx, ny );
				}
				// compte g
				temp = 0.0f;
				memset( g, 0x00, nx * sizeof( double ) );
				for( i = 0; i < nx; i++ ){
					pj = jacobian + i;
					for( j = 0; j < ny; j++, pj += nx ) g[i] += pj[0] * e[j];
					temp += SQR( g[i] );
				}
				// stop condition 1
				if( temp < SQR( e1 ) ){
					free( mem );
					return Resd( goal, x, y, nx, ny );
				}

				// update u, v
				temp = 2.0 * q - 1.0;
				v = 1.0 - temp * SQR( temp );
				u *= ( 0.33333 > v ? 0.33333 : v );
				v = 2.0;
				break;
			}
			else u *= v; v *= 2;
		}
	}
	free( mem );
	return Resd( goal, x, y, nx, ny );
}

Errf LevMar2f( LevMarFunc2f goal, Vecf_unknown x, Mat2f_value y, Dim_in nx, Row_out my, Col_out ny, 
			  IterTime time, Paramf tau, Errf e1, Errf e2, Errf e3 ){
	int i, n;
	float ret;
	Vecf vecy, py;

	// set thread variables
	_levmar2f = goal;
	_yf = y;
	_row = my;
	_col = ny;
	// memory allocate and copy
	n = my * ny;
	py = vecy = ( float* )malloc( n * sizeof( float ) );
	for( i = 0; i < my; i++, py += ny ){
		memcpy( py, y[i], ny * sizeof( float ) );
	}
	// exccute
	ret = LevMarf( GoalFunc2f, x, vecy, nx, n, time, tau, e1, e2, e3 );
	// copy back
	py = vecy;
	for( i = 0; i < my; i++, py += ny ){
		memcpy( y[i], py, ny * sizeof( float ) );
	}
	free( vecy );
	return ret;
}

Errd LevMar2d( LevMarFunc2d goal, Vecd_unknown x, Mat2d_value y, Dim_in nx, Row_out my, Col_out ny, 
			  IterTime time, Paramd tau, Errd e1, Errd e2, Errd e3 ){
	int i, n;
	double ret;
	Vecd vecy, py;

	// set thread variables
	_levmar2d = goal;
	_yd = y;
	_row = my;
	_col = ny;
	// memory allocate and copy
	n = my * ny;
	py = vecy = ( double* )malloc( n * sizeof( double ) );
	for( i = 0; i < my; i++, py += ny ){
		memcpy( py, y[i], ny * sizeof( double ) );
	}
	// exccute
	ret = LevMard( GoalFunc2d, x, vecy, nx, n, time, tau, e1, e2, e3 );
	// copy back
	py = vecy;
	for( i = 0; i < my; i++, py += ny ){
		memcpy( y[i], py, ny * sizeof( double ) );
	}
	free( vecy );
	return ret;
}

//
// resudual
//

Errf Resf( LevMarFuncf goal, Vecf_in x, Vecf_out y, Dim_in nx, Dim_out ny ){
	float temp, e = 0.0f, *y_head;
	int i;

	y_head = ( float * )malloc( ny * sizeof( float ) );

	// compute
	goal( x, y_head );
	// residual
	for( i = 0; i < ny; i++ ){
		temp = y_head[i] - y[i];
		e += SQR( temp );
	}

	// finally
	free( y_head );
	return ( float )sqrt( e );
}

Errd Resd( LevMarFuncd goal, Vecd_in x, Vecd_out y, Dim_in nx, Dim_out ny ){
	double temp, e = 0.0f, *y_head;
	int i;

	y_head = ( double * )malloc( ny * sizeof( double ) );

	// compute
	goal( x, y_head );
	// residual
	for( i = 0; i < ny; i++ ){
		temp = y_head[i] - y[i];
		e += SQR( temp );
	}

	// finally
	free( y_head );
	return sqrt( e );
}

//
// Goal Function
//

void GoalFunc2f( Vecf_in x, Vecf_out y ){
	int i;
	_levmar2f( x, _yf );		// goal function
	// copy
	for( i = 0; i < _row; i++, y += _col ){
		memcpy( y, _yf[i], _col * sizeof( float ) );
	}
}

void GoalFunc2d( Vecd_in x, Vecd_out y ){
	int i;
	_levmar2d( x, _yd );		// goal function
	// copy
	for( i = 0; i < _row; i++, y += _col ){
		memcpy( y, _yf[i], _col * sizeof( double ) );
	}
}
