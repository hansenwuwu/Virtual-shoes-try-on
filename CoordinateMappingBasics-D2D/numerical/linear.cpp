#include "linear.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>

//
// macro
//

#define SQR( x ) ( (x) * (x) )
#define min( x, y ) ( ( (x) < (y) ) ? (x) : (y) )
#define max( x, y ) ( ( (x) > (y) ) ? (x) : (y) )
#define ERR 0.0000000000001

//
// Linear Algebra
//

double det( double *A, int n, void *mem ){	
	double ans = 1.0;
	double ratio, max, temp;
	bool interchange=false;
	int i, j, k;
	int row;
	double *M;

	// memory pool
	if( mem == NULL ) M = ( double* )malloc( n * n * sizeof( double ) );
	else M = ( double* )mem;
	for( i = 0; i < n * n; i++ ) M[i] = A[i];	// copy matrix
	for(i = 0; i < n-1; i++){
		row = i;
		max = 0.0;
		for( j = i; j < n; j++ ){
			if( fabs( M[j * n + i] ) > fabs( max ) ) {
				max = M[j * n + i];
				row = j;
			}
		}
		if( M[row * n + i] == 0.0 ) continue;	//all downside value is 0
		if( row != i ){
			for( j = 0; j < n; j++ ) {
				temp = M[i * n + j];
				M[i * n + j] = M[row * n + j];
				M[row * n + j] = temp;
			}
			interchange=!interchange;
		}
		for( j = i + 1; j < n; j++ ) {
			if( M[j * n + i] == 0.0 ) continue;
			ratio = M[j * n + i] / M[i * n + i];
			for( k = 0; k < n; k++ ) M[j * n + k] -= ( M[i * n + k] * ratio );
		}
	}
	for(i=0; i<n; i++) ans*=M[i * n + i];

	if( mem == NULL ) free( M );
	if(interchange) return -ans;
	else return ans;
}

int inverse( double *A, double *inv, int n, void *mem ){
	double *M, max, ratio;
	int i, j, k, row;

	// memory pool
	if( mem == NULL ) M = ( double* )malloc( n * n * sizeof( double ) );
	else M = ( double* )mem;
	for( i = 0; i < n; i++ ) {
		for( j = 0; j < n; j++ ){
			M[i * n + j] = A[i * n + j];	// copy matrix
			inv[i * n + j] = ( i == j ) ? 1.0 : 0.0;
		}
	}

	for(j=0; j<n; j++){
		max=0.0;
		// find the pivot value
		for(i=j; i<n; i++){
			if(fabs(M[i * n + j])>fabs(max)) {
				max=M[i * n + j];
				row=i;
			}
		}
		// check if equal to 0
		if( max == 0.0 ) return 1;
		// row exchange if pivot isn't a diogonal element
		if( row != j ){
			for( k = 0; k < n; k++ ) {
				// swap
				max = M[row * n + k];
				M[row * n + k] = M[j * n + k];
				M[j * n + k] = max;
				max = inv[row * n + k];
				inv[row * n + k] = inv[j * n + k];
				inv[j * n + k] = max;
			}
		}
		// row operations
		for( i = 0; i < n; i++ ){
			if( i == j ){
				ratio = 1.0 / M[i * n + j];
				for( k = 0; k < n; k++ ){
					M[i * n + k] *= ratio;
					inv[i * n + k] *= ratio;
				}
			}
			else{
				if( M[i * n + j] == 0.0 ) continue;
				ratio = M[i * n + j] / M[j * n + j];
				for( k = 0; k < n; k++ ){
					M[i * n + k] -= ratio * M[j * n + k];
					inv[i * n + k] -= ratio * inv[j * n + k];
				}
			}
		}
	}

	if( mem == NULL ) free( M );	// free memory pool
	return 0;
}

int inverse( double *inv, int n, void *mem ){
	double *M, max, ratio;
	int i, j, k, row;

	// memory pool
	if( mem == NULL ) M = ( double* )malloc( n * n * sizeof( double ) );
	else M = ( double* )mem;
	for( i = 0; i < n; i++ ) {
		for( j = 0; j < n; j++ ){
			M[i * n + j] = inv[i * n + j];	// copy matrix
			inv[i * n + j] = ( i == j ) ? 1.0 : 0.0;
		}
	}

	for(j=0; j<n; j++){
		max=0.0;
		// find the pivot value
		for(i=j; i<n; i++){
			if(fabs(M[i * n + j])>fabs(max)) {
				max=M[i * n + j];
				row=i;
			}
		}
		// check if equal to 0
		if( max == 0.0 ) return 1;
		// row exchange if pivot isn't a diogonal element
		if(row!=j){
			for(k=0; k<n; k++) {
				// swap
				max = M[row * n + k];
				M[row * n + k] = M[j * n + k];
				M[j * n + k] = max;
				max = inv[row * n + k];
				inv[row * n + k] = inv[j * n + k];
				inv[j * n + k] = max;
			}
		}
		// row operations
		for(i=0; i<n; i++){
			if(i==j){
				ratio=1.0/M[i * n + j];
				for(k=0; k<n; k++){
					M[i * n + k]*=ratio;
					inv[i * n + k]*=ratio;
				}
			}
			else{
				if(M[i * n + j]==0.0) continue;
				ratio=M[i * n + j]/M[j * n + j];
				for(k=0; k<n; k++){
					M[i * n + k]-=ratio*M[j * n + k];
					inv[i * n + k]-=ratio*inv[j * n + k];
				}
			}
		}
	}

	if( mem == NULL ) free( M );	// free memory pool
	return 0;
}

// Gauss-Jordan Elimination
int GaussJordan( double *M, double *b, int n ){
	double max, ratio;
	int i, j, k, row;

	for( j = 0; j < n; j++ ){		// for every column
		max=0.0;
		// find the pivot value
		for( i = j; i < n; i++ ){	// for every element in column
			if( fabs( M[i * n + j]) > fabs(max) ) {
				max = M[i * n + j];
				row = i;
			}
		}
		// check if equal to 0
		if( max == 0.0 ) return 1;
		// row exchange if pivot isn't a diogonal element
		if( row != j ){
			for(k=0; k<n; k++) {
				// swap
				max = M[row * n + k];
				M[row * n + k] = M[j * n + k];
				M[j * n + k] = max;
			}
			// swap
			max = b[row];
			b[row] = b[j];
			b[j] = max;
		}
		// row operations
		ratio = 1.0 / M[j * n + j];
		for( k = 0; k < n; k++ ) M[j * n + k] *= ratio;
		b[j] *= ratio;
		// elimination
		for( i = 0; i < n; i++ ){	// for every row below j
			if( i == j || M[i * n + j] == 0.0 ) continue;
			ratio = M[i * n + j];
			for( k = j; k < n; k++ ) M[i * n + k] -= ratio * M[j * n + k];
			b[i] -= ratio * b[j];
		}
	}

	return 0;
}

// Gram-Schmidt Algorithm

void Gram_Schmidt( double *A, double *Q, int row, int col ){
	double temp;
	int i, j, k;

	// initialize
	for( i = 0; i < row * col; i++ ) Q[i] = A[i];

	// Processing
	for( j = 0; j < col; j++ ){			// for each column
		for( k = j - 1; k >= 0; k-- ){	// for the other bases
			// inner product
			temp = 0.0;
			for( i = 0; i < row; i++ ) temp += Q[i * col + j] * Q[i * col + k];
			// subtract the projects to previer bases
			for( i = 0; i < row; i++ ) Q[i * col + j] -= temp * Q[i * col + k];
		}
		// unitize
		temp = 0.0;
		for( i = 0; i < row; i++ ) temp += SQR( Q[i * col + j] );
		if( temp == 0.0 ) continue;		// this column is spaned by previer bases
		temp = sqrt( temp );
		for( i = 0; i < row; i++ ) Q[i * col + j] /= temp;
	}
}


// SIngular Value Decomposition

void SVD( double *A, double *U, double *W, double *V, int row, int col, void *mem ){
	double *memory, *ATA, *eig, *eig_vec, *u, *v;	// u: row, v: col
	double temp, sum;
	int i, j, k;

	// memory pool allocate
	if( mem == NULL ){			// check is memory pool inputed
		memory = ( double * )malloc( ( 2 * col * col + 2 * col + row ) * sizeof( double ) );
	}
	else memory = ( double* )mem;
	ATA = memory;				// col * col
	eig_vec = ATA + col * col;	// col * col
	eig = eig_vec + col * col;	// col
	v = eig + col;				// the temp vector: col
	u = v + col;				// the temp vector: row

	// initialize ATA
	for( i = 0; i < col; i++ ){			// for every row of ATA
		for( j = i; j < col; j++ ){		// for column: [i, n]
			temp = 0.0;
			for( k = 0; k < row; k++ ) temp += A[k * col + i] * A[k * col + j];
			ATA[i * col + j] = ATA[j * col + i] = temp;
		}
	}

	// Solve eigen problem
	Jacobi( ATA, eig, eig_vec, col );

	// sort eigenvalue and eigenvector
	for( i = 0; i < col; i++ ){				// for every column vector
		for( j = i + 1; j < col; j++ ){		// for every column vector unless vector i 
			if( eig[i] < eig[j] ){			// compare the eigenvalue and swap the value and vector
				// swap( eig[i], eig[j] )
				temp = eig[i];
				eig[i] = eig[j];
				eig[j] = temp;
				for( k = 0; k < col; k++ ) {
					// swap( eig_vec[k][i], eig_vec[k][j] );
					temp = eig_vec[k * col + i];
					eig_vec[k * col + i] = eig_vec[k * col + j];
					eig_vec[k * col + j] = temp;
				}
			}
		}
	}
	// save V: col * col = eigenvalue matrix
	for( i = 0; i < col; i++ ){
		for( j = 0; j < col; j++ ) V[i * col + j] = eig_vec[i * col + j];
	}

	// save W: row * col
	for( i = 0; i < row; i++ ){			// for every row
		for( j = 0; j < col; j++ ){		// for every column
			if( i == j ){				// diogonal element
				if( eig[i] > ERR ) W[i * col + i] = sqrt( eig[i] );
				else W[i * col + i] = 0.0;		// too small
			}
			else W[i * col + j] = 0.0;			// set the other to 0
		}
	}

	// save U: row * row, size of u: row, size of v: col
	for( i = 0; i < row; i++ ){
		for( j = 0; j < row; j++ ){
			if( i > col ) v[j] = 0.0;
			else v[j] = V[j * col + i];
		}
		if( i < col ){
			sum = 0.0;						// initialize the length of vector u
			for( j = 0; j < row; j++ ){		// for every element of u
				temp = 0.0;
				for( k = 0; k < col; k++ )	temp += A[j * col + k] * v[k];	// multi matrix
				u[j] = temp;				// save solution
				sum += temp * temp;			// length square
			}
			sum = sqrt( sum );				// length
			for( j = 0; j < row; j++ ) U[j * row + i] = u[j] / sum;	// unitize vector u and save as column vector of U
		}
		else {
			for( j = 0; j < i; j++ ) U[ j * row + i ] = 0.0;
			U[ i * row + i ] = 1.0;
			for( j = i + 1; j < row; j++ ) U[ j * row + i ] = 0.0;
		}
	}

	// free memory pool if it's builded in this function
	if( mem == NULL ) free( memory );
}

void svd( double *x, double *u, double *w, double *v, int row, int col, void *mem ){
	double *xx, *xi, *xj, *px, *pv, *pu, temp;
	int n, i, j, k;

	// memory pool
	n = ( row < col ) ? row : col;			// n = min( row, col )
	if( mem == NULL ) xx = ( double* )calloc( n * n, sizeof( double ) );
	else{
		xx = ( double* )mem;
		for( i = 0, j = n * n; i < j; i++ ) xx[i] = 0.0;
	}
	
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
		Jacobi( xx, w, u, n );
		for( i = 0; i < n; i++ ) w[i] = ( w[i] > ERR ) ? sqrt( w[i] ) : 0.0;
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
		// XTX = XXT, n*n matrix
		if( row == col ){
			for( i = 0, j = n * n; i < j; i++ ) u[i] = v[i];
			if( mem == NULL ) free( xx );
			return;
		}
		// v = XTu
		pv = v;
		for( i = 0; i < col; i++ ){				// XT * ui
			for( j = 0; j < row; j++, v++ ){	// u1 ... un
				pu = u + j;
				px = x + i;
				pv[0] = 0.0;
				for( k = 0; k < row; k++, pu += n, px += col ) pv[0] += pu[0] * px[0];
				pv[0] /= w[j];			// normalize by eigenvalue
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
		Jacobi( xx, w, v, n );
		for( i = 0; i < n; i++ ) w[i] = ( w[i] > ERR ) ? sqrt( w[i] ) : 0.0;
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
		}
		// u = Xv
		pu = u;
		for( i = 0; i < row; i++ ){				// X * vi
			for( j = 0; j < col; j++, pu++ ){	// v1 ... vn
				pv = v + j * col;
				px = x + i * col;
				pu[0] = 0.0;
				for( k = 0; k < col; k++, pv++, px++ ) pu[0] += pv[0] * px[0];
				pu[0] /= w[j];			// normalize by eigenvalue
			}
		}
	}
	// free memory pool
	if( mem == NULL ) free( xx );
}

// QR factorization

void QRfact( double *A, double *Q, double *R, int row, int col ){
	double temp;
	int i, j, k;

	// initialize
	for( i = 0; i < row * col; i++ ) Q[i] = A[i];
	for( i = 0; i < col * col; i++ ) R[i] = 0.0;

	// Processing
	for( j = 0; j < col; j++ ){			// for each column
		for( k = j - 1; k >= 0; k-- ){	// for the other bases
			// inner product
			temp = 0.0;
			for( i = 0; i < row; i++ ) temp += Q[i * col + j] * Q[i * col + k];
			// subtract the projects to previer bases
			for( i = 0; i < row; i++ ) Q[i * col + j] -= temp * Q[i * col + k];
			// set to R
			R[k * col + j] = temp;
		}
		// unitize
		temp = 0.0;
		for( i = 0; i < row; i++ ) temp += SQR( Q[i * col + j] );
		if( temp == 0.0 ) continue;		// this column is spaned by previer bases
		temp = sqrt( temp );
		for( i = 0; i < row; i++ ) Q[i * col + j] /= temp;
		// set to R's diagonal element
		R[j * col + j] = temp;
	}
}

// LU factorization

int LUfact( double *A, int *index, int n, double *mem ){
	double sum, big, *v, *a, *b;
	int i, j, k, row;
	void *memory;

	// memory pool
	if( mem == NULL ) memory = malloc( n * sizeof( double ) );
	else memory = mem;
	v = ( double * )memory;

	// calculate largest value per row
	a = A;
	for( i = 0; i < n; i++ ){
		big = 0.0;
		for( j = 0; j < n; j++, a++ ){
			if( ( sum = fabs( *a ) ) > big ) big = sum;
		}
		if( big == 0.0 ) return 1;
		v[i] = 1.0 / big;
	}
	for( j = 0; j < n; j++ ){
		// calculate U when i < j
		for( i = 0; i < j; i++ ){
			sum = A[ i * n + j ];
			for( k = 0; k < i; k++ ) sum -= A[ i * n + k ] * A[ k * n + j ];
			A[ i * n + j ] = sum;
		}
		big = 0.0;
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
	// free memory pool
	if( mem == NULL ) free( memory );
	return 0;
}

void SolveLU( double *A, int *ind, double *b, int n ){
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

void LUcomb( double *A, int *ind, int n ){
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

// Least-Square Method

void LS_method( double *A, double *B, double *X, int row, int col, void *mem ){
	double *U, *W, *V, *UX, *VX;			// SVD matrix: A = U * W * V
	double *memory, temp;
	int i, j, n = min( row, col );
	
	// memory pool
	if( mem == NULL ){
		memory = ( double* )malloc( ( row + col + SQR( row ) + SQR( col ) + row * col +
			2 * col * col + 2 * col + row ) * sizeof( double ) );
	}
	else memory = ( double* )mem;
	UX = memory;			// row
	VX = memory + row;		// col
	U = VX + col;			// row * row
	V = U + row * row;		// col * col
	W = V + col * col;		// row * col

	SVD( A, U, W, V, row, col, W + row * col );
	for( i = 0; i < n; i++ ){
		if( fabs( W[i * col + i] ) < ERR || fabs( W[i * col + i] ) > 1.0 / ERR )
			W[i * col + i] = 0.0;		// set 0 is too small or too big
		else W[i * col + i] = 1.0 / W[i * col + i];
	}
	// compute UX = UT * B
	for( i = 0; i < row; i++ ){
		temp = 0.0;
		for( j = 0; j < row; j++ ) temp += U[j * row + i] * B[j];
		UX[i] = temp;
	}
	// compute VX = WT * UX
	for( i = 0; i < row; i++ ){
		temp = 0.0;
		for( j = 0; j < col; j++ ) temp += W[j * col + i] * UX[j];
		VX[i] = temp;
	}
	// compute X = V * VX
	for( i = 0; i < col; i++ ){
		temp = 0.0;
		for( j = 0; j < col; j++ ) temp += V[i * col + j] * VX[j];
		X[i] = temp;
	}
	
	// free memory pool if it's builded in this function
	if( mem == NULL ) free( memory );
}

//
//	Eigensystem
//

void Jacobi( double *A, double *eigen, double *v, int n, int time ){
	int i, j, p, q;
	double prev=0, sum, tangent, cosine, sine, tau, ptemp, qtemp,
		*pp, *pq, *qq, *a, *b;

	if(n<=1) return;

	//initialize A, V, eigen
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			if(i==j) v[i * n + j]=1.0;
			else v[i * n + j]=0.0;
		}
		eigen[i]=A[i * n + i];
	}
	//jacobi
	for( i = 0; i < time; i++, prev = sum ){
//		printf( "%d-th times. ", i + 1 );
		sum = 0.0;
		for( p = 0; p < n-1; p++ ){
			a = A + p * n + p + 1;
			for( q = p + 1; q < n; q++, a++ ){
				sum += fabs( *a );
			}
		}
//		printf( "sum: %lf\n", sum );
		if( sum == 0.0 || ( sum >= prev && i > 1 ) ) break;
		for( p = 0; p < n - 1; p++ ){
			pp = A + p * n + p;
			for( q = p + 1; q < n; q++ ){
				qq = A + q * n + q;
				pq = A + p * n + q;
				if( *pq == 0.0 && *qq - *pp == 0.0 ) tangent = tan( atan( 2.0 ) / 2.0 );
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

double PowerIter( double *A, double *u, int n, int time, double *memory ){
	double *a, *v, max = 0.0, prev = DBL_MAX;
	int i, j, k;
	void *mem;

	// initialize guess and memory pool
	if( memory == NULL ) mem = malloc( n * sizeof( double ) );
	else mem = memory;
	v = ( double * )mem;
	for( i = 0; i < n; i++ ) u[i] = ( double )rand() / ( double )RAND_MAX;

	for( k = 0; k < time; k++ ){
		if( fabs( max - prev ) < ERR ) break;
		// calculate v
		for( i = 0; i < n; i++ ){
			v[i] = 0.0;
			a = A + i * n;
			for( j = 0; j < n; j++ ) v[i] += u[j] * a[j];
		}
		// calculate u
		prev = max;
		max = 0.0;
		for( i = 0; i < n; i++ ) {
			if( fabs( v[i] ) > max ) max = fabs( v[i] );
		}
		if( max > ERR ) for( i = 0; i < n; i++ ) u[i] = v[i] / max;
	}
	if( memory == NULL ) free( mem );
	return max;
}

double InvIter( double *A, double *b0, double shift, int n, int time, double *memory ){
	double *b1, sum, g, prev;
	int *ind, i, k;
	void *mem;

	// memory pool
	if( memory == NULL ) mem = malloc( n * ( sizeof( double ) + sizeof( int ) ) );
	else mem = memory;
	b1 = ( double * )mem;
	ind = ( int* )( b1 + n );

	// initialize
	for( i = 0; i < n; i++ ) A[ i * n + i ] -= shift;
	LUfact( A, ind, n, b1 );

	// keep guess when grown factor too small
	do{
		for( i = 0; i < n; i++ ) b0[i] = b1[i] = ( double )rand() / ( double )RAND_MAX;
		// normalize for both
		sum = 0.0;
		for( i = 0; i < n; i++ ) sum += SQR( b0[i] );
		sum = sqrt( sum );
		for( i = 0; i < n; i++ ) b0[i] = b1[i] = b0[i] / sum;
		// solve ( A- aI )b1 = b0
		SolveLU( A, ind, b1, n );
		// normalized b1
		g = 0.0;					// grown factor: | b0 |
		for( i = 0; i < n; i++ ) g += SQR( b1[i] );
		g = sqrt( g );
	} while( g < ERR );	
	
	// normalize
	prev = 0.0;							// | b1 - b0 |
	for( i = 0; i < n; i++ ) {	
		b1[i] /= g;						// nomralize
		prev += SQR( b0[i] - b1[i] );	
		b0[i] = b1[i];					// updata b0
	}

	for( k = 0; k < time; k++ ){
		// solve ( A- aI )y = b
		SolveLU( A, ind, b1, n );
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
		if( g < ERR ) break;
		// updata shift value when shift of | b1 - b0 | < error bound
		if( prev - g < ERR ){
			prev = g;
			g = 0.0;
			for( i = 0; i < n; i++ ) g += b0[i] * b1[i];
			g = 1.0 / ( g * sum );	// unnormalized y 
			LUcomb( A, ind, n );
			for( i = 0; i < n; i++ ) A[ i * n + i ] += shift;
			shift += g;
			for( i = 0; i < n; i++ ) A[ i * n + i ] -= shift;
			LUfact( A, ind, n, b0 );
		}
		else prev = g;
		// updata b0
		for( i = 0; i < n; i++ ) b0[i] = b1[i];
	}

	if( memory == NULL ) free( mem );
	return shift;
}

void ReductHessen( double *A, int n ){		// reduce a balanced matrix to Hessenberg form
	int m, i, j;
	double x, y, *a, *b;

	for( m = 1; m < n - 1; m++ ){			// r+1 
		x = 0.0;
		i = m;
		a = A + m * n + m - 1;
		for( j = m; j < n; j++, a += n ){			// find pivot
			if( fabs( *a ) > fabs( x ) ){
				x = *a;
				i = j;
			}
		}
		if( i != m ){
			// row interchange
			a = A + i * n + m - 1;
			b = A + m * n + m - 1;
			for( j = m - 1; j < n; j++, a++, b++ ){
				y = *a;
				*a = *b;
				*b = y;
			}
			// column interchange
			a = A + i;
			b = A + m;
			for( j = 0; j < n; j++, a += n, b += n ){
				y = *a;
				*a = *b;
				*b = y;
			}
		}
		if( x ){		// elimination
			for( i = m + 1; i < n; i++ ){
				if( ( y = A[ i * n + m - 1 ] ) != 0.0 ){
					y /= x;
					// subtract row
					a = A + i * n + m - 1;
					b = A + m * n + m - 1;
					for( j = m - 1; j < n; j++, a++, b++ ) *a -= y * *b;
					// add column
					a = A + m;
					b = A + i;
					for( j = 0; j < n; j++, a += n, b += n ) *a += y * *b;
				}
			}
		}
	}
}

void HessenQR( double *A, int n, double *wr, int time ){
	int nn, m, l, k, j, its, i, mmin;
	double x, y, z, w, u, v, t, s, r, q, p, anorm,
		*a, *b;

	anorm = 0.0;
	for( i = 0; i < n; i++ ){	// comput norm of matrix
		a = A + i * n;			// start from i - 1 or 0 ( hessenberg form )
		if( i > 0 ) {
			j = i - 1;
			a += j;
		}
		else j = 0;
		for( ; j < n; j++, a++ ) anorm += fabs( *a );
	}
	nn = n - 1;
	t = 0.0;				// Gets changed only by an exceptional shift
	while( nn >= 0 ){		// Search for next eigenvalue
		its = 0;
		do{
			// begin itaration: look for single small subdiagonal element
			l = nn;
			a = A + l * n + l;
			b = A + ( l - 1 ) * n + l - 1;
			for( /*l = nn*/; l >= 1; l--, a -= n + 1, b -= n + 1 ){
				s = fabs( *a ) + fabs( *b );
				if( s == 0.0 ) s = anorm;
				if( fabs( *( a - 1 ) ) + s == s ) break;
			}
			x = a[ nn * n + nn ];
			if( l == nn ) {				// one root found
				wr[ nn-- ] = x + t;
			}
			else{
				y = A[ ( nn - 1 ) * n + ( nn - 1 ) ];
				w = A[ nn * n + nn - 1 ] * A[ ( nn - 1 ) * n + nn ];
				if( l == nn - 1 ){		// two root found...
					p = 0.5 * (y - x );
					q = p * p + w;
					z = sqrt( fabs( q ) );
					x += t;
					if( q >= 0.0 ){		// a real pair
						z = p + ( ( p > 0.0 ) ? fabs( z ) : -fabs( z ) );
						wr[ nn - 1 ] = wr[ nn ] = x + z;
						if( z ) wr[ nn ] = x - w / z;
					}
					else{				// a complex pair
						wr[ nn - 1 ] = wr[ nn ] = x + p;
					}
					nn -= 2;
				}
				else{								// no roots found, continue iteration
					if( its == time ) return;
					if( its % 10 == 0 ){			// form exceptional shift
						t += x;
						a = A;
						for( i = 0; i < nn; i++, a += n + 1 ) *a -= x;
						s = fabs( A[ nn * n + nn - 1 ] ) + fabs( A[ ( nn - 1 ) * n + nn - 2 ] );
						y = x = 0.75 * s;
						w = -0.4375 * s * s;
					}
					its++;
					m = nn - 2;
					a = A + m * n + m;
					for( /*m = nn - 2*/; m >= l; m--, a -= n + 1 ){	// form shift and then 
						z = *a;										// look for 2 consecutive small subdiagonal element
						r = x - z;
						s = y - z;
						p = ( r * s - w ) / *( a + n ) + *( a + 1 );
						q = *( a + n + 1 ) - z - r - s;
						r = *( a + 2 * n + 1 );
					}
				}
			}
		} while( l < nn - 1 );
	}
}
