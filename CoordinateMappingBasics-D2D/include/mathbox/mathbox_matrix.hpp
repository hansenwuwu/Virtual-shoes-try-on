
#ifndef __MATHBOX_MATRIX__
#define __MATHBOX_MATRIX__

#include <math.h>
#include <string>


int		mElementMultiple	( const int *matrix_a, const int *matrix_b, int n );
float	mElementMultiple	( const float *matrix_a, const float *matrix_b, int n );
void	mMultiple			( int *martix, float multiple, int n );
void	mMultiple			( float *martix, float multiple, int n );
void	mMultiple			( const int *martix_left, const int *matrix_right, int *matrix_out, int n );
void	mMultiple			( const float *martix_left, const float *matrix_right, float *matrix_out, int n );
void	mPow				( const float *matrix_in, float *matrix_out, int n, int pow );
int		mDet				( int *matrix, int n );
float	mDet				( float *matrix, int n );
int		mDet				( int *const matrix, int n, int row, int col );
float	mDet				( float *const matrix, int n, int row, int col );
void	mTranspose			( int *matrix_in, int *matrix_out, int n );
void	mTranspose			( float *matrix_in, float *matrix_out, int n );
bool	mInverse			( int *const matrix_in, int *matrix_out, int n );
bool	mInverse			( float *const matrix_in, float *matrix_out, int n );


/*	Do
________________________________________________________________*/

int mRowMultCol( const int *_row, const int *_col, int _n )
{
	static int temp, i;
	temp = 0;
	for( i = 0; i < _n; i++, _row++, _col += _n )
		temp += _row[0]*_col[0];
	return temp;
}

float mRowMultCol( const float *_row, const float *_col, int _n )
{
	static int i;
	static float temp;
	temp = 0.f;
	for( i = 0; i < _n; i++, _row++, _col += _n )
		temp += _row[0]*_col[0];
	return temp;
}

int mElementMultiple( const int *_a, const int *_b, int _n )
{
	static int temp;
	temp = 0;
	for( _n *= _n; _n > 0; _n--, _a++, _b++ )
		temp += _a[0]*_b[0];
	return temp;
}

float mElementMultiple( const float *_a, const float *_b, int _n )
{
	static float temp;
	temp = 0.f;
	for( _n *= _n; _n > 0; _n--, _a++, _b++ )
		temp += _a[0]*_b[0];
	return temp;
}

void mMultiple( int *_m, float _r, int _n )
{
	_n *= _n;
	for( ; _n > 0; _n--, _m++ )
		_m[0] = (int)( _m[0]*_r );
}

void mMultiple( float *_m, float _r, int _n )
{
	_n *= _n;
	for( ; _n > 0; _n--, _m++ )
		_m[0] *= _r;
}

void mMultiple( const int *_l, const int *_r, int *_out, int _n )
{
	static int r, c;
	memset( _out, 0, _n*_n*sizeof(int) );
	for( r = 0; r < _n; r++ )
		for( c = 0; c < _n; c++, _out++ )
			_out[0] = mRowMultCol( _r + _n*r, _l + c, _n );
}

void mMultiple( const float *_l, const float *_r, float *_out, int _n )
{
	static int r, c;
	memset( _out, 0, _n*_n*sizeof(float) );
	for( r = 0; r < _n; r++ )
		for( c = 0; c < _n; c++, _out++ )
			_out[0] = mRowMultCol( _r + _n*r, _l + c, _n );
}

void mPow( const float *_s, float *_out, int _n, int _pow )
{
	int n = _n*_n;
	static float *m;
	m = (float*)malloc( n*sizeof(float) );
	memcpy( m, _s, n*sizeof(float) );
	for( _pow -= 1; _pow > 0; _pow-- )
	{
		mMultiple( _s, m, _out, _n );
		memcpy( m, _out, n*sizeof(float) );
	}
}

int mDet( int *_m, int _n )
{
	static int sum;
	sum = 0;
	switch( _n )
	{
		case 1:
			return _m[0];

		case 2:
			sum += _m[0]*_m[3];
			sum += -_m[1]*_m[2];
			return sum;

		case 3:
			sum += _m[0]*_m[4]*_m[8];
			sum += _m[1]*_m[5]*_m[6];
			sum += _m[2]*_m[3]*_m[7];
			sum += -_m[2]*_m[4]*_m[6];
			sum += -_m[1]*_m[3]*_m[8];
			sum += -_m[0]*_m[5]*_m[7];
			return sum;

		default:
			static int i;
			for( i = 0; i < _n; i++ )
				sum += _m[i]*mDet( _m, _n, 0, i );
			return sum;
	}
}

float mDet( float *_m, int _n )
{
	static float sum;
	sum = 0.f;
	switch( _n )
	{
		case 1:
			return _m[0];

		case 2:
			sum += _m[0]*_m[3];
			sum += -_m[1]*_m[2];
			return sum;

		case 3:
			sum += _m[0]*_m[4]*_m[8];
			sum += _m[1]*_m[5]*_m[6];
			sum += _m[2]*_m[3]*_m[7];
			sum += -_m[2]*_m[4]*_m[6];
			sum += -_m[1]*_m[3]*_m[8];
			sum += -_m[0]*_m[5]*_m[7];
			return sum;

		default:
			static int i;
			for( i = 0; i < _n; i++ )
				sum += _m[i]*mDet( _m, _n, 0, i );
			return sum;
	}
}

int mDet( int *const _m, int _n, int _r, int _c )
{
	static int row, col, size;
	static int i, j, temp;
	static int *matrix, *ptr;

	row = ( _r + 1 )%_n;
	col = ( _c + 1 )%_n;
	size = _n - 1;

	matrix = (int*)malloc( size*size*sizeof(int) );
	ptr = matrix;
	for( i = row; i < row + size; i++ )
		for( j = col; j < col + size; j++, ptr++ )
			ptr[0] = (int)powf( -1.f, (float)( i + j ) )*_m[ _n*( i%_n ) + j%_n ];

	temp = mDet( matrix, size );
	free( matrix );
	return temp;
}

float mDet( float *const _m, int _n, int _r, int _c )
{
	static int row, col, size;
	static int i, j;
	static float *matrix, *ptr, temp;

	row = ( _r + 1 )%_n;
	col = ( _c + 1 )%_n;
	size = _n - 1;

	matrix = (float*)malloc( size*size*sizeof(float) );
	ptr = matrix;
	for( i = row; i < row + size; i++ )
		for( j = col; j < col + size; j++, ptr++ )
			ptr[0] = powf( -1.f, (float)( i + j ) )*_m[ _n*( i%_n ) + j%_n ];

	temp = mDet( matrix, size );
	free( matrix );
	return temp;
}

void mTranspose( int *_s, int *_m, int _n )
{
	static int n, i;
	n = _n*_n;
	for( i = 0; i < n; i++, _m++ )
		if( i % ( _n + 1 ) > 0 )
			_m[0] = _s[ _n*( i % _n ) + ( i / _n ) ];
		else
			_m[0] = _s[i];
}

void mTranspose( float *_s, float *_m, int _n )
{
	static int n, i;
	n = _n*_n;
	for( i = 0; i < n; i++, _m++ )
		if( i % ( _n + 1 ) > 0 )
			_m[0] = _s[ _n*( i % _n ) + ( i / _n ) ];
		else
			_m[0] = _s[i];
}

bool mInverse( int *const _s, int *_m, int _n )
{
	static float det;
	det = (float)mDet( _s, _n );
	if( det == 0.f )
		return false;

	static int *matrix, *ptr;
	static int r, c;

	matrix = (int*)malloc( _n*_n*sizeof(int) );
	ptr = matrix;
	for( r = 0; r < _n; r++ )
		for( c = 0; c < _n; c++, ptr++ )
			ptr[0] = mDet( _s, _n, r, c );

	mTranspose( matrix, _m, _n );
	mMultiple( _m, 1.f/det, _n );
	free( matrix );
	return true;
}

bool mInverse( float *const _s, float *_m, int _n )
{
	static float det;
	det = mDet( _s, _n );
	if( det == 0.f )
		return false;

	static int r, c;
	static float *matrix, *ptr;
	matrix = (float*)malloc( _n*_n*sizeof(float) );
	ptr = matrix;
	for( r = 0; r < _n; r++ )
		for( c = 0; c < _n; c++, ptr++ )
			ptr[0] = mDet( _s, _n, r, c );

	mTranspose( matrix, _m, _n );
	mMultiple( _m, 1.f/det, _n );
	free( matrix );
	return true;
}


#endif