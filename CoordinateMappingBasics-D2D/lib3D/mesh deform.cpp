#include <windows.h>
#include <GL\gl.h>

#include "mesh deform.h"
#include "track.h"
#include "3dconfig.h"
#include "numerical\linear2.h"
#include "numerical\nonlinear.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stack>

using namespace MeshDeform;

//
// Macro
// 

#define SQR(x)	( (x) * (x) )
#define Uint8	unsigned char

//
// static global variables
//

extern const unsigned char *RainbowTexture;

namespace { 

	inline void DrawRainbowLine( void ){
		int i, viewport[4];
		float mp[16], mm[16];

		// save matrices
		glGetIntegerv( GL_MATRIX_MODE, &i );
		glGetIntegerv( GL_VIEWPORT, viewport );
		glGetFloatv( GL_PROJECTION_MATRIX, mp );
		glGetFloatv( GL_MODELVIEW_MATRIX, mm );
		// Set orthogonal projection
		glDisable( GL_DEPTH_TEST );
		glMatrixMode( GL_PROJECTION );
		glLoadIdentity();
		glMatrixMode( GL_MODELVIEW );
		glLoadIdentity();
		// Draw background image
		glRasterPos2f( -1.0f, -0.9f );
		glPixelZoom( 20.0f, ( float )( viewport[3] - viewport[1] ) * 0.9f * 0.00078215f );
		glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );
		glDrawPixels( 1, 256 * 5, GL_RGB, GL_UNSIGNED_BYTE, RainbowTexture );
		// Setting back
		glEnable( GL_DEPTH_TEST );
		glLoadMatrixf( mm );
		glMatrixMode( GL_PROJECTION );
		glLoadMatrixf( mp );
		glMatrixMode( i );
	}

	inline void SetColor( float x, float max ){
		x /= max;
		if( x < 0.2f ){				// blue to bluegreen
			glColor3f( 0.f, x * 5.0f, 1.f );
		}
		else if( x < 0.4f ){		// bluegreen to green
			x -= 0.2f;
			glColor3f( 0.f, 1.f, 1.f - x * 5.0f );
		}
		else if( x < 0.6f ){		// green to orange
			x -= 0.4f;
			glColor3f( x * 5.0f, 1.0f, 0.0f );
		}
		else if( x < 0.8f ){		// orange to red
			x -= 0.6f;
			glColor3f( 1.0f, 1.f - x * 5.0f, 0.0f );
		}
		else if( x < 1.0f ){		// red to white
			x -= 0.8f;
			glColor3f( 1.0f, x * 5.0f, x * 5.0f );
		}
		else{						// white
			glColor3f( 1.0f, 1.0f, 1.0f );
		}
	}

	inline void CompAxis( Vecf w, Parameters param, float axis[9], const int ns ){
		RT_Paramf rt = { 0.f, 0.f, 0.f, 0.f, 0.f, 0.f };
		Vector3f x, y, z;
		float *pw, &temp = x.x;
		int i;

		// initialize
		pw = w;
		rt.a = pw[0];
		rt.b = pw[1];
		rt.c = pw[2];
		x.y = x.z = y.x = y.z = z.x = z.y = 0.f;
		x.x = y.y = z.z = pw[3];
		TransformPointf( &x, rt );
		TransformPointf( &y, rt );
		TransformPointf( &z, rt );
		axis[0] = x.x;
		axis[1] = x.y;
		axis[2] = x.z;
		axis[3] = y.x;
		axis[4] = y.y;
		axis[5] = y.z;
		axis[6] = z.x;
		axis[7] = z.y;
		axis[8] = z.z;
		// add axices
		pw += 4;
		for( i = 0; i < ns; i++, pw += 4 ){
			rt.a = pw[0];
			rt.b = pw[1];
			rt.c = pw[2];
			x.y = x.z = y.x = y.z = z.x = z.y = 0.f;
			x.x = y.y = z.z = pw[3];
			TransformPointf( &x, rt );
			TransformPointf( &y, rt );
			TransformPointf( &z, rt );
			axis[0] += x.x * param[i];
			axis[1] += x.y * param[i];
			axis[2] += x.z * param[i];
			axis[3] += y.x * param[i];
			axis[4] += y.y * param[i];
			axis[5] += y.z * param[i];
			axis[6] += z.x * param[i];
			axis[7] += z.y * param[i];
			axis[8] += z.z * param[i];
		}
		// normalize
		temp = 1.f / sqrt( SQR( axis[0] ) + SQR( axis[1] ) + SQR( axis[2] ) );
		axis[0] *= temp;
		axis[1] *= temp;
		axis[2] *= temp;
		temp = 1.f / sqrt( SQR( axis[3] ) + SQR( axis[4] ) + SQR( axis[5] ) );
		axis[3] *= temp;
		axis[4] *= temp;
		axis[5] *= temp;
		temp = 1.f / sqrt( SQR( axis[6] ) + SQR( axis[7] ) + SQR( axis[8] ) );
		axis[6] *= temp;
		axis[7] *= temp;
		axis[8] *= temp;
	}

	namespace ML{										// for ML estimate
		__declspec( thread ) float *param;				// sementic parameters
		__declspec( thread ) int n;						// mesh number
		__declspec( thread ) int ns;					// sementic parameter dimension

		inline void LinearAxis( Vecf w, Mat2f axis ){
			float *pparam;
			int i;

			pparam = param;
			for( i = 0; i < n; i++, pparam += ns ){
				CompAxis( w, pparam, axis[i], ns );
			}
		}
	};
};

//
// Base Deformation Class
//

float BaseDeform::Residual( Mesh &m, Parameters p ){
	const int n = m.GetVertexSize();
	float e, temp, *p1, *p2;
	int i;

	// deform
	Mesh m2( m );
	Deform( m2, p );
	// residual
	e = 0.f;
	p1 = m.GetVertex(0);
	p2 = m2.GetVertex(0);
	for( i = 0; i < n; i++ ){
		temp = ( *p1 ) - ( *p2 );	// x
		e += SQR( temp );
		p1++;	p2++;
		temp = ( *p1 ) - ( *p2 );	// y
		e += SQR( temp );
		p1++;	p2++;
		temp = ( *p1 ) - ( *p2 );	// z
		e += SQR( temp );
		p1++;	p2++;
	}

	// normalize
	temp = 0.3333333f / ( float )n;
	return e * temp;
}

bool BaseDeform::Import( wchar_t path[] ){
	bool ret;
	FILE *infile;
	_wfopen_s( &infile, path, L"rb" );
	if( infile == NULL ) return false;
	ret = Import( infile );
	fclose( infile );
	return ret;
}

bool BaseDeform::Export( wchar_t path[] ){
	bool ret;
	FILE *outfile;
	_wfopen_s( &outfile, path, L"wb" );
	if( outfile == NULL ) return false;
	ret = Export( outfile );
	fclose( outfile );
	return ret;
}

//
// linear deformation
//

LinearDeform::LinearDeform( size_t vertex, size_t dim ):
_dim( dim )
{
	size_t buf;
	_nv = vertex;

	if( ( buf = 3 * vertex * ( dim + 1 ) ) ){
		_coef = ( float * )malloc( buf * sizeof( float ) );
	}
	else{
		_coef = NULL;
	}
}

LinearDeform::~LinearDeform( void ){
	Clear();
}

LinearDeform & LinearDeform::operator =( const LinearDeform & src ){
	size_t buf;
	_nv = src._nv;
	_dim = src._dim;
	
	if( ( buf = 3 * _nv * ( _dim + 1 ) ) ){
		_coef = ( float * )realloc( _coef, buf * sizeof( float ) );
	}
	else{
		if( _coef ) {
			free( _coef );
			_coef = NULL;
		}
	}
	return ( *this );
}

// Modifiers

void LinearDeform::Clear( void ){
	this->Memory::Clear();
	if( _coef ) {
		free( _coef );
		_coef = NULL;
	}
	_nv = _dim = 0;
}

// Functions

int LinearDeform::Training( vector< Mesh* > &mv, Parameters param ){
	const size_t n = mv.size(), d = _dim + 1, v3 = 3 * _nv;
	Vecf y, pinv;
	Matf x, inv, pparam, px, pc;
	size_t i, j, k;

	// check
	if( n == 0 || mv[0]->GetVertexSize() != _nv ){
		return 0;
	}

	// memory pool
	x = ( Matf )Allocate( ( 2 * n * d + n ) * sizeof( float ) );	// n x ( dim + 1 )
	y = x + n * d;													// n x 1	
	inv = y + n;													// ( dim + 1 ) x n
	if( x == NULL ) return 0;

	// build constant x
	px = x;
	pparam = param;
	for( i = 0; i < n; i++, px += d, pparam += _dim ){
		memcpy( px, pparam, _dim * sizeof( float ) );	// copy from param
		px[ _dim ] = 1.0;								// 1 for constant coefficient
	}
	// build psudo-inverse of x
	i = pinvf( inv, x, n, d, 0.0000000001f );
	if( i == 0 ) return 0;

	// linear regression
	memset( _coef, 0x00, v3 * d * sizeof( float ) );
	pc = _coef;
	for( i = 0; i < _nv; i++ ){
		// build y for p[i].x
		pinv = inv;
		for( j = 0; j < d; j++, pc++ ){
			for( k = 0; k < n; k++, pinv++ ){
				pc[0] += pinv[0] * mv[k]->GetVertex(i)[0];
			}
		}
		// build y for p[i].y
		pinv = inv;
		for( j = 0; j < d; j++, pc++ ){
			for( k = 0; k < n; k++, pinv++ ){
				pc[0] += pinv[0] * mv[k]->GetVertex(i)[1];
			}
		}
		// build y for p[i].z
		pinv = inv;
		for( j = 0; j < d; j++, pc++ ){
			for( k = 0; k < n; k++, pinv++ ){
				pc[0] += pinv[0] * mv[k]->GetVertex(i)[2];
			}
		}
	}
	return 1;
}

int LinearDeform::Training( vector< Mesh > &m, vector< Parameters > &p ){
	if( m.size() != p.size() ) return 0;

	const int n = m.size();
	vector< Mesh * > mesh( n );
	float *param = ( float * )malloc( n * _dim * sizeof( float ) );
	int i;
	
	for( i = 0; i < n; i++ ){
		mesh[i] = &( m[i] );
		memcpy( param + i * _dim, p[i], _dim * sizeof( float ) );
	}
	i = Training( mesh, param );
	free( param );
	return i;
}

int LinearDeform::Deform( Mesh &out, Parameters param ){	
	const size_t v3 = 3 * _nv, d = _dim + 1;
	float *pc, *vertex_array = out.GetVertex(0);
	size_t i, j;

	if( out.GetVertexSize() != _nv ) return 0;

	// coefficient * param
	pc = _coef;
	for( i = 0; i < v3; i++ ){
		// const cooeficient
		vertex_array[i] = pc[_dim];
		for( j = 0; j < _dim; j++, pc++ ){
			vertex_array[i] += pc[0] * param[j];
		}
		pc++;
	}
	return 1;
}

// file operation

bool LinearDeform::Import( FILE *infile ){
	fread( &_nv, sizeof( size_t ), 1, infile );
	fread( &_dim, sizeof( size_t ), 1, infile );
	_coef = ( float * )malloc( 3 * _nv * ( _dim + 1 ) * sizeof( float ) );
	fread( _coef, sizeof( float ), 3 * _nv * ( _dim + 1 ), infile );
	return true;
}

bool LinearDeform::Export( FILE *outfile ){
	fwrite( &_nv, sizeof( size_t ), 1, outfile );
	fwrite( &_dim, sizeof( size_t ), 1, outfile );
	fwrite( _coef, sizeof( float ), 3 * _nv * ( _dim + 1 ), outfile );
	return true;
}

// 
// Linear Deformation with Constant Geometry
//

ConstGeometryDeform::ConstGeometryDeform( size_t vertex, size_t dim ): 
LinearDeform( vertex, dim ){		
	CGD_Group g0 = { vector< int >(0), NULL, NULL };
	
	// linear group
	for( size_t i = 0; i < vertex; i++ ){
		g0.vertex.push_back( i );
	}
	_group.push_back( g0 );				

	if( vertex ){
		_v = ( Point3f * )malloc( 3 * vertex * sizeof( float ) );	// mesh verteces
	}
	else{
		_v = NULL;
	}
}

ConstGeometryDeform::~ConstGeometryDeform( void ){
	Clear();
}

ConstGeometryDeform & ConstGeometryDeform::operator =( const ConstGeometryDeform &src ){
	size_t size;
	this->LinearDeform::operator =( src );						// copy linear deform

	// copy verteces
	size = 3 * _nv * sizeof( float );
	_v = ( Point3f * )realloc( _v, size );
	memcpy( _v, src._v, size );

	// copy groups
	list< CGD_Group >::iterator it;
	float *buf;
	_group = src._group;
	size = ( 4 * ( _dim + 1 ) + 3 * ( _dim + 2 ) ) * sizeof( float );
	for( it = ++_group.begin(); it != _group.end(); it++ ){
		buf = ( float * )malloc( size );	// copy ra
		memcpy( buf, it->w, size );
		it->w = buf;						// 4 * ( d + 1 )
		it->coef = buf + 4 * ( _dim + 1 );	// 3 * ( d + 2 )
	}
	return (*this);
}

// protected

void ConstGeometryDeform::CheckGroup( void ){
	list< CGD_Group >::iterator it;

	for( it = ++_group.begin(); it != _group.end(); ){
		if( it->vertex.size() < 3 ) {
			// delete this group
			free( it->w );
			it = _group.erase( it );
		}
		else{
			it++;
		}
	}
}

// Modifiers

int ConstGeometryDeform::SetMesh( Mesh &m ){
	if( _nv != m.GetVertexSize() ) return 0;					// vertex count error
	if( _v == NULL ){
		_v = ( Point3f* )malloc( _nv * sizeof( Point3f ) );		// memory allocate
	}
	memcpy( _v, m.GetVertex(0), _nv * sizeof( Point3f ) );		// copy all verteces
	return _nv;
}

void ConstGeometryDeform::PushGroup( vector< int > &g ){
	CGD_Group gi = { g, ( float * )calloc( 4 * ( _dim + 1 ) + 3 * ( _dim + 2 ), sizeof( float ) ) };
	gi.coef = gi.w + 4 * ( _dim + 1 );
	_group.push_back( gi );
}

int ConstGeometryDeform::PopGroup( vector< int > &g ){
	list< CGD_Group >::iterator it = _group.begin();
	for( it++; it != _group.end(); it++ ){
		if( it->vertex == g ) {
			free( it->w );
			_group.erase( it );
			return 1;
		}
	}
	return 0;
}

void ConstGeometryDeform::Clear(){
	this->LinearDeform::Clear();
	if( _v ){
		free( _v );
		_v = NULL;
	}
	// free group linkedlist
	list< CGD_Group >::iterator it = _group.begin();
	for( it++; it != _group.end(); it++ ){
		free( it->w );
	}
	_group.clear();
}

GroupList ConstGeometryDeform::GetGroupList( void ){
	list< CGD_Group >::iterator it;
	GroupList ret(0);

	for( it = ++_group.begin(); it != _group.end(); it++ ){
		ret.push_back( it->vertex );
	}
	return ret;
}

// Functions

int ConstGeometryDeform::Training( vector< Mesh* > &mv, Parameters param ){	
	// Check group
	CheckGroup();

	list< CGD_Group >::iterator it;
	vector< int > &g0 = _group.begin()->vertex;
	const size_t n = mv.size(), d = _dim + 1, d2 = _dim + 2, v3 = 3 * _nv;
	const size_t gn = _group.size();
	Mat2f p;													// projected linear vertecex in each mesh
	Vector3f x, y, z, *c;
	float buf, *p1, *p2, *p3, *ptr, *w, *inv, *mat, *pparam;
	size_t i, j, nv;

	// check
	if( _coef == NULL || _group.size() < 1 ) return 0;
	for( it = ++_group.begin(); it != _group.end(); it++ ){
		w = it->w;
		if( it->w == NULL ) return 0;
	}

	// initialize
	SetMesh( *mv[0] );
	w = _group.begin()->w;

	// linear deformation
	this->LinearDeform::Training( mv, param );					// _Training_Linear( mv, param );
	
	// memory pool for rotation axices in each mesh
	p = ( Mat2f )Allocate( n * sizeof( Vecf ) + n * ( 9 + d2 + d2 ) * sizeof( float ) + n * sizeof( Point3f ) );
	p[0] = ( Vecf )( p + n );							// 2d rotation angle matrix
	for( i = 1; i < n; i++ ) p[i] = p[ i - 1 ] + 9;		// rotation angles
	inv = p[ i - 1 ] + 9;								// psudo-inverse matrix
	mat = inv + n * d2;									// linear matrix
	c = ( Point3f* )( mat + n * d2 );					// center array

	// MLE algorithm variables
	ML::param = param;
	ML::ns = _dim;
	ML::n = n;

	// solve for each group
	for( it = ++_group.begin(); it != _group.end(); it++ ){
		memset( c, 0x00, n * sizeof( Point3f ) );
		// axis setting for each sample
		for( i = 0; i < n; i++ ){	// for each sample
			p1 = mv[i]->GetVertex( it->vertex[0] );
			p2 = mv[i]->GetVertex( it->vertex[1] );
			p3 = mv[i]->GetVertex( it->vertex[2] );
			// x
			x.x = p2[0] - p1[0];
			x.y = p2[1] - p1[1];
			x.z = p2[2] - p1[2];
			buf = 1.f / sqrt( SQR( x.x ) + SQR( x.y ) + SQR( x.z ) );
			x.x *= buf;	x.y *= buf;	x.z *= buf;
			// y
			y.x = p3[0] - p1[0];
			y.y = p3[1] - p1[1];
			y.z = p3[2] - p1[2];
			buf = Dotf( x, y );
			y.x -= buf * x.x;
			y.y -= buf * x.y;
			y.z -= buf * x.z;
			buf = 1.f / sqrt( SQR( y.x ) + SQR( y.y ) + SQR( y.z ) );
			y.x *= buf;	y.y *= buf;	y.z *= buf;
			// z
			Crossf( x, y, &z );
			// project for each reference points
			ptr = p[i];											// i-th mesh
			memcpy( ptr + 0, &x, 3 * sizeof( float ) );
			memcpy( ptr + 3, &y, 3 * sizeof( float ) );
			memcpy( ptr + 6, &z, 3 * sizeof( float ) );

			// compute center
			nv = it->vertex.size();
			for( j = 0; j < nv; j++ ){
				p1 = mv[i]->GetVertex( it->vertex[j] );
				c[i].x += p1[0];	c[i].y += p1[1];	c[i].z += p1[2];
			}
			buf = 1.f / ( float )nv;
			c[i].x *= buf;	c[i].y *= buf;	c[i].z *= buf;
		}

		// apply LevMar algorithm	
		w = it->w;
		memset( w, 0x00, 4 * d * sizeof( float ) );	
		for( i = 0; i < d; i++ ){
			w[ i * 4 + 3 ] =  1.f;
		}
		buf = LevMar2f( ML::LinearAxis, w, p, 4 * d, n, 9, 30, 2.f, 
			0.000000001f, 0.000000001f, 0.000000001f );

		// zerolize linear coefficient
		memset( it->coef, 0x00, 3 * d2 * sizeof( float ) );

		// linear regression
		w = it->coef;
		for( int k = 0; k < 3; k++ ){								// for x, y, z-value
			ptr = mat;
			pparam = param;
			for( i = 0; i < n; i++, ptr += /*4*/2, pparam += _dim ){		// build matrix
				memcpy( ptr, pparam, _dim * sizeof( float ) ); 
				ptr += _dim;
				// x.k, y.k, z.k value of i-th mesh
				ptr[0] = p[i][ 0 + k ] + p[i][ 3 + k ] + p[i][ 6 + k ];
				ptr[1] = 1.f;
			}
			// compute variance of x.k, y.k, z.k value
//			for( j = 0; j < 3; j++ ){
				ptr = mat + _dim;
				memset( &x, 0x00, sizeof( Point3f ) );
				for( i = 0; i < n; i++, ptr += d2 ){
					x.x += SQR( ptr[0] );							// sun(X^2)
					x.y += ptr[0];									// sum(X)
				}
				x.z = 1.0f / ( float )n;							// 1/n
				x.y *= x.z;											// E(X)
				buf = x.x * x.z - SQR( x.y );						// E[X^2] - [ E(X) ]^2
				if( buf < 0.01f ){									// if this threatment is constant
					// set to zero
					ptr = mat + _dim;
					for( i = 0; i < n; i++, ptr += /*d4*/d2 ) ptr[0] = 0.0f;
				}
//			}
 			pinvf( inv, mat, n, d2, 0.000001f );					// compute psudo-inverse matrix
			ptr = inv;
			for( i = 0; i < d2; i++, w++ ){
				for( j = 0; j < n; j++, ptr++ ){
					switch( k ){
						case 0:	w[0] += ptr[0] * c[j].x; break;
						case 1:	w[0] += ptr[0] * c[j].y; break;
						case 2:	w[0] += ptr[0] * c[j].z; break;
					}
				}
			}
		}
	}
#ifdef _DEBUG
	w = (++_group.begin())->w;	// for debug
#endif

	return gn;
}

int ConstGeometryDeform::Deform( Mesh &mesh, Parameters param ){
	list< CGD_Group >::iterator it;
	Vector3f x, y, z, p1, p2, p3;
	float m[3][3][4], buf[9], &temp = buf[0], *px, *py, *pz;
	size_t i, n;

	// linear deform
	if( LinearDeform::Deform( mesh, param ) == 0 ) return 0;
	
	// constant deform
	for( it = ++_group.begin(); it != _group.end(); it++ ){	// algorithm 4
		// compute old rotation
		p1 = _v[ it->vertex[0] ];
		p2 = _v[ it->vertex[1] ];
		p3 = _v[ it->vertex[2] ];
		// x
		x.x = p2.x - p1.x;
		x.y = p2.y - p1.y;
		x.z = p2.z - p1.z;
		temp = 1.f / sqrt( SQR( x.x ) + SQR( x.y ) + SQR( x.z ) );
		x.x *= temp;	x.y *= temp;	x.z *= temp;
		// y
		y.x = p3.x - p1.x;
		y.y = p3.y - p1.y;
		y.z = p3.z - p1.z;
		temp = Dotf( x, y );
		y.x -= temp * x.x;
		y.y -= temp * x.y;
		y.z -= temp * x.z;
		temp = 1.f / sqrt( SQR( y.x ) + SQR( y.y ) + SQR( y.z ) );
		y.x *= temp;	y.y *= temp;	y.z *= temp;
		// z
		Crossf( x, y, &z );
		// rotation matrix
		m[0][0][0] = x.x;	m[0][0][1] = y.x;	m[0][0][2] = z.x;
		m[0][1][0] = x.y;	m[0][1][1] = y.y;	m[0][1][2] = z.y;
		m[0][2][0] = x.z;	m[0][2][1] = y.z;	m[0][2][2] = z.z;
		
		// compute old center
		n = it->vertex.size();
		temp = 1.f / ( float )n;
		m[0][0][3] = m[0][1][3] = m[0][2][3] = 0.f;
		for( i = 0; i < n; i++ ){
			p1 = _v[ it->vertex[i] ];
			m[0][0][3] += p1.x;
			m[0][1][3] += p1.y;
			m[0][2][3] += p1.z;
		}
		m[0][0][3] *= temp;	m[0][1][3] *= temp;	m[0][2][3] *= temp;

		// compute new rotation
		CompAxis( it->w, param, buf, _dim );
		m[1][0][0] = buf[0];	m[1][0][1] = buf[3];	m[1][0][2] = buf[6];
		m[1][1][0] = buf[1];	m[1][1][1] = buf[4];	m[1][1][2] = buf[7];
		m[1][2][0] = buf[2];	m[1][2][1] = buf[5];	m[1][2][2] = buf[8];

		// compute new center
		px = it->coef;
		py = px + _dim + 2;//4;
		pz = py + _dim + 2;//4;
		m[1][0][3] = m[1][1][3] = m[1][2][3] = 0.f;
		for( i = 0; i < _dim; i++, px++, py++, pz++ ){	// linear regression for x, y, z, parameter coefficient
			m[1][0][3] += param[i] * px[0];
			m[1][1][3] += param[i] * py[0];
			m[1][2][3] += param[i] * pz[0];
		}
		// linear regression for rotation axis coefficient
		m[1][0][3] += buf[0] * px[0] + buf[3] * px[0] + buf[6] * px[0] + px[1];
		m[1][1][3] += buf[1] * py[0] + buf[4] * py[0] + buf[7] * py[0] + py[1];
		m[1][2][3] += buf[2] * pz[0] + buf[5] * pz[0] + buf[8] * pz[0] + pz[1];
		/* average translation 
		m[1][0][3] = m[1][1][3] = m[1][2][3] = 0.f;
		n = it->vertex.size();
		temp = 1.f / ( float )n;
		for( i = 0; i < n; i++ ){
			px = mesh.GetVertex( it->vertex[i] );
			m[1][0][3] += px[0];
			m[1][1][3] += px[1];
			m[1][2][3] += px[2];
		}
		m[1][0][3] *= temp;	m[1][1][3] *= temp;	m[1][2][3] *= temp;
		/**/

		// inverse old matrix
		InvRTMatf( m[2], m[0] );
		// new * inv(old)
		MultiMatf( m[0], m[1], m[2] );
		
		// do transforming
		n = it->vertex.size();
		for( i = 0; i < n; i++ ){
			x = _v[ it->vertex[i] ];
			MultiPointf( &x, m[0] );
			memcpy( mesh.GetVertex( it->vertex[i] ), &x, 3 * sizeof( float ) );
		}
	}

	return 1;
}
/*
float ConstGeometryDeform::Residual( Mesh &m, Parameters p ){
	const int n = m.GetVertexSize();
	float e, temp, *p1, *p2, rt[3][4];
	Point3fv_Pair pair;
	Transform_Paramf transform_param = { 10, 0.00001f, 0.01f, 0.5f };
	bool *flag;
	list< CGD_Group >::iterator it;
	int i;

	// deform
	Mesh m2( m );
	LinearDeform::
		Deform( m2, p );
	e = 0.f;

	// build linear group
	flag = ( bool* )calloc( n, sizeof( bool ) );	// flag array
	for( it = ++_group.begin(); it != _group.end(); it++ ){
		// const geometry pair
		pair.n = it->vertex.size();
		pair.pt1 = ( Point3fv* )Allocate( ( pair.n + pair.n ) * sizeof( Point3fv ) );
		pair.pt2 = pair.pt1 + pair.n;
		// record
		for( i = 0; i < pair.n; i++ ){
			flag[ it->vertex[i] ] = true;
			pair.pt1[i].ptr = m.GetVertex( it->vertex[i] );
			pair.pt2[i].ptr = m2.GetVertex( it->vertex[i] );
		}
		// solving rigid transform
		RigidTransform_SVDMatfv( pair, rt );
		temp = RigidTransformMatfv( pair, rt, transform_param );
		// residual
		e += temp;
	}

	// residual
	p1 = m.GetVertex(0);
	p2 = m2.GetVertex(0);
	for( i = 0; i < n; i++ ){
		if( flag[i] ) continue;
		temp = ( *p1 ) - ( *p2 );	// x
		e += SQR( temp );
		p1++;	p2++;
		temp = ( *p1 ) - ( *p2 );	// y
		e += SQR( temp );
		p1++;	p2++;
		temp = ( *p1 ) - ( *p2 );	// z
		e += SQR( temp );
		p1++;	p2++;
	}

	// normalize
	free( flag );
//	temp = 0.333333f / ( float )n;
//	return e * temp;
	return e;
}/**/

// file operation

bool ConstGeometryDeform::Import( FILE *infile ){
	size_t i, j, n, temp;
	CGD_Group g;
	if( !this->LinearDeform::Import( infile ) ) return false;
	
	// mesh
	_v = ( Point3f * )malloc( _nv * sizeof( Point3f ) );
	fread( _v, sizeof( Point3f ), _nv, infile );

	// group count
	fread( &n, sizeof( size_t ), 1, infile );
	if( !_group.empty() ) _group.clear();		// clear elements

	// each group
	for( i = 0; i < n; i++ ){
		// vertex index count
		fread( &temp, sizeof( size_t ), 1, infile );
		g.vertex.resize( temp );
		// vertex index array
		for( j = 0; j < temp; j++ ){
			fread( &( g.vertex[j] ), sizeof( int ), 1, infile );
		}
		if( i == 0 ) {
			g.w = g.coef = NULL;
			_group.push_back( g );
			continue;
		}
		// roate angle
		g.w = ( float * )malloc( 4 * ( _dim + 1 ) * sizeof( float ) );
		fread( g.w, sizeof( float ), 4 * ( _dim + 1 ), infile );
		// translate regression
		g.coef = ( float * )malloc( 3 * ( _dim + 2 ) * sizeof( float ) );
		fread( g.coef, sizeof( float ), 3 * ( _dim + 2 ), infile );

		// add to list
		_group.push_back( g );
	}
	return true;
}

bool ConstGeometryDeform::Export( FILE *outfile ){
	size_t i = 0, n = 0;
	list< CGD_Group >::iterator it;
	if( !this->LinearDeform::Export( outfile ) ) return false;
	
	// mesh
	fwrite( _v, sizeof( Point3f ), _nv, outfile );
	
	// count group
	for( it = _group.begin(); it != _group.end(); it++, n++ ){}
	fwrite( &n, sizeof( size_t ), 1, outfile );

	// each group
	for( it = _group.begin(); it != _group.end(); it++ ){
		// vertex index count
		n = it->vertex.size();
		fwrite( &n, sizeof( size_t ), 1, outfile );
		// vertex index array
		for( i = 0; i < n; i++ ){
			fwrite( &( it->vertex[i] ), sizeof( int ), 1, outfile );
		}
		// rotate angle
		if( it == _group.begin() ) continue;
		fwrite( it->w, sizeof( float ), 4 * ( _dim + 1 ), outfile );
		// translate regression
		fwrite( it->coef, sizeof( float ), 3 * ( _dim + 2 ), outfile );
	}
	return true;
}

//
// Global Functions
//

float MeshDeform::ShapeSSE( BaseDeform &deform, GroupList &group, Mesh &m, Parameters p ){
	const int n = m.GetVertexSize();
	float e, temp, rt[3][4];
	Point3fv_Pair pair;
	Transform_Paramf transform_param = { 100, 0.00001f, 0.01f, 0.5f };
	GroupList::iterator it;
	int i;

	// deform
	Mesh m2( m );
	deform.Deform( m2, p );
	e = 0.f;

	if( m2.GetVertexSize() != n ) return FLT_MAX;

	// compute for each group
	for( it = group.begin(); it != group.end(); it++ ){
		// const geometry pair
		pair.n = it->size();
		pair.pt1 = ( Point3fv* )deform.Allocate( ( pair.n + pair.n ) * sizeof( Point3fv ) );
		pair.pt2 = pair.pt1 + pair.n;
		// record
		for( i = 0; i < pair.n; i++ ){
			pair.pt1[i].ptr = m.GetVertex( (*it)[i] );
			pair.pt2[i].ptr = m2.GetVertex( (*it)[i] );
		}
		// solving rigid transform
		RigidTransform_SVDMatfv( pair, rt );
		temp = RigidTransformMatfv( pair, rt, transform_param );
		// residual
		e += temp;
	}
	// normalize
	return e;
}

float MeshDeform::ShapeMSE( BaseDeform &deform, GroupList &group, Mesh &m, Parameters p ){
	GroupList::iterator it;
	float ret;
	int n = 0;

	for( it = group.begin(); it != group.end(); it++ ){
		n += it->size();
	}
	ret = ShapeSSE( deform, group, m, p );
	return ( ret == FLT_MAX ) ? FLT_MAX : ( ret * 0.333f / ( float )n );
}

float MeshDeform::DrawError( BaseDeform &deform, Mesh &m, Parameters p, float max ){
	const int n = m.GetVertexSize(), fn = m.GetTriangleSize();
	float *e, e2, temp, *p1, *p2, ret = 0.0f;
	int i, j;

	// deform
	Mesh m2( m );
	deform.Deform( m2, p );

	if( m2.GetVertexSize() != n ) return FLT_MAX;
	p1 = m.GetVertex(0);
	p2 = m2.GetVertex(0);

	// memory allocate
	e = ( float * )calloc( n, sizeof( float ) );

	// compute error for each vertex
	for( i = 0; i < n; i++, p1 += 3, p2 += 3 ){
		temp = p1[0] - p2[0];
		e2 = SQR( temp );
		temp = p1[1] - p2[1];
		e2 += SQR( temp );
		temp = p1[2] - p2[2];
		e2 += SQR( temp );
		e[i] = sqrt( e2 );
		ret += e2;
	}

	// rainbow texture ( 256 * 5 )
	unsigned char tex2d, clr_mat, light;
	GLuint tex[1];

	glGetBooleanv( GL_TEXTURE_2D, &tex2d );
	glEnable( GL_TEXTURE_2D );
	glGenTextures( 1, tex );
	glBindTexture( GL_TEXTURE_2D, tex[0] );
	glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D( GL_TEXTURE_2D, 0, GL_RGB, 1280, 1,
		0, GL_RGB, GL_UNSIGNED_BYTE, RainbowTexture );

	// drawing
	glGetBooleanv( GL_COLOR_MATERIAL, &clr_mat );
	glGetBooleanv( GL_LIGHTING, &light );
	glDisable( GL_COLOR_MATERIAL );
	glDisable( GL_LIGHTING );
	glBegin( GL_TRIANGLES );
	glColor4f( 1.0f, 1.0f, 1.0f, 1.0f );
	for( i = 0; i < fn; i++ ){
		j = m2.GetTriangle(i)->vertex_index[0];
		temp = 0.00078125f + e[j] / max;
		glTexCoord2f( ( temp > 0.99921875f ) ? 0.99921875f : temp, 0.f );
		glVertex3fv( m2.GetVertex( j ) );

		j = m2.GetTriangle(i)->vertex_index[1];
		temp = 0.00078125f + e[j] / max;
		glTexCoord2f( ( temp> 0.99921875f ) ? 0.99921875f : temp, 0.f );
		glVertex3fv( m2.GetVertex( j ) );

		j = m2.GetTriangle(i)->vertex_index[2];
		temp = 0.00078125f + e[j] / max;
		glTexCoord2f( ( temp> 0.99921875f ) ? 0.99921875f : temp, 0.f );
		glVertex3fv( m2.GetVertex( j ) );
	}
	glEnd();
	if( tex2d == GL_FALSE ) glDisable( GL_TEXTURE_2D );

	// wire frame
	glEnable( GL_BLEND );
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
	glColor4f( 1.0f, 1.0f, 1.0f, .2f );
	m2.Draw( MESH_NORMAL | MESH_WIRE );
	glDisable( GL_BLEND );

	if( clr_mat == GL_TRUE ) glEnable( GL_COLOR_MATERIAL );
	if( light == GL_TRUE ) glEnable( GL_LIGHTING );

	DrawRainbowLine();

	free( e );
	return ret;
}

float MeshDeform::DrawShapeError( BaseDeform &deform, GroupList &group, Mesh &m, Parameters p, float max ){
	const int n = m.GetVertexSize(), fn = m.GetTriangleSize();
	float *e, e2, temp, rt[3][4], pt[3], ret;
	Point3fv_Pair pair;
	Transform_Paramf transform_param = { 100, 0.00001f, 0.01f, 0.5f };
	GroupList::iterator it;
	int i, j;

	// deform
	Mesh m2( m );
	deform.Deform( m2, p );

	if( m2.GetVertexSize() != n ) return FLT_MAX;

	// memory allocate
	e = ( float * )calloc( n, sizeof( float ) );
	ret = 0.0f;
	// compute for each group
	for( it = group.begin(); it != group.end(); it++ ){
		// const geometry pair
		pair.n = it->size();
		pair.pt1 = ( Point3fv* )deform.Allocate( ( pair.n + pair.n ) * sizeof( Point3fv ) );
		pair.pt2 = pair.pt1 + pair.n;
		// record
		for( i = 0; i < pair.n; i++ ){
			pair.pt1[i].ptr = m.GetVertex( (*it)[i] );
			pair.pt2[i].ptr = m2.GetVertex( (*it)[i] );
		}
		// solving rigid transform
		RigidTransform_SVDMatfv( pair, rt );
		ret += RigidTransformMatfv( pair, rt, transform_param );
		// residual
		for( i = 0; i < pair.n; i++ ){
			memcpy( pt, pair.pt1[i].ptr, 3 * sizeof( float ) );
			MultiPointfv( pair.pt1 + i, rt );
			temp = pair.pt1[i].ptr[0] - pair.pt2[i].ptr[0];
			e2 = SQR( temp );
			temp = pair.pt1[i].ptr[1] - pair.pt2[i].ptr[1];
			e2 += SQR( temp );
			temp = pair.pt1[i].ptr[2] - pair.pt2[i].ptr[2];
			e2 += SQR( temp );
			e[ (*it)[i] ] += sqrt( e2 );
			memcpy( pair.pt1[i].ptr, pt, 3 * sizeof( float ) );
		}
	}

	// rainbow texture ( 256 * 5 )
	unsigned char tex2d, clr_mat, light;
	GLuint tex[1];

	glGetBooleanv( GL_TEXTURE_2D, &tex2d );
	glEnable( GL_TEXTURE_2D );
	glGenTextures( 1, tex );
	glBindTexture( GL_TEXTURE_2D, tex[0] );
	glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexImage2D( GL_TEXTURE_2D, 0, GL_RGB, 1280, 1,
		0, GL_RGB, GL_UNSIGNED_BYTE, RainbowTexture );

	// drawing
	glGetBooleanv( GL_COLOR_MATERIAL, &clr_mat );
	glGetBooleanv( GL_LIGHTING, &light );
	glDisable( GL_COLOR_MATERIAL );
	glDisable( GL_LIGHTING );
	glBegin( GL_TRIANGLES );
	glColor4f( 1.0f, 1.0f, 1.0f, 1.0f );
	for( i = 0; i < fn; i++ ){
		j = m2.GetTriangle(i)->vertex_index[0];
		temp = 0.00078125f + e[j] / max;
		glTexCoord2f( ( temp > 0.99921875f ) ? 0.99921875f : temp, 0.f );
		glVertex3fv( m2.GetVertex( j ) );

		j = m2.GetTriangle(i)->vertex_index[1];
		temp = 0.00078125f + e[j] / max;
		glTexCoord2f( ( temp> 0.99921875f ) ? 0.99921875f : temp, 0.f );
		glVertex3fv( m2.GetVertex( j ) );

		j = m2.GetTriangle(i)->vertex_index[2];
		temp = 0.00078125f + e[j] / max;
		glTexCoord2f( ( temp> 0.99921875f ) ? 0.99921875f : temp, 0.f );
		glVertex3fv( m2.GetVertex( j ) );
	}
	glEnd();
	if( tex2d == GL_FALSE ) glDisable( GL_TEXTURE_2D );

	// wire frame
	glEnable( GL_BLEND );
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
	glColor4f( 1.0f, 1.0f, 1.0f, .2f );
	m2.Draw( MESH_NORMAL | MESH_WIRE );
	glDisable( GL_BLEND );

	if( clr_mat == GL_TRUE ) glEnable( GL_COLOR_MATERIAL );
	if( light == GL_TRUE ) glEnable( GL_LIGHTING );

	DrawRainbowLine();

	free( e );
	return ret;
}