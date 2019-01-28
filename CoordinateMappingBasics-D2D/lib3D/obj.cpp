
//#pragma comment( lib, "opengl32.lib" )
//#pragma comment( lib, "glu32.lib" )
//#pragma comment( lib, "glaux.lib" )

#include <string.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <time.h>
//#include <image\bitmap.h>
//#include <image\imgproc.h>
#include "image\image.h"
#include "numerical\algo.h"
#include "numerical\linear2.h"
#include "obj.h"
#include <GL\gl.h>

#include <include\version.h>

// global variables
static	float gv[16];

//
// macro
//

#define Cross( a1, a2, a3, b1, b2, b3, c1, c2, c3 ) \
	(c1) = (a2) * (b3) - (a3) * (b2); \
	(c2) = (a3) * (b1) - (a1) * (b3); \
	(c3) = (a1) * (b2) - (a2) * (b1);

#define Dot( a, b ) \
	( (a)[0] * (b)[0] + (a)[1] * (b)[1] + (a)[2] * (b)[2] )

#define SQR( x ) ( (x) * (x) )

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#if _VC6_
#define fprintf_s fprintf
#endif

//
// inline function
//

inline void LoadIdentity( float m[16] ){
	int i;
	for( i = 0; i < 16; i++ ) m[i] = 0.0;
	for( i = 0; i < 16; i += 5 ) m[i] = 1.0;
}

inline void Translate( float m[16], float x, float y, float z ){
	::LoadIdentity( m );
	m[3] = x;
	m[7] = y;
	m[11] = z;
}

inline void MultiMatrix( float a[16], float b[16], float c[16] ){
	int i, j, k;
	for( i = 0; i < 4; i++ ){
		for( j = 0; j < 4; j++ ){
			c[i * 4 + j] = 0.0;
			for( k = 0; k < 4; k++ ) c[i * 4 + j] += a[i * 4 + k] * b[k * 4 + j];
		}
	}
}

inline void MultiVector( float a[16], float b[4] ){
	int i, j;

	for( i = 0; i < 4; i++ ) gv[i] = b[i];
	for( i = 0; i < 4; i++ ){
		b[i] = 0.0;
		for( j = 0; j < 4; j++ ) b[i] += a[i * 4 + j] * gv[j];
	}
}

inline void Rotate( float m[16], float deg, float x, float y, float z ){
	float theta = deg / 180.0f * ( float )M_PI,
		uvw[16] = {
			0.0f, 0.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f },
		r[16] = {
			cos( theta ), -sin( theta ), 0.0f, 0.0f,
			sin( theta ), cos( theta ), 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f };
	;

	if( y == 0.0f && z == 0.0f ) uvw[1] = uvw[6] = uvw[8] = 1.0f;
	else{
		// normalize Z-axis
		gv[0] = sqrt( SQR(x) + SQR(y) + SQR(z) );
		uvw[8] = x / gv[0];
		uvw[9] = y / gv[0];
		uvw[10] = z / gv[0];

		// generate X-axis, X * Z = Zx
		x = 1.0f - uvw[8] * uvw[8];	// x = e1 - ( e1 * z ) * z
		y = -uvw[9] * uvw[8];
		z = -uvw[10] * uvw[8];
		gv[0] = sqrt( SQR(x) + SQR(y) + SQR(z) );
		uvw[0] = x / gv[0];
		uvw[1] = y / gv[0];
		uvw[2] = z / gv[0];

		// y = z cross x
		Cross( uvw[8], uvw[9], uvw[10], uvw[0], uvw[1], uvw[2], uvw[4], uvw[5], uvw[6] );
	}

	MultiMatrix( r, uvw, gv );
	swap( uvw[1], uvw[4] );	
	swap( uvw[2], uvw[8] );	
	swap( uvw[6], uvw[9] );
	MultiMatrix( uvw, gv, m );
}

// 
// Global functions
//

int Copy( Mesh &dst, const Mesh &src ){
	int size;

	// clear dstination
	dst.Clear();

	// copy vertex
	dst.vertex_size = dst.vertex_max = src.vertex_size;
	size = 3 * dst.vertex_size * sizeof( float );
	dst.vertex_array = ( float * )malloc( size );
	memcpy( dst.vertex_array, src.vertex_array, size );

	// copy normal array
	dst.normal_size = dst.normal_max = src.normal_size;
	size = 3 * dst.normal_size * sizeof( float );
	dst.normal_array = ( float * )malloc( size );
	memcpy( dst.normal_array, src.normal_array, size );

	// copy texcoord array
	dst.tex_coord_size = dst.tex_coord_max = src.tex_coord_size;
	size = 2 * dst.tex_coord_size * sizeof( float );
	dst.tex_coord_array = ( float * )malloc( size );
	memcpy( dst.tex_coord_array, src.tex_coord_array, size );

	// copy triangle array
	dst.triangle_size = dst.triangle_max = src.triangle_size;
	size = dst.triangle_size * sizeof( Triangle );
	dst.triangle_array = ( Triangle * )malloc( size );
	memcpy( dst.triangle_array, src.triangle_array, size );

	// copy material array
	dst.material_size = dst.material_max = src.material_size;
	size = dst.material_size * sizeof( Material );
	dst.material_array = ( Material * )malloc( size );
	memcpy( dst.material_array, src.material_array, size );

	// copy group array
	dst.group_size = dst.group_max = src.group_size;
	size = dst.group_size * sizeof( Group );
	dst.group_array = ( Group * )malloc( size );
	memcpy( dst.group_array, src.group_array, size );

	// copy group triangle array
	dst.group_triangle_size = dst.group_triangle_max = src.group_triangle_size;
	size = dst.group_triangle_size * sizeof( int );
	dst.group_triangle_array = ( int* )malloc( size );
	memcpy( dst.group_triangle_array, src.group_triangle_array, size );

	// copy neighbor neighbor information
	dst.neighbor_array = NULL;
	dst.neighbor_vertex_array = NULL;
	dst.UpdateNeighbor();

	dst.Memory::Clear();
	return 1;
}

//
// Mesh
//

Mesh::Mesh( void ): Memory(){
	// all element is null
	vertex_array = NULL;
	normal_array = NULL;
	tex_coord_array = NULL;
	triangle_array = NULL;
	material_array = NULL;
	group_array = NULL;
	neighbor_array = NULL;
	neighbor_vertex_array = NULL;
	vertex_size = 
	normal_size = 
	tex_coord_size = 
	triangle_size = 
	material_size = 
	group_size = 
	group_triangle_size = 
	vertex_max = 
	neighbor_max = 
	normal_max = 
	tex_coord_max = 
	triangle_max = 
	material_max = 
	group_max = 
	group_triangle_max = 
	neighbor_vertex_max = 0;
	_texIsInit = false;
	// memory pool
//	_mem = NULL;
//	_m_size = 0;
}

Mesh::Mesh( const Mesh &m ): Memory(){
	// all element is null
	vertex_array = NULL;
	normal_array = NULL;
	tex_coord_array = NULL;
	triangle_array = NULL;
	material_array = NULL;
	group_array = NULL;
	neighbor_array = NULL;
	neighbor_vertex_array = NULL;
	vertex_size = 
	normal_size = 
	tex_coord_size = 
	triangle_size = 
	material_size = 
	group_size = 
	group_triangle_size = 
	vertex_max = 
	neighbor_max = 
	normal_max = 
	tex_coord_max = 
	triangle_max = 
	material_max = 
	group_max = 
	group_triangle_max = 
	neighbor_vertex_max = 0;
	_texIsInit = false;
	Copy( *this, m );
}	

Mesh::~Mesh( void ){
	Clear();
}

//
// File Processing
//

// decide array size and allacate memory for all array
void Mesh::FirstPass( FILE *objfile, char *filename ){
	char string[500];
	vertex_size = 
	normal_size = 
	tex_coord_size = 
	triangle_size = 
	material_size = 
	group_size = 
	group_triangle_size = 
	vertex_max = 
	neighbor_max = 
	normal_max = 
	tex_coord_max = 
	triangle_max = 
	material_max = 
	group_max = 
	group_triangle_max = 
	neighbor_vertex_max = 0;

	while( !feof( objfile ) ){
		fgets( string, 500, objfile );
		
		switch( string[0] ){	// check the first character
		case 'v':
			// check the second character
			switch( string[1] ){
			case ' ':		// v
				vertex_size++;
				break;
			case 'n':		// vn
				normal_size++;
				break;
			case 't':		// vt
				tex_coord_size++;
				break;
			}
			break;
		
		case 'f':	// f
			triangle_size++;
			break;

		case 'm':	// mtllib
			if( strcmp( string, "mtllib" ) > 0 ){	// if this string is mtllib
				char *ptr;
				ptr = strchr( string, '\n' );		// search the newline char
				if( ptr != NULL ) *ptr = '\0';		// add the null char
				ptr = strrchr( string, ':' );		// is it a local direction
				if( ptr != NULL ){					// global direction
					ReadMTL( string + 7 );
				}
				else{								// local direction
					ptr = strrchr( filename, '\\' );	// last direction
					if( ptr == NULL ){
						ReadMTL( string + 7 );
					}
					else{
						char path[500];
#if _NOT_VC6_
						strcpy_s( path, 500, filename );
#else
						strcpy( path, filename );
#endif
						ptr = strrchr( path, '\\' );
#if _NOT_VC6_
						sprintf_s( ptr + 1, 500 - ( ptr + 1 - path ), "%s", string + 7 );
#else	
						sprintf( ptr + 1, "%s", string + 7 );
#endif
						//sprintf( ptr + 1, "%s", string + 7 );
						ReadMTL( path );
					}
				}
			}
			break;

		case 'g':	// group
			group_size++;
			break;

		case '#':	// comment
		case 'u':	// usemtl
		default:
			// do nothing
			break;
		}
	}

	// memory allocate
	vertex_max = vertex_size;			// vertex
	normal_max = normal_size;			// normal
	tex_coord_max = tex_coord_size;		// texture
	triangle_max = triangle_size;		// tirangle
	group_max = group_size;				// group

	vertex_array = ( float* )malloc( vertex_max * 3 * sizeof( float ) );
	normal_array = ( float* )malloc( normal_max * 3 * sizeof( float ) );
	tex_coord_array = ( float* )malloc( tex_coord_max * 2 * sizeof( float ) );
	triangle_array = ( Triangle* )malloc( triangle_max * sizeof( Triangle ) );
	group_array = ( Group* )malloc( group_max * sizeof( Group ) );

	group_triangle_max = triangle_max;
	group_triangle_array = ( int* )malloc( group_triangle_max * sizeof( int ) );
}

// save all vertex, normal, texture coordinate, and triangle arrays
void Mesh::SecondPass( FILE *objfile ){
	char string[500], buf[10];
	int i;
	int triangle_index = 0;
	float *vertex = vertex_array, *normal = normal_array, *tex_coord = tex_coord_array;
	Triangle *triangle = triangle_array;
	int material_index = -1, group_index = -1;

	group_triangle_size = 0;

	while( !feof( objfile ) ){
		fgets( string, 500, objfile );
		
		switch( string[0] ){	// check the first character
		case 'v':
			// check the second character
			switch( string[1] ){
			case ' ':		// v
#if _IS_VC6_
				sscanf( string, "%s %f %f %f", buf, 
#else
				sscanf_s( string, "%s %f %f %f", buf, 10,
#endif
					vertex + 0, //vertex_array + vertex_index + 0,
					vertex + 1, //vertex_array + vertex_index + 1,
					vertex + 2 ); //vertex_array + vertex_index + 2 );
				vertex += 3;
				break;
			case 'n':		// vn
#if _IS_VC6_
				sscanf( string, "%s %f %f %f", buf,  
#else
				sscanf_s( string, "%s %f %f %f", buf, 10, 
#endif
					normal + 0, //normal_array + normal_index + 0,
					normal + 1, //normal_array + normal_index + 1,
					normal + 2 ); //normal_array + normal_index + 2 );
				normal += 3;
				break;
			case 't':		// vt
#if _IS_VC6_
				sscanf( string, "%s %f %f", buf, 
#else
				sscanf_s( string, "%s %f %f", buf, 10,
#endif
					tex_coord + 0, //tex_coord_array + tex_coord_index + 0,
					tex_coord + 1 ); //tex_coord_array + tex_coord_index + 1 );
				tex_coord += 2;
				break;
			}
			break;
		
		case 'f':	// f
		{
			// count information item
			char *ptr = strchr( string + 1, '/' );
			i = 0;
			while( ptr != NULL ){
				ptr = strchr( ptr + 1, '/' );
				i++;
			}

			// load face information
			if( i / 3 == 2 ){
				if( *( strchr( string, '/' ) + 1 ) == '/' ){		//  triangle//normal, missing tex-coord
#if _IS_VC6_
					sscanf( string, "%s %d//%d %d//%d %d//%d" , buf,
#else
					sscanf_s( string, "%s %d//%d %d//%d %d//%d" , buf, 10,
#endif
						&( triangle->vertex_index[0] ), //&( triangle_array[ triangle_index ].vertex_index[0] ),
						&( triangle->normal_index[0] ), //&( triangle_array[ triangle_index ].normal_index[0] ),
						&( triangle->vertex_index[1] ), //&( triangle_array[ triangle_index ].vertex_index[1] ),
						&( triangle->normal_index[1] ), //&( triangle_array[ triangle_index ].normal_index[1] ),
						&( triangle->vertex_index[2] ), //&( triangle_array[ triangle_index ].vertex_index[2] ),
						&( triangle->normal_index[2] ) ); //&( triangle_array[ triangle_index ].normal_index[2] ) );
					triangle->vertex_index[0]--; //triangle_array[ triangle_index ].vertex_index[0] --;
					triangle->normal_index[0]--; //triangle_array[ triangle_index ].normal_index[0] --;
					triangle->vertex_index[1]--; //triangle_array[ triangle_index ].vertex_index[1] --;
					triangle->normal_index[1]--; //triangle_array[ triangle_index ].normal_index[1] --;
					triangle->vertex_index[2]--; //triangle_array[ triangle_index ].vertex_index[2] --;
					triangle->normal_index[2]--; //triangle_array[ triangle_index ].normal_index[2] --;
					triangle->tex_coord_index[0] = -1;
				}
				else{
#if _IS_VC6_
					sscanf( string, "%s %d/%d/%d %d/%d/%d %d/%d/%d" , buf,
#else
					sscanf_s( string, "%s %d/%d/%d %d/%d/%d %d/%d/%d" , buf, 10,
#endif
						&( triangle->vertex_index[0] ), //&( triangle_array[ triangle_index ].vertex_index[0] ),
						&( triangle->tex_coord_index[0] ), //&( triangle_array[ triangle_index ].tex_coord_index[0] ),
						&( triangle->normal_index[0] ), //&( triangle_array[ triangle_index ].normal_index[0] ),
						&( triangle->vertex_index[1] ), //&( triangle_array[ triangle_index ].vertex_index[1] ),
						&( triangle->tex_coord_index[1] ), //&( triangle_array[ triangle_index ].tex_coord_index[1] ),
						&( triangle->normal_index[1] ), //&( triangle_array[ triangle_index ].normal_index[1] ),
						&( triangle->vertex_index[2] ), //&( triangle_array[ triangle_index ].vertex_index[2] ),
						&( triangle->tex_coord_index[2] ), //&( triangle_array[ triangle_index ].tex_coord_index[2] ),
						&( triangle->normal_index[2] ) ); //&( triangle_array[ triangle_index ].normal_index[2] ) );
					triangle->vertex_index[0]--; //triangle_array[ triangle_index ].vertex_index[0] --;
					triangle->tex_coord_index[0]--; //triangle_array[ triangle_index ].tex_coord_index[0] --;
					triangle->normal_index[0]--; //triangle_array[ triangle_index ].normal_index[0] --;
					triangle->vertex_index[1]--; //triangle_array[ triangle_index ].vertex_index[1] --;
					triangle->tex_coord_index[1]--; //triangle_array[ triangle_index ].tex_coord_index[1] --;
					triangle->normal_index[1]--; //triangle_array[ triangle_index ].normal_index[1] --;
					triangle->vertex_index[2]--; //triangle_array[ triangle_index ].vertex_index[2] --;
					triangle->tex_coord_index[2]--; //triangle_array[ triangle_index ].tex_coord_index[2] --;
					triangle->normal_index[2]--; //triangle_array[ triangle_index ].normal_index[2] --;
				}
			}
			else if( i / 3 == 1 ){
#if _IS_VC6_
				sscanf( string, "%s %d/%d %d/%d %d/%d" , buf,
#else
				sscanf_s( string, "%s %d/%d %d/%d %d/%d" , buf, 10,
#endif
					&( triangle->vertex_index[0] ), //&( triangle_array[ triangle_index ].vertex_index[0] ),
					&( triangle->tex_coord_index[0] ), //&( triangle_array[ triangle_index ].tex_coord_index[0] ),
					&( triangle->vertex_index[1] ), //&( triangle_array[ triangle_index ].vertex_index[1] ),
					&( triangle->tex_coord_index[1] ), //&( triangle_array[ triangle_index ].tex_coord_index[1] ),
					&( triangle->vertex_index[2] ), //&( triangle_array[ triangle_index ].vertex_index[2] ),
					&( triangle->tex_coord_index[2] ) ); //&( triangle_array[ triangle_index ].tex_coord_index[2] ) );
				triangle->vertex_index[0]--; //triangle_array[ triangle_index ].vertex_index[0] --;
				triangle->tex_coord_index[0]--; //triangle_array[ triangle_index ].tex_coord_index[0] --;
				triangle->vertex_index[1]--; //triangle_array[ triangle_index ].vertex_index[1] --;
				triangle->tex_coord_index[1]--; //triangle_array[ triangle_index ].tex_coord_index[1] --;
				triangle->vertex_index[2]--; //triangle_array[ triangle_index ].vertex_index[2] --;
				triangle->tex_coord_index[2]--; //triangle_array[ triangle_index ].tex_coord_index[2] --;
				triangle->normal_index[0] = -1;
			}
			else{
#if _IS_VC6_
				sscanf( string, "%s %d %d %d" , buf,
#else
				sscanf_s( string, "%s %d %d %d" , buf, 10,
#endif
					&( triangle->vertex_index[0] ), //&( triangle_array[ triangle_index ].vertex_index[0] ),
					&( triangle->vertex_index[1] ), //&( triangle_array[ triangle_index ].vertex_index[1] ),
					&( triangle->vertex_index[2] ) ); //&( triangle_array[ triangle_index ].vertex_index[2] ) );
				triangle->vertex_index[0]--; //triangle_array[ triangle_index ].vertex_index[0] --;
				triangle->vertex_index[1]--; //triangle_array[ triangle_index ].vertex_index[1] --;
				triangle->vertex_index[2]--; //triangle_array[ triangle_index ].vertex_index[2] --;
				triangle->tex_coord_index[0] = triangle->normal_index[0] = -1;
			}

			// load material index
			triangle->material_index = material_index;
			//triangle_array[ triangle_index ].material_index = material_index;

			// count group triangles
			if( group_size ) {
				group_array[ group_index ].triangle_size++;		// count triangles
				group_triangle_array[ group_triangle_size++ ] = triangle_index;	// record triangle
			}
			triangle_index++;
			triangle++;
			break;
		}

	case 'u':	// usemtl
		{
			char mtlname[500];
#if _NOT_VC6_
			sscanf_s( string, "%s %s", buf, 10, mtlname, 500 );
#else
			sscanf( string, "%s %s", buf, mtlname );
#endif
			for( i = 0; i < material_size; i++ )	// find the material index
				if( strcmp( mtlname, material_array[i].name ) == 0 )	break;	
			if( i < material_size ) material_index = i;
			break;
		}

		case 'g':	// group
			group_index ++;
#if _NOT_VC6_
			sscanf_s( string, "%s %s", buf, 10, group_array[ group_index ].name, 500 );
#else
			sscanf( string, "%s %s", buf, group_array[ group_index ].name );
#endif
			group_array[ group_index ].triangle_index = group_triangle_array + group_triangle_size;
			group_array[ group_index ].triangle_size = 0;
			break;

		case '#':	// comment
		default:
			// do nothing
			break;
		}
	}
}

//
// Update Processing
//

#define vertex( i, j ) triangle_array[(i)].vertex_index[(j)]
void Mesh::UpdateNeighbor( void ){
	int i, j, k, *size;

	// memory allocate
	if( neighbor_max < vertex_max ) {
		neighbor_max = vertex_max;
		neighbor_array = ( NeighborVertex* )realloc( 
			neighbor_array, neighbor_max * sizeof( NeighborVertex ) );
	}
	if( neighbor_vertex_max < triangle_size * 6 ){
		neighbor_vertex_max = triangle_size * 6;
		neighbor_array[0].vertex_index = neighbor_vertex_array =
			( int* )realloc( 
			neighbor_vertex_array, neighbor_vertex_max * sizeof( int ) );
	}

	// all element is 0	
	for( i = 0; i < neighbor_max; i++ ) neighbor_array[i].size = 0;

	// memory pool for tempory size array
	if( _m_size < vertex_size * sizeof( int ) ){
		_m_size = vertex_size * sizeof( int );
		_mem = realloc( _mem, _m_size );
	}
	size = ( int* )_mem;

	// initialize size array
	for( i = 0; i < vertex_size; i++ ) size[i] = 0;

	// first pass: counting
	for( i = 0; i < triangle_size; i++ ){
		// each triangle vertex has 2 neighbor verteces
		for( j = 0; j < 3; j++ )
			neighbor_array[ vertex( i, j ) ].size += 2;
	}

	// set offset of array of neighbor vertex indeces
	for( i = 1; i < vertex_size; i++ ){
		neighbor_array[i].vertex_index = neighbor_array[i - 1].vertex_index + neighbor_array[i - 1].size;
	}

	// second pass: recording
	for( i = 0; i < triangle_size; i++ ){
		for( j = 0; j < 3; j++ ){
			// check the repeat element
			for( k = 0; k < size[ vertex( i, j ) ]; k++ ){
				if( neighbor_array[ vertex( i, j ) ].vertex_index[k] == vertex( i, ( j + 1 ) % 3 ) ) break;
			}
			if( k == size[ vertex( i, j ) ] ){		// if not be recorded yet
				neighbor_array[ vertex( i, j ) ].vertex_index[k] = vertex( i, ( j + 1 ) % 3 );
				size[ vertex( i, j ) ]++;
			}
			else neighbor_array[ vertex( i, j ) ].size--;		// subtract 1 repeat element

			// check the repeat element
			for( k = 0; k < size[ vertex( i, j ) ]; k++ ){
				if( neighbor_array[ vertex( i, j ) ].vertex_index[k] == vertex( i, ( j + 2 ) % 3 ) ) break;
			}
			if( k == size[ vertex( i, j ) ] ){		// if not be recorded yet
				neighbor_array[ vertex( i, j ) ].vertex_index[k] = vertex( i, ( j + 2 ) % 3 );
				size[ vertex( i, j ) ]++;
			}
			else neighbor_array[ vertex( i, j ) ].size--;		// subtract 1 repeat element
		}
	}
}
#undef vertex

void Mesh::UpdateTriangle( int v ){
	int i, j, *match;
	Triangle *t, *t2;

	// memory pool
	i = triangle_size * sizeof( int );
	if( i > ( int )_m_size ) _mem = realloc( _mem, i );
	match = ( int* )_mem;

	// initialize match table
	for( i = 0; i < triangle_size; i++ ) match[i] = i;

	for( i = 0; i < triangle_size; ){	// for every triangles
		t = triangle_array + i;
		if( t->vertex_index[0] == v ||
			t->vertex_index[1] == v ||
			t->vertex_index[2] == v ){		// include deleted vertex
			triangle_size--;
			match[i] = -1;
			// delete triangle and offset
			for( j = i; j < triangle_size; j++ ){
				t = triangle_array + j;
				t2 = t + 1;
				t->material_index = t2->material_index;
				t->normal_index[0] = t2->normal_index[0]; 
				t->normal_index[1] = t2->normal_index[1]; 
				t->normal_index[2] = t2->normal_index[2];
				t->tex_coord_index[0] = t2->tex_coord_index[0]; 
				t->tex_coord_index[1] = t2->tex_coord_index[1]; 
				t->tex_coord_index[2] = t2->tex_coord_index[2];
				t->vertex_index[0] = t2->vertex_index[0]; 
				t->vertex_index[1] = t2->vertex_index[1]; 
				t->vertex_index[2] = t2->vertex_index[2];
				match[j + 1]--;
			}
			continue;
		}
		if( t->vertex_index[0] > v ) t->vertex_index[0]--;
		if( t->vertex_index[1] > v ) t->vertex_index[1]--;
		if( t->vertex_index[2] > v ) t->vertex_index[2]--;
		i++;
	}

	UpdateGroup( match );
}

void Mesh::UpdateGroup( int *m ){
}

//
// File Processing
//

#if _NOT_VC6_

bool Mesh::ReadOBJ( wchar_t *filename ){
	size_t total, bytes;
	char *buf = NULL;
	bool ret;
	// count bytes needed
	bytes = WideCharToMultiByte( 950, 0, filename, -1, buf, 0, 0, 0 );
	// mamory
	buf = ( char* )malloc( bytes );
	if( buf == NULL ) return false;
	// convert widechar to multibytes
	wcstombs_s( &total, buf, bytes, filename, _TRUNCATE );
	ret = ReadOBJ( buf );	// main function
	free( buf );
	return ret;
}

bool Mesh::ReadMTL( wchar_t *filename ){
	size_t total, bytes;
	char *buf = NULL;
	bool ret;
	// count bytes needed
	bytes = WideCharToMultiByte( 950, 0, filename, -1, buf, 0, 0, 0 );
	// mamory
	buf = ( char* )malloc( bytes );
	if( buf == NULL ) return false;
	// convert widechar to multibytes
	wcstombs_s( &total, buf, bytes, filename, _TRUNCATE );
	ret = ReadMTL( buf );	// main function
	free( buf );
	return ret;
}

bool Mesh::WriteOBJ( wchar_t *filename, int mode ){
	size_t total, bytes;
	char *buf = NULL;
	bool ret;
	// count bytes needed
	bytes = WideCharToMultiByte( 950, 0, filename, -1, buf, 0, 0, 0 );
	// mamory
	buf = ( char* )malloc( bytes );
	if( buf == NULL ) return false;
	// convert widechar to multibytes
	wcstombs_s( &total, buf, bytes, filename, _TRUNCATE );
	ret = WriteOBJ( buf, mode );	// main function
	free( buf );
	return ret;
}

#endif

bool Mesh::ReadOBJ( char *filename ){
	char mtlfile[500] = "";
	FILE *objfile;
	_texIsInit = false;
#if _NOT_VC6_
	fopen_s( &objfile, filename, "r+" );
#else
	objfile = fopen( filename, "r+" );
#endif
	if( objfile == NULL ) {
		printf( "Can't open OBJ file: %s\n", filename );
		return false;	// fail to open file
	}

	FirstPass( objfile, filename );
	rewind( objfile );
	SecondPass( objfile );

	fclose( objfile );

	UpdateNeighbor();
	return true;
}

bool Mesh::ReadOBJ(std::string filename)
{
	char tmpChar[1024];
	strcpy(tmpChar, filename.c_str());
	return ReadOBJ(tmpChar);
}

bool Mesh::ReadMTL( char *filename ){
	FILE *mtlfile;
#if _NOT_VC6_
	fopen_s( &mtlfile, filename, "r" );
#else
	mtlfile = fopen( filename, "r" );
#endif
	if( mtlfile == NULL ) {
		printf( "Can't open MTL file: %s\n", filename );
		return false;	// fail to open file
	}

	char string[500], buf[10];
	int index = 0;
	ImageR< BYTE > imgbuf;
	material_size = 0;

	// first pass
	while( !feof( mtlfile ) ){
		fgets( string, 500, mtlfile );
		if( string[0] != 'n' ) continue;

		// check if see "newmtl"
		string[6] = '\0';
		if( strcmp( string, "newmtl" ) == 0 ) material_size++;
	}

	// memory allocate to material array
	material_max = material_size;
	material_array = ( Material* )calloc( material_size, sizeof( Material ) );
	for( index = 0; index < material_size; index++ )
		material_array[ index ].map.bits = NULL;			// set empty image
	index = -1;

	// second pass
	rewind( mtlfile );
	while( !feof( mtlfile ) ){
		fgets( string, 500, mtlfile );

		// eliminate the whitespace character
		if( string[0] == '\t' ){
			for( int i = 0; i < (int)strlen( string ) + 1; i++ ) string[i] = string[i + 1];
		}

		switch( string[0] ){
		case 'n':	// "newmtl"
			index++;
#if _NOT_VC6_
			sscanf_s( string, "%s %s", buf, 10, material_array[ index ].name, 500 );
#else
			sscanf( string, "%s %s", buf, material_array[ index ].name );
#endif
			break;

		case 'N':
			switch( string[1] ){
			case 's':
#if _NOT_VC6_
				sscanf_s( string, "%s %f", buf, 10, &material_array[ index ].Ns );
#else
				sscanf( string, "%s %f", buf, &material_array[ index ].Ns );
#endif
				break;
			
			case 'i':
			default:
				;	// do nothing
			}
			break;

		case 'd':
#if _NOT_VC6_
			sscanf_s( string, "%s %f", buf, 10, &material_array[ index ].d );
#else
			sscanf( string, "%s %f", buf, &material_array[ index ].d );
#endif
			break;

		case 'K':
			switch( string[1] ){
			case 'a':
#if _NOT_VC6_
				sscanf_s( string, "%s %f %f %f", buf, 10,
#else
				sscanf( string, "%s %f %f %f", buf,
#endif
					&material_array[ index ].Ka[0],
					&material_array[ index ].Ka[1],
					&material_array[ index ].Ka[2] );
				material_array[ index ].Ka[3] = 1;
				break;

			case 'd':
#if _NOT_VC6_
				sscanf_s( string, "%s %f %f %f", buf, 10,
#else
				sscanf( string, "%s %f %f %f", buf,
#endif
					&material_array[ index ].Kd[0],
					&material_array[ index ].Kd[1],
					&material_array[ index ].Kd[2] );
				material_array[ index ].Kd[3] = 1;
				break;
				
			case 's':
#if _NOT_VC6_
				sscanf_s( string, "%s %f %f %f", buf, 10,
#else
				sscanf( string, "%s %f %f %f", buf,
#endif
					&material_array[ index ].Ks[0],
					&material_array[ index ].Ks[1],
					&material_array[ index ].Ks[2] );
				material_array[ index ].Ks[3] = 1;
				break;
				
			case 'e':
			default:
				;	// do nothing
			}
			break;

		case 'm':
			if( strcmp( string, "map_Kd" ) > 0 ){	// map_Kd 
				if( material_array[ index ].map.bits == NULL ){
					char *ptr;
					ptr = strchr( string, '\n' );		// search the newline char
					if( ptr != NULL ) *ptr = '\0';		// add the null char
					ptr = strrchr( string, ':' );		// is it a local direction
					if( ptr != NULL ){					// global direction
						if( imgbuf.ReadFile( string + 7 ) == 0 ) break;
//						material_array[ index ].map = new IMAGE( string + 7 );
						material_array[ index ].map.bits = imgbuf.GetImage();
						material_array[ index ].map.w = imgbuf.Width();
						material_array[ index ].map.h = imgbuf.Height();
						material_array[ index ].map.stride = imgbuf.Stride();
						material_array[ index ].map.color = imgbuf.Color();
					}
					else{								// local direction
						ptr = strrchr( filename, '\\' );	// last direction
						if( ptr == NULL ){
							if( imgbuf.ReadFile( string + 7 ) == 0 ) break;
//							material_array[ index ].map = new IMAGE( string + 7 );
							material_array[ index ].map.bits = imgbuf.GetImage();
							material_array[ index ].map.w = imgbuf.Width();
							material_array[ index ].map.h = imgbuf.Height();
							material_array[ index ].map.stride = imgbuf.Stride();
							material_array[ index ].map.color = imgbuf.Color();
						}
						else{
							char path[500];
#if _NOT_VC6_
							strcpy_s( path, 500, filename );
#else
							strcpy( path, filename );
#endif
							ptr = strrchr( path, '\\' );
#if _NOT_VC6_
							sprintf_s( ptr + 1, 500 - ( int )( ptr + 1 - path ), string + 7 );
#else
							sprintf( ptr + 1, "%s", string + 7 );
#endif
							if( imgbuf.ReadFile( path ) == 0 ) break;
//							material_array[ index ].map = new IMAGE( path );
							material_array[ index ].map.bits = imgbuf.GetImage();
							material_array[ index ].map.w = imgbuf.Width();
							material_array[ index ].map.h = imgbuf.Height();
							material_array[ index ].map.stride = imgbuf.Stride();
							material_array[ index ].map.color = imgbuf.Color();
						}
					}

					if( material_array[ index ].map.bits == NULL ) {
//						delete material_array[ index ].map;
						material_array[ index ].map.bits = NULL;
					}
				}
			}
			imgbuf.SetEmpty();
			break;
			
		case 'T':
		case 'i':
		default:
			;	// do nothing
		}
	}

	fclose( mtlfile );

	return true;
}

bool Mesh::WriteOBJ( char *filename, int mode ){
	FILE *outfile;
	float *p;
	Triangle *pt;
	Group *pg;
	int i, j, temp;

	// open file
	#if _NOT_VC6_
	fopen_s( &outfile, filename, "w" );
	#else 
	outfile = fopen( filename, "w" );
	#endif
	if( outfile == NULL ) return false;
	
	// output vertex
	p = vertex_array;
	fprintf_s( outfile, " # %d verteces\n", vertex_size );
	for( i = 0; i < vertex_size; i++, p += 3 )
		fprintf_s( outfile, "v %f %f %f\n", p[0], p[1], p[2] );
	fprintf_s( outfile, "\n\n" );
	
	// output normal
	if( mode & OBJ_NORMAL ){
		p = normal_array;
		fprintf_s( outfile, " # %d normal vectors\n", normal_size );
		for( i = 0; i < normal_size; i++, p += 3 )
			fprintf_s( outfile, "vn %f %f %f\n", p[0], p[1], p[2] );
	}
	fprintf_s( outfile, "\n\n" );

	// output texture coordinate
	if( mode & OBJ_TEXTURE_COORD ){
		p = tex_coord_array;
		fprintf_s( outfile, "# %d texture coordinate\n", tex_coord_size );
		for( i = 0; i < tex_coord_size; i++, p += 3 )
			fprintf_s( outfile, "vt %f %f %f\n", p[0], p[1], p[2] );
		if( i == 0 ){
			mode &= ~OBJ_TEXTURE_COORD;
		}
	}
	fprintf_s( outfile, "\n\n" );

	// output triangles
	if( mode & OBJ_TRIANGLE ){

		// grouping
		if( mode & OBJ_GROUP ){
			fprintf_s( outfile, "#grouping\n\n" );
			pg = group_array;
			for( j = 0; j < group_size; j++, pg++ ){
				fprintf_s( outfile, "#group %d %s: %d triangles\n", j + 1, pg->name, pg->triangle_size );
				fprintf_s( outfile, "g %s\n", pg->name );
				temp = -1;
				// face / texture / normal
				if( ( mode & OBJ_TEXTURE_COORD ) && ( mode & OBJ_NORMAL ) ){
					for( i = 0; i < pg->triangle_size; i++ ){
						pt = triangle_array + pg->triangle_index[i];
						// change material
						if( temp != pt->material_index ){
							temp = pt->material_index;
							fprintf_s( outfile, "ustmtl %s\n", material_array[temp].name );
						}
						fprintf_s( outfile, "f %d/%d/%d %d/%d/%d %d/%d/%d\n", 
						pt->vertex_index[0] + 1, pt->tex_coord_index[0] + 1, pt->normal_index[0] + 1,
						pt->vertex_index[1] + 1, pt->tex_coord_index[1] + 1, pt->normal_index[1] + 1,
						pt->vertex_index[2] + 1, pt->tex_coord_index[2] + 1, pt->normal_index[2] + 1 );
					}
				}
				// face / texture
				else if( mode & OBJ_TEXTURE_COORD ){
					for( i = 0; i < pg->triangle_size; i++ ){
						pt = triangle_array + pg->triangle_index[i];
						// change material
						if( temp != pt->material_index ){
							temp = pt->material_index;
							fprintf_s( outfile, "ustmtl %s\n", material_array[temp].name );
						}
						fprintf_s( outfile, "f %d/%d %d/%d %d/%d\n", 
						pt->vertex_index[0] + 1, pt->tex_coord_index[0] + 1,
						pt->vertex_index[1] + 1, pt->tex_coord_index[1] + 1, 
						pt->vertex_index[2] + 1, pt->tex_coord_index[2] + 1 );
					}
				}
				// face / / normal 
				else if( mode & OBJ_NORMAL ){
					for( i = 0; i < pg->triangle_size; i++ ){
						pt = triangle_array + pg->triangle_index[i];
						// change material
						if( temp != pt->material_index ){
							temp = pt->material_index;
							fprintf_s( outfile, "ustmtl %s\n", material_array[temp].name );
						}
						fprintf_s( outfile, "f %d//%d %d//%d %d//%d\n", 
						pt->vertex_index[0] + 1, pt->normal_index[0] + 1,
						pt->vertex_index[1] + 1, pt->normal_index[1] + 1, 
						pt->vertex_index[2] + 1, pt->normal_index[2] + 1 );
					}
				}
				// face 
				else{
					for( i = 0; i < pg->triangle_size; i++ ){
						fprintf_s( outfile, "f %d %d %d\n",	pt->vertex_index[0] + 1, pt->vertex_index[1] + 1, pt->vertex_index[2] + 1 );
					}
				}
			}
		}

		// non grouping
		else{
			pt = triangle_array;
			fprintf_s( outfile, "# %d triangles\n", triangle_size );
			// face / texture / normal
			if( ( mode & OBJ_TEXTURE_COORD ) && ( mode & OBJ_NORMAL ) ){
				for( i = 0; i < triangle_size; i++, pt++ ){
					fprintf_s( outfile, "f %d/%d/%d %d/%d/%d %d/%d/%d\n", 
					pt->vertex_index[0] + 1, pt->tex_coord_index[0] + 1, pt->normal_index[0] + 1,
					pt->vertex_index[1] + 1, pt->tex_coord_index[1] + 1, pt->normal_index[1] + 1,
					pt->vertex_index[2] + 1, pt->tex_coord_index[2] + 1, pt->normal_index[2] + 1 );
				}
			}
			// face / texture 
			else if( mode & OBJ_TEXTURE_COORD ){
				for( i = 0; i < triangle_size; i++, pt++ ){
					fprintf_s( outfile, "f %d/%d %d/%d %d/%d\n", 
					pt->vertex_index[0] + 1, pt->tex_coord_index[0] + 1,
					pt->vertex_index[1] + 1, pt->tex_coord_index[1] + 1, 
					pt->vertex_index[2] + 1, pt->tex_coord_index[2] + 1 );
				}
			}
			// face / / normal 
			else if( mode & OBJ_NORMAL ){
				for( i = 0; i < triangle_size; i++, pt++ ){
					fprintf_s( outfile, "f %d//%d %d//%d %d//%d\n", 
					pt->vertex_index[0] + 1, pt->normal_index[0] + 1,
					pt->vertex_index[1] + 1, pt->normal_index[1] + 1, 
					pt->vertex_index[2] + 1, pt->normal_index[2] + 1 );
				}
			}
			// face 
			else{
				for( i = 0; i < triangle_size; i++, pt++ ){
					fprintf_s( outfile, "f %d %d %d\n",
						pt->vertex_index[0] + 1, pt->vertex_index[1] + 1, pt->vertex_index[2] + 1 );
				}
			}
		}
	}
	fclose( outfile );
	return true;
}

//
// Mesh Functions
//

int Mesh::DeleteVertex( int v ){
	if( v < 0 || v >= vertex_size ) return 1;

	int i;
	UpdateTriangle( v );

	// delete vertex
	vertex_size--;
	for( i = v; i < vertex_size; i++ ) {
		vertex_array[ i * 3 + 0 ] = vertex_array[ (i + 1 ) * 3 + 0 ];
		vertex_array[ i * 3 + 1 ] = vertex_array[ (i + 1 ) * 3 + 1 ];
		vertex_array[ i * 3 + 2 ] = vertex_array[ (i + 1 ) * 3 + 2 ];
	}
	
	UpdateNeighbor();
	return 0;
}

int Mesh::DeleteTriangle( int t ){
	if( t < 0 || t >= triangle_size ) return 1;

	int i;
	Triangle *t1, *t2;

	// triangle offset
	triangle_size--;
	for( i = t; i < triangle_size; i++ ){
		t1 = triangle_array + i;
		t2 = t1 + 1;
		t1->material_index = t2->material_index;
		t1->normal_index[0] = t2->normal_index[0]; 
		t1->normal_index[1] = t2->normal_index[1]; 
		t1->normal_index[2] = t2->normal_index[2];
		t1->tex_coord_index[0] = t2->tex_coord_index[0]; 
		t1->tex_coord_index[1] = t2->tex_coord_index[1]; 
		t1->tex_coord_index[2] = t2->tex_coord_index[2];
		t1->vertex_index[0] = t2->vertex_index[0]; 
		t1->vertex_index[1] = t2->vertex_index[1]; 
		t1->vertex_index[2] = t2->vertex_index[2];
	}

	UpdateNeighbor();
	return 0;
}
/**/
int Mesh::DeleteTriangle( Group &g ){
	if( g.triangle_size == 0 || g.triangle_index == NULL ) return 1;

	int i;
	for( i = 0; i < g.triangle_size; i++ )
		triangle_array[ g.triangle_index[i] ].vertex_index[0] = -1;

	return 0;
}
/**/
//
// Advanced Functions
//

std::vector< Group > Mesh::ContinueTriangles( void ){
	int i, j, k, a, b, c, n, 
		*group, *edge;
	std::vector< Group > g(1);
	NeighborVertex temp;

	// memory pool setting
	a = max( vertex_size, triangle_size );
	i = 2 * a * sizeof( int );
	if( ( int )_m_size < i ){
		_m_size = ( unsigned char )i;
		_mem = realloc( _mem, _m_size );
	}
	// initialize
	group = ( int* )_mem;
	edge = group + a;
	for( i = 0; i < vertex_size; i++ ) group[i] = 0;
	n = 0;

	for( k = 0; k < vertex_size; k++ ){
		if( group[k] ) continue;
		// another group
		group[i] = ++n;
		edge[0] = k;
		a = 0;
		b = 1;
		do{	
			c = 0;
			for( i = a; i < b; i++ ){		// for every edge vertex
				temp = neighbor_array[ edge[i] ];
				// check edge's neighbor
				for( j = 0; j < temp.size; j++ ){
					if( group[ temp.vertex_index[j] ] == 0 ) {
						// select to group n
						group[ temp.vertex_index[j] ] = n;
						edge[ b + c ] = temp.vertex_index[j]; 
						c++;
					}
				}
			}
			a = b;
			b += c;
		} while( c != 0 );		// still have edge
	}

	g.resize( n );

	// translate result to triangle
	for( i = 0; i < vertex_size; i++ ) edge[i] = group[i];
	for( i = 0; i < triangle_size; i++ )
		group[i] = edge[ triangle_array[i].vertex_index[0] ];

	// set to grop vector
	for( i = 1; i <= n; i++ ){
		for( j = a = 0; j < triangle_size; j++ ){
			if( group[j] == i ) a++;
		}
		g[i - 1].triangle_index = ( int* )malloc( a * sizeof( int ) );
		g[i - 1].triangle_size = a;
		for( j = a = 0; j < triangle_size; j++ ){
			if( group[j] == i ) g[i - 1].triangle_index[a++] = j;
			if( a == g[i - 1].triangle_size ) break;
		}
	}

	return g;
}

std::vector< int > Mesh::NeighborTriangles( const Triangle &t ){
	int i, j;
	bool *neighbor;
	std::vector< int > trivec;

	// memory pool
	if( _m_size < vertex_size * sizeof( bool ) ){
		_m_size = vertex_size * sizeof( bool );
		_mem = realloc( _mem, _m_size );
	}
	neighbor = ( bool* )_mem;

	// initialize
	for( i = 0; i < vertex_size; i++ ) neighbor[i] = false;
	for( i = 0; i < 3; i++ ) {
		neighbor[ t.vertex_index[i] ] = true;	// verteces in triangle
		for( j = 0; j < neighbor_array[ t.vertex_index[i] ].size; j++ ){	// neighbor of each vertex
			neighbor[ neighbor_array[ t.vertex_index[i] ].vertex_index[j] ] = true;
		}
	}

	// searching every triangle
	for( i = 0; i < triangle_size; i++ ){
		for( j = 0; j < 3; j++ ){
			if( !neighbor[ triangle_array[i].vertex_index[j] ] ) break;
		}
		if( j == 3 )		// is a neighbor of this triangle
			trivec.push_back( i );
	}
	return trivec;
}

std::vector< int > Mesh::NeighborTriangles( int t ){
	if( t < 0 || t >= triangle_size ) return std::vector< int >();
	return NeighborTriangles( triangle_array[t] );
}

int Mesh::DeleteContinueTriangles( int threshold ){
	int i, j, n;
	vector< Group >g = ContinueTriangles();

	// deleting
	n = g.size();
	for( i = j = 0; i < n; i++ ){
		if( g[i].triangle_size < threshold ) {
			DeleteTriangle( g[i] );
			j++;
		}
	}

	// free memory space
	for( i = 0; i < n; i++ ) free( g[i].triangle_index );

	// updating
	if( j > 0 ) {			// if did delete triangles
		UpdateTriangle();
		return 0;
	}
	else return 1;			// if didn't delete any triangle, return 1
}

//
// Updating
//

int Mesh::UpdateTriangle( void ){
	int i, jmp = 0;
	Triangle *t1, *t2;
		
	for( i = 0; i < triangle_size; i++ ){
		if( jmp == 0 && triangle_array[i + jmp].vertex_index[0] != -1 ) continue;	// don't need to delete
		while( triangle_array[i + jmp].vertex_index[0] == -1 ){							// if jump to a deleted vertex
			jmp++; 
			triangle_size--;
		}
		t1 = triangle_array + i;
		t2 = t1 + jmp;
		t1->material_index = t2->material_index;
		t1->normal_index[0] = t2->normal_index[0]; 
		t1->normal_index[1] = t2->normal_index[1]; 
		t1->normal_index[2] = t2->normal_index[2];
		t1->tex_coord_index[0] = t2->tex_coord_index[0]; 
		t1->tex_coord_index[1] = t2->tex_coord_index[1]; 
		t1->tex_coord_index[2] = t2->tex_coord_index[2];
		t1->vertex_index[0] = t2->vertex_index[0]; 
		t1->vertex_index[1] = t2->vertex_index[1]; 
		t1->vertex_index[2] = t2->vertex_index[2];
	}
	UpdateNeighbor();

	if( jmp ) return 0;
	else return 1;
}

//
// Basic Function
//

void Mesh::Clear( void ){
	int i;

	// free vertex array
	if( vertex_array ) {
		free( vertex_array );
		vertex_array = NULL;
		vertex_size = vertex_max = 0;
	}

	// free normal array
	if( normal_array ) {
		free( normal_array );
		normal_array = NULL;
		normal_size = normal_max = 0;
	}

	// free texture coordinate array
	if( tex_coord_array ) {
		free( tex_coord_array );
		tex_coord_array = NULL;
		tex_coord_size = tex_coord_max = 0;
	}

	// free triangle array
	if( triangle_array ) {
		free( triangle_array );
		triangle_array = NULL;
		triangle_size = triangle_max = 0;
	}

	// free material array
	if( material_array ) {
		for( i = 0; i < material_max; i++ ){
			// free texture image
			if( material_array[i].map.bits != NULL ) free( material_array[i].map.bits );
		}
		free( material_array );
		material_array = NULL;
		material_size = material_max = 0;
	}

	// free group array
	if( group_array ) {
		free( group_array );
		free( group_triangle_array );
		group_array = NULL;
		group_triangle_array = NULL;
		group_size = group_max = group_triangle_size = group_triangle_max = 0;
	}

	// free neighbor array
	if( neighbor_array ){
		free( neighbor_vertex_array );
		free( neighbor_array );
		neighbor_vertex_array = NULL;
		neighbor_array = NULL;
		neighbor_max = neighbor_vertex_max = 0;
	}

	this->Memory::Clear();
}

void Mesh::Scale( float x, float y, float z ){
	int i;
	float *p = vertex_array;

	for( i = 0; i < vertex_size; i++, p += 3 ){
		p[0] *= ( float )x;
		p[1] *= ( float )y;
		p[2] *= ( float )z;
	}
}

void Mesh::Translate( float x, float y, float z ){
	int i;
	for( i = 0; i < vertex_size; i++ ){
		vertex_array[i * 3 + 0] += ( float )x;
		vertex_array[i * 3 + 1] += ( float )y;
		vertex_array[i * 3 + 2] += ( float )z;
	}
}

void Mesh::Rotate( float deg, float x, float y, float z ){
	float m[16], v[4] = { 0.0, 0.0, 0.0, 1.0 };
	::Rotate( m, deg, x, y, z );
	MultiMatrix( m );
}

void Mesh::CenterRotate( float deg, float x, float y, float z ){
	float m[3][16], v[4] = { 0.0, 0.0, 0.0, 1.0 }, *pv;
	int i;

	pv = vertex_array;
	for( i = 0; i < vertex_size; i++, pv += 3 ){
		v[0] += pv[0];	v[1] += pv[1];	v[2] += pv[2];
	}
	v[3] = 1.0f / ( float )vertex_size;
	v[0] *= v[3];
	v[1] *= v[3];
	v[2] *= v[3];

	// T(v) * R * T(-v)
	::Translate( m[0], v[0], v[1], v[2] );
	::Rotate( m[1], deg, x, y, z );
	::MultiMatrix( m[0], m[1], m[2] );
	::Translate( m[0], -v[0], -v[1], -v[2] );
	::MultiMatrix( m[2], m[0], m[1] );
	MultiMatrix( m[1] );
}

void Mesh::MultiMatrix( float m[16] ){
	int i, j;
	float *p;
	float inv[16], temp[16], v[4] = { 1.0, 1.0, 1.0, 1.0 };

	// inverse of m
	memcpy( temp, m, 16 * sizeof( float ) );	// copy
	invf( inv, m, 4 );
	memcpy( m, temp, 16 * sizeof( float ) );	// copy back
	// transpose
	swap( inv[1], inv[4] );	swap( inv[2], inv[8] );	swap( inv[3], inv[12] );
	swap( inv[6], inv[9] );	swap( inv[7], inv[13] );
	swap( inv[11], inv[14] );
	// transformation
	p = vertex_array;
	for( i = 0; i < vertex_size; i++, p += 3 ){
		for( j = 0; j < 3; j++ ) v[j] = p[j];
		MultiVector( m, v );
		for( j = 0; j < 3; j++ ) p[j] = v[j];
	}
	// normal transformation
	p = normal_array;
	for( i = 0; i < normal_size; i++, p += 3 ){
		v[3] = 1.f;
		for( j = 0; j < 3; j++ ) v[j] = p[j];
		MultiVector( inv, v );
		v[3] = 1.f / v[3];
		for( j = 0; j < 3; j++ ) p[j] = v[j];// * v[3];
	}
}

void Mesh::Smooth( void ){
	int i, j, k;
	Double3D *nv;
	
	// memory pool updating
	k = vertex_size * sizeof( Double3D );
	if( k > ( int )_m_size ){
		_m_size = ( unsigned int )k;
		_mem = realloc( _mem, _m_size );
	}
	nv = ( Double3D* )_mem;

	for( i = 0; i < vertex_size; i++ ){
		// initialize
		for( j = 0; j < 3; j++ ) nv[i].p[j] = ( double )vertex_array[i * 3 + j];
		// sum coordinate for every neighbor vertex
		for( j = 0; j < neighbor_array[i].size; j++ ){
			for( k = 0; k < 3; k++ )
				nv[i].p[k] += ( double )vertex_array[ neighbor_array[i].vertex_index[j] * 3 + k ];
		}
		// averaging
		for( j = 0; j < 3; j++ ) nv[i].p[j] /= (double)( neighbor_array[i].size + 1 );
	}

	// updating smoothed verteces
	for( i = 0; i < vertex_size; i++ ) {
		for( j = 0; j < 3; j++ ) vertex_array[i * 3 + j] = ( float )nv[i].p[j];
	}
}

//
// OpenGL Drawing Functions
//

void Mesh::Draw( int mode ){
	unsigned char tex2d, clr_mat;
	int shade;
	int i, j, image = -1;
	//GLuint tex[32];
	float buf[4];
	Triangle *triangle;

	// mode checking
	if( vertex_size == 0 ) return;
	if( triangle_size == 0 ) mode &= MESH_VERTEX;
	if( normal_size == 0 ) mode &= ~MESH_NORMAL;
	if( tex_coord_size == 0 ) mode &= ~MESH_TEXTURE;
	if( material_size == 0 ) mode &= ~MESH_MATERIAL;

	if( mode & MESH_DEPTH_MAP ) {
		DrawDepthMap();
		return;
	}
	glGetIntegerv( GL_SHADE_MODEL, &shade );
	if( mode & MESH_SMOOTH ) glShadeModel( GL_SMOOTH );
	if( mode & MESH_FLAT ) glShadeModel( GL_FLAT );

	glGetBooleanv( GL_COLOR_MATERIAL, &clr_mat );
	glDisable( GL_COLOR_MATERIAL );

	glGetBooleanv( GL_TEXTURE_2D, &tex2d );
	if( mode & MESH_TEXTURE ){
		glEnable( GL_TEXTURE_2D );
		if (!_texIsInit)
		{
			glGenTextures(material_size, _tex);
			for (i = 0; i < material_size; i++){
				glBindTexture(GL_TEXTURE_2D, _tex[i]);
				glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
				glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB,
					material_array[i].map.w,
					material_array[i].map.h,
					0, GL_RGB, GL_UNSIGNED_BYTE,
					material_array[i].map.bits);
			}
			_texIsInit = true;
		}
	}

	triangle = triangle_array;
	for( i = 0; i < triangle_size; i++, triangle++ ){
		// change texture
		if( mode & MESH_TEXTURE && image != triangle->material_index ) {
			image = triangle->material_index;
			glBindTexture( GL_TEXTURE_2D, _tex[ triangle->material_index ] );
		}
		// Setting material
		if( mode & MESH_MATERIAL && triangle->material_index >= 0 ){
			glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT, material_array[ triangle->material_index ].Ka );
			glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, material_array[ triangle->material_index ].Kd );
			glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, material_array[ triangle->material_index ].Ks );
			glMaterialf( GL_FRONT_AND_BACK, GL_SHININESS, material_array[ triangle->material_index ].Ns );
		}
		if( mode & MESH_VERTEX ) glBegin( GL_POINTS );
		else if( mode & MESH_WIRE ) glBegin( GL_LINE_LOOP );
		else glBegin( GL_TRIANGLES );
		for( j = 0; j < 3; j++ ){
			// Draw Verteces
			if( triangle->tex_coord_index[0] >= 0 && mode & MESH_TEXTURE ) {
				glTexCoord2fv( tex_coord_array + triangle->tex_coord_index[j] * 2 );
			}
			if( triangle->normal_index[0] >= 0 && mode & MESH_NORMAL ) {
				glNormal3fv( normal_array + triangle->normal_index[j] * 3 );
			}
			if( mode & MESH_MATERIAL ){
				memcpy( buf, vertex_array + triangle->vertex_index[j] * 3, 3 * sizeof( float ) );
//				buf[0] = vertex_array[ triangle_array->vertex_index[j] * 3 + 0 ];
//				buf[1] = vertex_array[ triangle_array->vertex_index[j] * 3 + 1 ];
//				buf[2] = vertex_array[ triangle_array->vertex_index[j] * 3 + 2 ];
				buf[3] = 
					( triangle->material_index < 0 ) ?
					1.0f : material_array[ triangle->material_index ].d;
				glVertex4fv( buf );
			}
			else glVertex3fv( vertex_array + triangle->vertex_index[j] * 3 );
		}
		glEnd();
	}
	if( clr_mat == GL_TRUE ) glEnable( GL_COLOR_MATERIAL );
	if( tex2d == GL_FALSE ) glDisable( GL_TEXTURE_2D );
	glShadeModel( shade );
	//glDeleteTextures( material_size, _tex );
}

void Mesh::DrawGroup( int g, int mode ){
	unsigned char tex2d, clr_mat;
	int shade;
	int i, j, image = -1;
	GLuint tex[32];
	float buf[4];
	Triangle *triangle;

	if( g >= group_size ) return;

	// mode checking
	if( vertex_size == 0 ) return;
	if( triangle_size == 0 ) mode = mode | MESH_VERTEX;
	if( normal_size == 0 ) mode = mode & ~MESH_NORMAL;
	if( tex_coord_size == 0 ) mode = mode & ~MESH_TEXTURE;
	if( material_size == 0 ) mode = mode & ~MESH_MATERIAL;

	if( mode & MESH_DEPTH_MAP ) {
		DrawDepthMap();
		return;
	}
	glGetIntegerv( GL_SHADE_MODEL, &shade );
	if( mode & MESH_SMOOTH ) glShadeModel( GL_SMOOTH );
	if( mode & MESH_FLAT ) glShadeModel( GL_FLAT );

	glGetBooleanv( GL_COLOR_MATERIAL, &clr_mat );
	glDisable( GL_COLOR_MATERIAL );

	glGetBooleanv( GL_TEXTURE_2D, &tex2d );
	if( mode & MESH_TEXTURE ){
		glEnable( GL_TEXTURE_2D );
		glGenTextures( material_size, tex );
		for( i = 0; i < material_size; i++ ){
			glBindTexture( GL_TEXTURE_2D, tex[i] );
			glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			glTexImage2D( GL_TEXTURE_2D, 0, GL_RGB,
				material_array[i].map.w,
				material_array[i].map.h,
				0, GL_RGB, GL_UNSIGNED_BYTE,
				material_array[i].map.bits );
		}
	}

	for( i = 0; i < group_array[g].triangle_size; i++ ){
		triangle = triangle_array + group_array[g].triangle_index[i];
		// change texture
		if( mode & MESH_TEXTURE && image != triangle->material_index ) {
			image = triangle->material_index;
			glBindTexture( GL_TEXTURE_2D, tex[ triangle->material_index ] );
		}
		// Setting material
		if( mode & MESH_MATERIAL && triangle->material_index >= 0 ){
			glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT, material_array[ triangle->material_index ].Ka );
			glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, material_array[ triangle->material_index ].Kd );
			glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, material_array[ triangle->material_index ].Ks );
			glMaterialf( GL_FRONT_AND_BACK, GL_SHININESS, material_array[ triangle->material_index ].Ns );
		}
		if( mode & MESH_VERTEX ) glBegin( GL_POINTS );
		else if( mode & MESH_WIRE ) glBegin( GL_LINE_LOOP );
		else glBegin( GL_TRIANGLES );
		for( j = 0; j < 3; j++ ){
			// Draw Verteces
			if( triangle->tex_coord_index[0] >= 0 && mode & MESH_TEXTURE ) {
				glTexCoord2fv( tex_coord_array + triangle->tex_coord_index[j] * 2 );
			}
			if( triangle->normal_index[0] >= 0 && mode & MESH_NORMAL ) {
				glNormal3fv( normal_array + triangle->normal_index[j] * 3 );
			}
			if( mode & MESH_MATERIAL ){
				memcpy( buf, vertex_array + triangle->vertex_index[j] * 3, 3 * sizeof( float ) );
//				buf[0] = vertex_array[ triangle->vertex_index[j] * 3 + 0 ];
//				buf[1] = vertex_array[ triangle->vertex_index[j] * 3 + 1 ];
//				buf[2] = vertex_array[ triangle->vertex_index[j] * 3 + 2 ];
				buf[3] = 
					( triangle->material_index < 0 ) ?
					1.0f : material_array[ triangle->material_index ].d;
				glVertex4fv( buf );
			}
			else glVertex3fv( vertex_array + triangle->vertex_index[j] * 3 );
		}
		glEnd();
	}
	if( clr_mat == GL_TRUE ) glEnable( GL_COLOR_MATERIAL );
	if( tex2d == GL_FALSE ) glDisable( GL_TEXTURE_2D );
	glShadeModel( shade );
	glDeleteTextures( material_size, tex );
}
#undef UseImage

void Mesh::DrawDepthMap( void ){
	unsigned char light, color;
	double m[16], max = -DBL_MAX, min = DBL_MAX, depth;
	int i, j, k;

	glGetDoublev( GL_MODELVIEW_MATRIX, m );		// get transform matrix
	glGetBooleanv( GL_COLOR_ARRAY, &color );
	glGetBooleanv( GL_LIGHTING, &light );
	// Matrix = [ [m0 m1 m2 m3]T [m4 m5 m6 m7]T [m8 m9 m10 m11]T [m12 m13 m14 m15]T ]

	// Calculate Depth Boundary
	for( i = 0; i < vertex_size; i++ ){
		depth = 0.0;
		for( j = 0; j < 3; j++ ){
			depth -= m[2 + j * 4] * ( double )vertex_array[i * 3 + j];	// depth after transform in -Z axis
		}
		depth -= m[14];
		if( depth > max ) max = depth;
		if( depth < min ) min = depth;
	}
	if( min < 0.0 ) min = 0.0;

	// Drawing
	glEnable( GL_COLOR );
	glDisable( GL_LIGHTING );
	glBegin( GL_TRIANGLES );
	for( i = 0; i < triangle_size; i++ ){
		for( j = 0; j < 3; j++ ){
			depth = 0.0;
			for( k = 0; k < 3; k++ ) {
				// depth after transform
				depth -= m[2 + k * 4] * ( double )vertex_array[ triangle_array[i].vertex_index[j] * 3 + k];		
			}
			depth -= m[14];
			depth = ( depth < 0.0 ) ? 0.0 : ( ( max - depth ) / ( max - min ) );

			glColor3d( depth, depth, depth );	// Set Depth Collor
			glVertex3fv( vertex_array + triangle_array[i].vertex_index[j] * 3 );	// Draw Vertex
		}
	}
	glEnd();
	if( light )	glEnable( GL_LIGHTING );
	if( color ) glEnable( GL_COLOR );
}

//
// Hiden Surface Mesh
//


void HidenMesh::Draw( DrawMode dump ){
	float inv[16], m[16], v[4], *pn;
	int i, j;

	glGetFloatv( GL_MODELVIEW_MATRIX, m );
	// transpose
	swap( m[1], m[4] );	swap( m[2], m[8] );	swap( m[3], m[12] );
	swap( m[6], m[9] );	swap( m[7], m[13] );
	swap( m[11], m[14] );
	invf( inv, m, 4 );	// inverse
	// transpose
	swap( inv[1], inv[4] );	swap( inv[2], inv[8] );	swap( inv[3], inv[12] );
	swap( inv[6], inv[9] );	swap( inv[7], inv[13] );
	swap( inv[11], inv[14] );

	// drawing
	glEnable( GL_BLEND );
	glBlendFunc( GL_ZERO, GL_ONE );
	glBegin( GL_TRIANGLES );
	for( i = 0; i < triangle_size; i++ ){
		memset( v, 0x00, 3 * sizeof( float ) );
		if( normal_size > 0 ){
			for( j = 0; j < 3; j++ ){
				pn = normal_array + triangle_array[i].normal_index[j] * 3;
				v[0] += pn[0];
				v[1] += pn[1];
				v[2] += pn[2];
			}
			// normal transform
			v[3] = 1.f;
			MultiVector( inv, v );
			v[2] = 1.f / v[3];
		}
		if( normal_size == 0 ||v[2] < 0.0f ){
			glVertex3fv( vertex_array + triangle_array[i].vertex_index[0] * 3 );
			glVertex3fv( vertex_array + triangle_array[i].vertex_index[1] * 3 );
			glVertex3fv( vertex_array + triangle_array[i].vertex_index[2] * 3 );
		}
	}
	glEnd();
	glBlendFunc( GL_ONE, GL_ZERO );
	glDisable( GL_BLEND );
}