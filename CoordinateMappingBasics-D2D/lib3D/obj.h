#ifndef _OBJ_MESH_H_
#define _OBJ_MESH_H_

#include <stdio.h>
#include <vector>
#include "memory\mem.h"
#include "3dconfig.h"

#define MESH_NONE			0x00

//
// Type Define
//

typedef int		DrawMode, GroupIndex;
typedef float	Axis_X, Axis_Y, Axis_Z;

//
// Display Mode
//

#define MESH_TRIANGLE		0x00
#define MESH_VERTEX			0x01
#define MESH_WIRE			0x02
#define MESH_SMOOTH			0x04
#define MESH_FLAT			0X08
#define MESH_NORMAL			0x10
#define MESH_MATERIAL		0x20
#define MESH_TEXTURE		0x40
#define MESH_DEPTH_MAP		0x80

//
// Output Mode
//

#define OBJ_VERTEX_ONLY		0x00
#define OBJ_USEMTL			0x01
#define OBJ_TRIANGLE		0x02
#define OBJ_TEXTURE_COORD	0x04
#define OBJ_NORMAL			0x08
#define OBJ_GROUP			0x10

//
// FaceMesh Type
//

#define FM_FACE				2
#define FM_FACE_AND_BODY	3

//
// LandMark Enumination
//

#define LM_NULL				0
#define LM_INNER_EYE_LEFT	1
#define LM_INNER_EYE_RIGHT	2
#define LM_OUTER_EYE_LEFT	3
#define LM_OUTER_EYE_RIGHT	4
#define LM_PRONASION		5
#define LM_MOUTH_LEFT		6
#define LM_MOUTH_RIGHT		7
#define LM_NASION			8
#define LM_STOMION			9
#define LM_ALA_LEFT			10
#define LM_ALA_RIGHT		11
#define LM_PROGION			12

//
// class referenced
//

//class IMAGE;
//class Mesh;

//
// Basic Structure
//

typedef struct _float3d{
	float p[3];
} Float3D;

typedef struct _float2d{
	float p[2];
} Float2D;

typedef struct _double3d{
	double p[3];
} Double3D;

typedef struct _texture{
	unsigned char *bits;
	int w;
	int h;
	int stride;
	int color;
} Texture;

typedef struct _material{
	float Ns;
	float d;
	float Ka[4];
	float Kd[4];
	float Ks[4];
	Texture map;
	char name[500];
} Material;

typedef struct _triangle{
	int vertex_index[3];
	int tex_coord_index[3];
	int normal_index[3];
	int material_index;
} Triangle;

typedef struct _group{
	int *triangle_index;
	int triangle_size;
	char name[500];
} Group;

typedef struct _neighbor{
	int *vertex_index;
	int size;
} NeighborVertex;

//
// Mesh Class
//

class Mesh: virtual protected Memory{
protected:
	float *vertex_array;
	int vertex_size, vertex_max;

	float *normal_array;
	int normal_size, normal_max;

	float *tex_coord_array;
	int tex_coord_size, tex_coord_max;

	Triangle *triangle_array;
	int triangle_size, triangle_max;

	Material *material_array;
	int material_size, material_max;

	Group *group_array;
	int group_size, group_max;
	
	NeighborVertex* neighbor_array;
	int *neighbor_vertex_array;
	int neighbor_max, neighbor_vertex_max;

	int *group_triangle_array;
	int group_triangle_size, group_triangle_max;
	
	bool _texIsInit;
	unsigned int _tex[32];
protected:
	// File Processing
	void FirstPass( FILE *, char * );
	void SecondPass( FILE * );

	// Update Processing
	void UpdateNeighbor( void );
	void UpdateTriangle( int );
	void UpdateGroup( int* );

public:
	Mesh( void );
	Mesh( const Mesh & );			// copy constructor, dosn't copy image
	~Mesh( void );
	
	// File Processing
	virtual bool ReadOBJ( char *filename );
	virtual bool ReadOBJ( std::string filename);
	virtual bool ReadMTL( char *filename );
	virtual bool WriteOBJ( char *filename, int mode );
//	virtual bool WriteMTL( char *filename );
	// VS2008 extension
	virtual bool ReadOBJ( wchar_t *filename );
	virtual bool ReadMTL( wchar_t *filename );
	virtual bool WriteOBJ( wchar_t *filename, int mode );
//	virtual bool WriteMTL( wchar_t *filename );
	
	// Mesh functions
	inline float* GetVertex( int );
	inline int GetVertex( int, float[] );
	inline Triangle* GetTriangle( int );
	inline int GetTriangle( int, Triangle * );
	inline float* GetNormal( int );
	inline int GetNormal( int, float[] );
	inline float* GetTexCoord( int );
	inline int GetTexCoord( int, float[] );
	inline Group* GetGroup( int );
	inline int GetGroup( int, Group * );
	inline int GetVertexSize( void );
	inline int GetTriangleSize( void );
	inline int GetNormalSize( void );
	inline int GetTexCoordSize( void );
	inline int GetGroupSize( void );
	int DeleteVertex( int );
	int DeleteTriangle( int );
	int DeleteTriangle( Group & );

	// Advanced functions
	std::vector< Group > ContinueTriangles( void );
	std::vector< int > NeighborTriangles( const Triangle& );
	std::vector< int > NeighborTriangles( int );
	int DeleteContinueTriangles( int );

	// Updating
	int UpdateTriangle( void );

	// Basic Function
	virtual void Clear( void );
	virtual void Scale( float, float, float );
	virtual void Translate( float, float, float );
	virtual void Rotate( float, Axis_X, Axis_Y, Axis_Z );
	virtual void CenterRotate( float, Axis_X, Axis_Y, Axis_Z );
	virtual void MultiMatrix( float m[16] );
	virtual void Smooth( void );

	// OpenGL Drawing Functions
	virtual void Draw( DrawMode );
	virtual void DrawGroup( GroupIndex, DrawMode );
	virtual void DrawDepthMap( void );

	friend class Ray;
	friend class Viewing;
	friend int Copy( Mesh &, const Mesh & );
};

//
// Hiden Surface Mesh
//

class HidenMesh: public Mesh{
public:
	virtual void Draw( DrawMode );
};

//
// Global functions
//

int Copy( Mesh &, const Mesh & );	// dosn't copy texture image

//
// Face 3D Mesh Structure
//

class FaceMesh: public Mesh{

	typedef struct _landmark{
		double inner_eye_left[3];
		double inner_eye_right[3];
		double outer_eye_left[3];
		double outer_eye_right[3];
		double nose_tip[3];
		double mouth_cor_left[3];
		double mouth_cor_right[3];
		double nasion[3];
		double lip[3];
		double ala_left[3];
		double ala_right[3];
		double progion[3];
	} LandMark;

protected:
	LandMark _lmk;

public:
	FaceMesh( void );
	~FaceMesh( void );
	
	// Basic Function
	virtual void Clean( void );
	virtual void Scale( double );
	virtual void Translate( double x, double y, double z );
	virtual void Rotate( double, double, double, double );
	virtual void MultiMatrix( double[16] );
	virtual void Smooth( void );

	// feature detect functions
	void DetectLandmark( int type = 3 );
	void DetectLandmark( float vertex[], int triangle[], int nv, int nt );
	int RefineLandmark( int lmk );
	void RefineLandmark( int lmk, float vertex[], int nv );
	int GetLandmark( int, double [] );
};

//
// inline functions
//

inline float* Mesh::GetVertex( int i ){
	if( i < 0 || i >= vertex_size ) return NULL;
	else return vertex_array + i * 3;
}
inline int Mesh::GetVertex( int i, float p[] ){
	if( i < 0 || i >= vertex_size ) return 1;
	for( int k = 0; k < 3; k++ ) p[k] = vertex_array[ i * 3 + k ];
	return 0;
}
inline float* Mesh::GetNormal( int i ){
	if( i < 0 || i >= normal_size ) return NULL;
	else return normal_array + i * 3;
}
inline int Mesh::GetNormal( int i, float p[] ){
	if( i < 0 || i >= normal_size ) return 1;
	for( int k = 0; k < 3; k++ ) p[k] = normal_array[ i * 3 + k ];
	return 0;
}
inline float* Mesh::GetTexCoord( int i ){
	if( i < 0 || i >= tex_coord_size ) return NULL;
	else return tex_coord_array + i * 2;
}
inline int Mesh::GetTexCoord( int i, float p[] ){
	if( i < 0 || i >= tex_coord_size ) return 1;
	for( int k = 0; k < 3; k++ ) p[k] = tex_coord_array[ i * 2 + k ];
	return 0;
}
inline Triangle* Mesh::GetTriangle( int i ){
	if( i < 0 || i >= triangle_size ) return NULL;
	else return triangle_array + i;
}
inline int Mesh::GetTriangle( int i, Triangle *p ){
	if( i < 0 || i >= triangle_size ) return 1;
	for( int k = 0; k < 3; k++ ) {
		p->vertex_index[k] = triangle_array[i].vertex_index[k];
		p->normal_index[k] = triangle_array[i].normal_index[k];
		p->tex_coord_index[k] = triangle_array[i].tex_coord_index[k];
	}
	p->material_index = triangle_array[i].material_index;
	return 0;
}

inline Group* Mesh::GetGroup( int i ){
	if( i < 0 || i >= group_size ) return NULL;
	else return group_array + i;
}

inline int Mesh::GetGroup( int i, Group *p ){
	if( i < 0 || i >= group_size ) return 1;
	memcpy( p, group_array + i, sizeof( Group ) );
	return 0;
}

inline int Mesh::GetVertexSize( void ){
	return vertex_size;
}

inline int Mesh::GetNormalSize( void ){
	return normal_size;
}

inline int Mesh::GetTriangleSize( void ){
	return triangle_size;
}

inline int Mesh::GetTexCoordSize( void ){
	return tex_coord_size;
}

inline int Mesh::GetGroupSize( void ){
	return group_size;
}

#endif