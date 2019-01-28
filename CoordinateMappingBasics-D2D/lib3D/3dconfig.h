#ifndef _3D_CONFIG_H_
#define _3D_CONFIG_H_

#ifdef __cplusplus
extern "C"{
#endif

//
// structures
//

// 2D Point

// Point 2D (int)
typedef struct _point_2d_int{
	int x, y;
} Point2i, Vector2i;

// Point 2D (float)
typedef struct _point_2d_float{
	float x, y;
} Point2f, Vector2f;

// Point 2D (double)
typedef struct _point_2d_double{
	double x, y;
} Point2d, Vector2d;

// Point 2D in pointer (int)
typedef struct _point_2d_int_ptr{
	int *ptr;
} Point2iv, Vector2iv;

// Point 2D in pointer (float)
typedef struct _point_2d_float_ptr{
	float *ptr;
} Point2fv, Vector2fv;

// Point 2D in pointer (double)
typedef struct _point_2d_double_ptr{
	double *ptr;
} Point2dv, Vector2dv;

// Array of Points 2D (int)
typedef struct _point_2d_int_array{
	Point2i *pt;
	int n;
} Point2i_Array, Vector2i_Array;

// Array of Points 2D (float)
typedef struct _point_2d_float_array{
	Point2f *pt;
	int n;
} Point2f_Array, Vector2f_Array;

// Array of Points 2D (double)
typedef struct _point_2d_double_array{
	Point2d *pt;
	int n;
} Point2d_Array, Vector2d_Array;

// Array of Points 2D in pointer (int)
typedef struct _point_2d_int_array_ptr{
	Point2iv *pt;
	int n;
} Point2iv_Array, Vector2iv_Array;

// Array of Points 2D in pointer (float)
typedef struct _point_2d_float_array_ptr{
	Point2fv *pt;
	int n;
} Point2fv_Array, Vector2fv_Array;

// Array of Points 2D in pointer (double)
typedef struct _point_2d_double_array_ptr{
	Point2dv *pt;
	int n;
} Point2dv_Array, Vector2dv_Array;

// 3D Point

// Point 3D (int)
typedef struct _point_3d_int{
	int x, y, z;
} Point3i, Vector3i;

// Point 3D (float)
typedef struct _point_3d_float{
	float x, y, z;
} Point3f, Vector3f;

// Point 3D (double)
typedef struct _point_3d_double{
	double x, y, z;
} Point3d, Vector3d;

// Point 3D in pointer (int)
typedef struct _point_3d_int_ptr{
	int *ptr;
} Point3iv, Vector3iv;

// Point 3D in pointer (float)
typedef struct _point_3d_float_ptr{
	float *ptr;
} Point3fv, Vector3fv;

// Point 3D in pointer (double)
typedef struct _point_3d_double_ptr{
	double *ptr;
} Point3dv, Vector3dv;

// Array of Points 3D (int)
typedef struct _point_3d_int_array{
	Point3i *pt;
	int n;
} Point3i_Array, Vector3i_Array;

// Array of Points 3D (float)
typedef struct _point_3d_float_array{
	Point3f *pt;
	int n;
} Point3f_Array, Vector3f_Array;

// Array of Points 3D (double)
typedef struct _point_3d_double_array{
	Point3d *pt;
	int n;
} Point3d_Array, Vector3d_Array;

// Array of Points 3D in pointer (int)
typedef struct _point_3d_int_array_ptr{
	Point3iv *pt;
	int n;
} Point3iv_Array, Vector3iv_Array;

// Array of Points 3D in pointer (float)
typedef struct _point_3d_float_array_ptr{
	Point3fv *pt;
	int n;
} Point3fv_Array, Vector3fv_Array;

// Array of Points 3D in pointer (double)
typedef struct _point_3d_double_array_ptr{
	Point3dv *pt;
	int n;
} Point3dv_Array, Vector3dv_Array;

// Array of Paired 3D Points (int)
typedef struct _point_3d_int_pair_array{
	Point3i *pt1, *pt2;
	int n;
} Point3i_Pair;

// Array of Paired 3D Points (float)
typedef struct _point_3d_float_pair_array{
	Point3f *pt1, *pt2;
	int n;
} Point3f_Pair;

// Array of Paired 3D Points (double)
typedef struct _point_3d_double_pair_array{
	Point3d *pt1, *pt2;
	int n;
} Point3d_Pair;

// Array of Paired 3D Points in pointer (int)
typedef struct _point_3d_int_ptr_pair_array{
	Point3iv *pt1, *pt2;
	int n;
} Point3iv_Pair;

// Array of Paired 3D Points in pointer (float)
typedef struct _point_3d_float_ptr_pair_array{
	Point3fv *pt1, *pt2;
	int n;
} Point3fv_Pair;

// Array of Paired 3D Points in pointer (double)
typedef struct _point_3d_double_ptr_pair_array{
	Point3dv *pt1, *pt2;
	int n;
} Point3dv_Pair;

//
// Rotate and Translate
//

// rotate and translate parameter (int)
typedef struct _rotate_translate_param_int{
	int x, y, z, a, b, c;
} RT_Parami, RP_Parami;

// rotate and translate parameter (float)
typedef struct _rotate_translate_param_float{
	float x, y, z, a, b, c;
} RT_Paramf, RP_Paramf;

// rotate and translate parameter (double)
typedef struct _rotate_translate_param_double{
	double x, y, z, a, b, c;
} RT_Paramd, RP_Paramd;

// axis-rotation parameter (int)
typedef struct _rotate_axis_param_int{
	Point3i p;
	Vector3i v;
	int a;
} RA_Parami;

// axis-rotation parameter (float)
typedef struct _rotate_axis_param_float{
	Point3f p;
	Vector3f v;
	float a;
} RA_Paramf;

// axis-rotation parameter (double)
typedef struct _rotate_axis_param_double{
	Point3d p;
	Vector3d v;
	double a;
} RA_Paramd;

//
// Basic Linear Algebra
//

// Dot: a * b

int Doti( Point3i a, Point3i b );
float Dotf( Point3f a, Point3f b );
double Dotd( Point3d a, Point3d b );
int Dotiptr( Point3iv a, Point3iv b );
float Dotfptr( Point3fv a, Point3fv b );
double Dotdptr( Point3dv a, Point3dv b );

// Cross: c = a x b

void Crossi( Point3i a, Point3i b, Point3i *c );
void Crossf( Point3f a, Point3f b, Point3f *c );
void Crossd( Point3d a, Point3d b, Point3d *c );
void Crossiv( Point3iv a, Point3iv b, Point3iv c );
void Crossfv( Point3fv a, Point3fv b, Point3fv c );
void Crossdv( Point3dv a, Point3dv b, Point3dv c );

//
// Functions decleration
//

// Solve Rotate and Translation

void CompRTMati( float [3][4], RT_Parami );
void CompRTMatf( float [3][4], RT_Paramf );
void CompRTMatd( double [3][4], RT_Paramd );

void CompRPMati( float [3][4], RP_Parami );
void CompRPMatf( float [3][4], RP_Paramf );
void CompRPMatd( double [3][4], RP_Paramd );

void CompRAMati( float [3][4], RA_Parami );
void CompRAMatf( float [3][4], RA_Paramf );
void CompRAMatd( double [3][4], RA_Paramd );

void CompRTFromMati( RT_Parami *, float [3][4] );
void CompRTFromMatf( RT_Paramf *, float [3][4] );
void CompRTFromMatd( RT_Paramd *, double [3][4] );

//
// Transformation
//

void MultiMatf( float [3][4], const float [3][4], const float [3][4] );
void MultiMatd( double [3][4], const double [3][4], const double [3][4] );

void MultiPointf( Point3f *, float [3][4] );
void MultiPointd( Point3d *, double [3][4] );
void MultiPointfv( Point3fv *, float [3][4] );
void MultiPointdv( Point3dv *, double [3][4] );

void TransformPointf( Point3f *, RT_Paramf );
void TransformPointd( Point3d *, RT_Paramd );
void TransformPointfv( Point3fv *, RT_Paramf );
void TransformPointdv( Point3dv *, RT_Paramd );

// Solve inverse of transform

void InvRTMatf( float [3][4], const float [3][4] );
void InvRTMatd( double [3][4], const double [3][4] );

//Kai
void Mult16Matf(float out[16],float A[16],float B[16]);
void Trans16Matf(float m[16]);


#ifdef __cplusplus
}
#endif

#endif