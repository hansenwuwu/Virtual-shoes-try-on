
#ifndef _3D_TRACK_H_
#define _3D_TRACK_H_

#include "3dconfig.h"

#ifdef __cplusplus
extern "C"{
#endif

//
// structures
//

// transform paramter (float)
typedef struct _transform_parameter_float{
	int iter;
	float err;
	float lr1;
	float lr2;
} Transform_Paramf;

// transform paramter (double)
typedef struct _transform_parameter_double{
	int iter;
	double err;
	double lr1;
	double lr2;
} Transform_Paramd;

// ICP paramter (float)
typedef struct _icp_parameter_float{
	int iter;
	struct _transform_parameter_float trans;
} ICP_Paramf;

// ICP paramter (double)
typedef struct _icp_parameter_double{
	int iter;
	struct _transform_parameter_double trans;
} ICP_Paramd;

// Motion vector (float)
typedef struct _motion_vector_float{
	Point3f *p;
	Vector3f *v;
	int n;
} MotionVecf;

// Motion vector (double)
typedef struct _motion_vector_double{
	Point3d *p;
	Vector3d *v;
	int n;
} MotionVecd;

//
// Data Translation
//

// Solve paired data, by gradient descent
float RigidTransformf( Point3f_Pair, RT_Paramf *, const Transform_Paramf );
double RigidTransformd( Point3d_Pair, RT_Paramd *, const Transform_Paramd );
float RigidTransformfv( Point3fv_Pair, RT_Paramf *, const Transform_Paramf );
double RigidTransformdv( Point3dv_Pair, RT_Paramd *, const Transform_Paramd );

float RigidTransformMatf( Point3f_Pair, float [3][4], const Transform_Paramf );
double RigidTransformMatd( Point3d_Pair, double [3][4], const Transform_Paramd );
float RigidTransformMatfv( Point3fv_Pair, float [3][4], const Transform_Paramf );
double RigidTransformMatdv( Point3dv_Pair, double [3][4], const Transform_Paramd );

// Solve un-paired data, by Iterative Closest Point( ICP )
float ICPf( Point3f_Array, Point3f_Array, RT_Paramf *, const ICP_Paramf, const Point3f_Array, float icpe[100]  );
double ICPd( Point3d_Array, Point3d_Array, RT_Paramd *, const ICP_Paramd, const Point3d_Array );
float ICPfv( Point3fv_Array, Point3fv_Array, RT_Paramf *, const ICP_Paramf, const Point3fv_Array );
double ICPdv( Point3dv_Array, Point3dv_Array, RT_Paramd *, const ICP_Paramd, const Point3dv_Array );
float TrICPf( Point3f_Array, Point3f_Array, RT_Paramf *, const ICP_Paramf, const Point3f_Array, float icpe[100] );
float ICP_kdtree(Point3f_Array, Point3f_Array, RT_Paramf *, const ICP_Paramf, const Point3f_Array, float icpe[100]);

// Solve cloud data, by PCA 
void CompPCACoordf( Point3f_Array, float m[4][4] );
void CompPCACoordd( Point3d_Array, double m[4][4] );

// Rigid Transform by SVD, p1 * rt = p2
void RigidTransform_SVDMatf( Point3f_Pair, float [3][4] );
void RigidTransform_SVDMatfv( Point3fv_Pair, float [3][4] );


#ifdef __cplusplus
}
#endif

#endif