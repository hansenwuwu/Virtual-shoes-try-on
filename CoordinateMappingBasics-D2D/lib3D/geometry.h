
#ifndef _PARAMETRIC_GEOMETRY_H_
#define _PARAMETRIC_GEOMETRY_H_

#include "3dconfig.h"

#ifdef __cplusplus
extern "C"{
#endif

//
// Line Definition
//

typedef struct _line_3d_float{
	Point3f p, v;	// start point and direction
} Line3f;

typedef struct _line_3d_double{
	Point3d p, v;	// start point and direction
} Line3d;

Point3f LinePointf( Line3f, float );
Point3d LinePointd( Line3d, double );

//
// Cylinder Definition
//

typedef struct _cylinder_float{
	Point3f p;		// center point of bottom
	Point3f v;		// direction of cylinder
	float r, h;		// radius and height
} Cylinderf;

typedef struct _cylinder_double{
	Point3d p;		// center point of bottom
	Point3d v;		// direction of cylinder
	double r, h;	// radius and height
} Cylinderd;

Point3f CylinderPointf( Cylinderf, float angle, float height );
Point3d CylinderPointd( Cylinderd, double angle, double height );

//
// intersection
//

Point3f LineCylinderIntersect( Line3f, Cylinderf );


#ifdef __cplusplus
}
#endif

#endif
