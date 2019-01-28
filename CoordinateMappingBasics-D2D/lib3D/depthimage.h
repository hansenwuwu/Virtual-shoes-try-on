#ifndef _DEPTH_IMAGE_H_
#define _DEPTH_IMAGE_H_

#include <image\image.h>
#include "3dconfig.h"
#include "surface.h"
#include "calibration.h"

//
// Structures
//
/*
typedef struct _curvature_map_buffer{
	double *k1, *k2;
	int w, h;
} CurvatureBuf;
/**/
//
// Depth Image Class
//

class DepthImage: public ImageR< int >, virtual public CAMERA_PARAM{
protected:
	// some addition data...
	
public:
	DepthImage( void );

	// Basic Functions
	virtual int Initialize( int w, int h, int color );
	virtual void Clear( void );
	
	// File Processing
	virtual int ReadFile( const char * );
	virtual int ReadBMP( const char * );
	virtual int ReadJPEG( const char * );
	virtual int ReadPPM( const char * );

	// Editting Functions	
	virtual int SetImage( int *bits, int w, int h, int color );	// assign the image only

	// 3D information
	virtual Point3f Comp3DPosf( int x, int y );
	virtual Point3d Comp3DPosd( int x, int y );
	virtual Point3f CompNormalf( int x, int y );
	virtual Point3d CompNormald( int x, int y );
	virtual int Curvature( double *k1, double *k2, int x, int y, int r );	// principle curvature
	virtual double GaussCurvature( int x, int y );
	virtual double MeanCurvature( int x, int y );
	virtual int CompRotateTranslate( Point2i_Array pv, double Rt[16] );	// compute local coordinate
};

#endif