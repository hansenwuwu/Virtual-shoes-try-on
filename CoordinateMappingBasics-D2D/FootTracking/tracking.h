
#pragma once

//#include <windows.h>
//#include "include\multiple_kinect_opengl.hpp"
//#include "kinect\Kinect.h"
//#include "kinect\MultiKinect.h"
#include "../lib3D/3dconfig.h"
#include "../include/scene.h"
#include "../memory/mem.h"
#include <XnTypes.h>
class   Memory;

//
// Macro
//

#define FT_BEGIN		1
#define FT_END			2

// Left/Right
#define FT_LEFT			1
#define FT_RIGHT		2

// color
#define FT_RED			0x00ff0000
#define FT_GREEN		0x0000ff00
#define FT_BLUE			0x000000ff
#define FT_SKIN			0xff000000
namespace FootTracker{

//
// Buffer Structures
//

struct _face_buffer{
	int vertex[3];
	int tex[3];
	int normal[3];
	int reserved;
};

typedef struct _mesh_buffer{
	float *v;					// vertex array
	float *vn;					// normal vector array
	float *vt;					// texture coordinate array
	struct _face_buffer *f;		// face array
	int n;						// vertex count
	int nf;						// face count
} MeshBuffer;

//
// Setting functions
//

int   InitTracker( point3f *, int width, int height );
// InitTracker:
// Initialize memory block and variables
// need to be call after you connect to KINECT
// and call it before any else tracking function

int   ExitTracker( void );
// ExitTracker:
// free all memory block and variables
// call it after all tracking function or the end of program

int   SetTrackerDepth( float * );
// SetDepth:
// set depth image to foot tracker
// it should be call every frame

int   SetTrackerTexture( color3uc * );
// SetTexture:
// set texture image to foot tracker
// it should be call every frame

int   SetTrackerBG( XnDepthPixel * );
// SetTrackerBG:
// Set depth image as back-ground image to tracker
// should be call between TrackerBGTraining( FT_BEGIN ) and TrackerBGTraining( FT_END )

int   TrackerBGTraining( int );
// TrackerBGTrain:
// Training back ground of enviroment
// param:
//	FT_BEGIN or FT_END, mean begin of training and the end of training

Point3f		DepthOfAnkle(float, float, float, float, float, float, float, float, float, float, float);

//
// Tracking functions
//

class   Tracker: protected Memory{
protected:
	float _prev[12];
	float _e;
	int _CountColorPoints( int &color, int &threshold );
	int _CompColorPoints( int &color, int &threshold, void *ptr );
	int _CountRadiusPoints( float &x, float &y, float &z, float r, int jmp = 5 );
	int _CompRadiusPoints( float &x, float &y, float &z, float r, void *, int &n, int jmp = 5 );

	void _DepthImgSegment( float &x, float &y, float &z );
public:
	int test(float m[16]);
	Tracker( void );
	int TrackInitSolution(float m[16], point3f center);
	int TrackSelect( Point3f_Array , float m[16] );
	int TrackColor( int color, int threshold, float m[16] );
	// TrackColor:
	// Tracking green marker's position and orientation,
	// param:
	//	color: the color be tracked
	//	threshold: the threshold of number of green pixels
	//	m: return the rotate and translate matrix

	float TrackICP( MeshBuffer foot, float m[16], float icpe[100], int LR);
	//int TrackICP( MeshBuffer foot, float m[16], float s[3] );
	// TrackFoot:
	// Tracking foot by a initial guess transform matrix m

	int TrackColorICP( MeshBuffer foot, int color, float m[16], float icpe[100] );
	// TrackFoot:
	// Tracking foot by a initial guess transform matrix m
	
	// Two cam
	//void GetDepthArray( float m[16], Point3f_Array *DArray );
	//float TrackMergeICP( MeshBuffer foot, float m[16], Point3f_Array *DArray1, Point3f_Array *DArray2 );


	// measurement
	float CompFootLength( int color, int threshold );
	float CompFootWidth( int color, int threshold );

	// Tracking Residual
	float Residual( void );
	// Camera Calibration
	
	void InitialCam( float mcam[16] );

};

//
// OpenGL Relating Operation
//

void   Draw( void );
void   CompGLMatf( float gl[16], float ft[16] );
void   CompGLMatd( double gl[16], double ft[16] );
void   DrawPit( int r );



};
