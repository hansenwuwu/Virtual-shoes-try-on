#ifndef _CALIBRATION_H_
#define _CALIBRATION_H_

#include "3dconfig.h"
#include <memory\mem.h>

class CAMERA_PARAM: virtual protected Memory{
protected:
	int _n_img;			// number of images
	int _n_pt;			// number of points
	double *_H;			// homography of each image
	double _u0, _v0;	// coordinates of the principal point
	double _alpha, _beta;	// scale factors in image u and v axis
	double _gamma;		// parameter describing the skewness of the two image axis
	double _r1[3], _r2[3], _r3[3];	// the roation matrix R=[r1 r2 r3]T
	double _t[3];		// the translation vector
	double _k1, _k2;	// coefficients of distortion
	double _Rt[12];
	int _h;
	int _w;
//	int _m_size;
//	void *_mem;			// memory pool

private:
	// Solve Homography
	bool SolveHomo( Point3d* m, Point3d* M );	
	bool SolveHomo( double *m, double *M );
	// Solve Camera Calibration
	bool SolveCalib( double *m, double *M );
	// Solve Camera Distortion
	bool SolveDistor( double *m, double *M );

public:
	CAMERA_PARAM( void );
	CAMERA_PARAM( int img, int pt, int width, int height );
	CAMERA_PARAM( char *filename );
	~CAMERA_PARAM(void);
	
	// Calibration
	virtual bool Calib( Point3d *uv, Point3d *xy );
	virtual bool Calib( double *uv, double *xy );
	
	// Get 3x3 Intrinsic matrix
	virtual void GetIntrinsic( double *A, int width = 0, int height = 0 );
	virtual void GetIntrinsicInv( double *A, int width = 0, int height = 0 );
	
	// Compute Extrinsic matrix compute from intrinsic matrix
	virtual void GetExtrinsic( double *Rt );
	
	// Get 3x4 Extrinsic matrix from stored data
	virtual void GetExtrinsic( double *Rt, int image_index );
	
	// Get 3x4 Projection Perspective Matrix
	virtual void GetPPM( double *P, int image_index, int height, int width );
	virtual void GetPPM( double *P, int height, int width );

	// Setting
	virtual void SetIntrinsic( double a, double b, double r, double u0, double v0 );
	virtual void SetExtrinsic( double* m );
	virtual void SetSize( int width, int height );
	virtual void SetNum( int img_num, int pt_num );
	virtual bool WriteParam( char *filename );
	virtual bool ReadParam( char *filename );
	virtual void Clear( void );
};

#endif