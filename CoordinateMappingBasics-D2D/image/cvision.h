#ifndef _CVISION_H_
#define _CVISION_H_

#include <windows.h>

//
// MeanShift Parameter
//

// Kernal Functions

#define EPANECHNIKOV_KERNAL	0X01
#define UNIFORM_KERNAL		0x02
#define NORMAL_KERNAL		0x04

// Dimension and color space

#define RGB_COLOR			0x08
#define LUV_COLOR			0x10
#define GRAY_COLOR			0x20
#define SPACE				0x40

// Feature Detector
void EdgeDetect( BITMAP &img, double threshold, bool onlyedge = true );
void CornerDetect( BITMAP &img, POINT **corner, int *n_corner, int win_size, double threshold ); 
void CornerDetect( BITMAP &img, POINT **corner, int *n_corner, int win_size, int left, int buttom, int right, int up, double threshold ); 
void sift( BITMAP &img, int oct, int s, double sig, double **xy, double **desc, int *n );
void EllipseDetect( BITMAP &img, int left, int bottom, int right, int top, 
				   float *x0, float *y0, float *a0, float *b0, float *angle, int threshold, float min = 10.0f );

// Segmentation
void KmeanClustering( BITMAP *img, int cluster, int color_space );
void KmeanClustering( BITMAP *img, int center[], int cluster, int color_space );
void KmeanSegment( BITMAP *img, int cluster, int color_space );
void KmeanSegment( BITMAP *img, int center[], int cluster, int color_space );
void KmeanSegment( BITMAP *img, int cluster, int x, int y, int color_space );
void KmeanSegment( BITMAP *img, int center[], int cluster, int x, int y, int color_space );
void MeanShift( BITMAP &img, int ***cluster_table, int *cluster, double bandwidth, int space ); 

#endif