#include "../image/imgproc.h"

//
// imgproc2.cpp
//

#include "imgconf.h"

//
// Image Painting
//

int DrawLine1b( ImageBuffer buf, ImagePoint, ImagePoint, Color1b );
int DrawLine3b( ImageBuffer buf, ImagePoint, ImagePoint, Color3b );
int DrawLine4b( ImageBuffer buf, ImagePoint, ImagePoint, Color4b );
int DrawLine1i( ImageBuffer buf, ImagePoint, ImagePoint, Color1i );
int DrawLine3i( ImageBuffer buf, ImagePoint, ImagePoint, Color3i );
int DrawLine4i( ImageBuffer buf, ImagePoint, ImagePoint, Color4i );
int DrawLine1f( ImageBuffer buf, ImagePoint, ImagePoint, Color1f );
int DrawLine3f( ImageBuffer buf, ImagePoint, ImagePoint, Color3f );
int DrawLine4f( ImageBuffer buf, ImagePoint, ImagePoint, Color4f );
int DrawLine1d( ImageBuffer buf, ImagePoint, ImagePoint, Color1d );
int DrawLine3d( ImageBuffer buf, ImagePoint, ImagePoint, Color3d );
int DrawLine4d( ImageBuffer buf, ImagePoint, ImagePoint, Color4d );

int DrawTriangle1b( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color1b, Color1b, Color1b );
int DrawTriangle3b( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color3b, Color3b, Color3b );
int DrawTriangle4b( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color4b, Color4b, Color4b );
int DrawTriangle1i( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color1i, Color1i, Color1i );
int DrawTriangle3i( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color3i, Color3i, Color3i );
int DrawTriangle4i( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color4i, Color4i, Color4i );
int DrawTriangle1f( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color1f, Color1f, Color1f );
int DrawTriangle3f( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color3f, Color3f, Color3f );
int DrawTriangle4f( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color4f, Color4f, Color4f );
int DrawTriangle1d( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color1d, Color1d, Color1d );
int DrawTriangle3d( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color3d, Color3d, Color3d );
int DrawTriangle4d( ImageBuffer, ImagePoint, ImagePoint, ImagePoint, Color4d, Color4d, Color4d );

int Filling1b( ImageBuffer, Region, Color1b );
int Filling3b( ImageBuffer, Region, Color3b );
int Filling4b( ImageBuffer, Region, Color4b );
int Filling1i( ImageBuffer, Region, Color1i );
int Filling3i( ImageBuffer, Region, Color3i );
int Filling4i( ImageBuffer, Region, Color4i );
int Filling1f( ImageBuffer, Region, Color1f );
int Filling3f( ImageBuffer, Region, Color3f );
int Filling4f( ImageBuffer, Region, Color4f );
int Filling1d( ImageBuffer, Region, Color1d );
int Filling3d( ImageBuffer, Region, Color3d );
int Filling4d( ImageBuffer, Region, Color4d );

//
// Editing
//

int Selectb( ImageBuffer, Selected, ImagePoint, float threshold );
int Selecti( ImageBuffer, Selected, ImagePoint, float threshold );
int Selectf( ImageBuffer, Selected, ImagePoint, float threshold );
int Selectd( ImageBuffer, Selected, ImagePoint, float threshold );

//
// Image Processing
//

int Convolutionb( ImageBuffer buf, int win, float *kernel );
int Convolutioni( ImageBuffer buf, int win, float *kernel );
int Convolutionf( ImageBuffer buf, int win, float *kernel );
int Convolutiond( ImageBuffer buf, int win, float *kernel );

int GaussFilterb( ImageBuffer buf, int win );
int GaussFilteri( ImageBuffer buf, int win );
int GaussFilterf( ImageBuffer buf, int win );
int GaussFilterd( ImageBuffer buf, int win );

int Scaleb( ImageBuffer src, ImageBuffer tar );
int Scalei( ImageBuffer src, ImageBuffer tar );
int Scalef( ImageBuffer src, ImageBuffer tar );
int Scaled( ImageBuffer src, ImageBuffer tar );

void Ditherb( ImageBuffer, unsigned char * );
void Ditheri( ImageBuffer, int *, int max, int min );
void Ditherf( ImageBuffer, float *, float max, float min );
void Ditherd( ImageBuffer, double *, double max, double min );

//
// Color Converting
//

int RGBtoNTSCi( ImageBuffer );
int RGBtoNTSCf( ImageBuffer );
int RGBtoNTSCd( ImageBuffer );

int NTSCtoRGBi( ImageBuffer );
int NTSCtoRGBf( ImageBuffer );
int NTSCtoRGBd( ImageBuffer );
