#ifndef _IMAGE_CONFIG_H_
#define _IMAGE_CONFIG_H_

#ifdef __cplusplus
extern "C"{
#endif

//
// Image Buffer
//

typedef struct _image_buffer{
	void *img;
	int w;
	int h;
	int stride;
	int color;
} ImageBuffer, GradientMap;

typedef struct _region_buffer{
	char *img;
	int w;
	int h;
} Region, Selected;

//
// Coordinate
//

typedef struct _img_point{
	int x, y;
} ImagePoint;

//
// Color or Pixel
//

typedef struct _pixel_byte_1d{
	unsigned char x;
} ImageByte1, Color1b;

typedef struct _pixel_byte_3d{
	unsigned char r, g, b;
} ImageByte3, Color3b;

typedef struct _pixel_byte_4d{
	unsigned char r, g, b, a;
} ImageByte4, Color4b;

typedef struct _pixel_int_1d{
	int x;
} ImageInt1, Color1i;

typedef struct _pixel_int_3d{
	int r, g, b;
} ImageInt3, Color3i;

typedef struct _pixel_int_4d{
	int r, g, b, a;
} ImageInt4, Color4i;

typedef struct _pixel_float_1d{
	float x;
} ImageFloat1, Color1f;

typedef struct _pixel_float_3d{
	float r, g, b;
} ImageFloat3, Color3f;

typedef struct _pixel_float_4d{
	float r, g, b, a;
} ImageFloat4, Color4f;

typedef struct _pixel_double_1d{
	double x;
} ImageDouble1, Color1d;

typedef struct _pixel_double_3d{
	double r, g, b;
} ImageDouble3, Color3d;

typedef struct _pixel_double_4d{
	double r, g, b, a;
} ImageDouble4, Color4d;

#ifdef __cplusplus
}
#endif

#endif