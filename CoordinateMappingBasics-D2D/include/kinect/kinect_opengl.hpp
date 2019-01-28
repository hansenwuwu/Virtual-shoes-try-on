
#ifndef __KINECT_OPENGL__
#define __KINECT_OPENGL__

#include <scene.h>
#include <gl/glew.h>
#include "kinect.hpp"


/* reduce time costing when 3D display, must be 1, 2, 4, ..., 
   the higher value the faster speed but lower clearness */
#define			_JUMP_SCROLL_			2
static int		JUMP					= 1;


static point2i	DISPLAY_SIZE			;		// display size
static point2f	WORLD_TO_WINDOW			= point2f( 1.f, 1.f );
static point2f	WINDOW_TO_WORLD			= point2f( 1.f, 1.f );
static color3f	DEPTH_RGB				;
static color3f	DIST_RGB				;


#ifndef		__R_PLUS_PLUS__
	#define		__R_PLUS_PLUS__
	#define		SQR( x )	( (x)*(x) )
	#define		ERR			0.0000000000001
#endif

int				K_LUfact				( double *A, int *index, int n, double *mem );
void			K_SolveLU				( double *A, int *ind, double *b, int n );
void			K_LUcomb				( double *A, int *ind, int n );
double			K_InvIter				( double *A, double *b0, double shift, int n, int time, double *memory );

bool			K_SetDisplaySize		( int width, int height );

static KINECT*	K_CAM					= NULL;
bool			K_SetDevice				( KINECT *device );
bool			K_Resize				( int width = DISPLAY_SIZE.x, int height = DISPLAY_SIZE.y );
bool			K_InitGLView			( void );
bool			K_SetDepthColor			( float in_red, float in_green, float in_blue );

void			K_DrawGray				( const float *in_gray );
void			K_DrawBackground		( const color *in_color );
void			K_DrawWithAlpha			( const color4f *in_color4f );
void			K_GetDepthBackground	( const point *in_point, color *out_color );
void			K_GetAlphaBackground	( const color *in_color, color4f *out_color4f, float in_alpha );
void			K_DrawObjectDepth		( const point *in_point, const float in_alpha );
void			K_DrawObjectColor		( const point *in_point, const color *in_color, float in_alpha );
void			K_DrawObjectPixel		( const point *in_point, const color *in_color, float in_alpha );


/*===============================================================
	Do
===============================================================*/

bool K_SetDevice( KINECT *_device )
{
	if( _device == NULL )
		return 0;

	K_CAM = _device;
	K_Resize();
	return 1;
}

bool K_Resize( int _width, int _height )
{
	if( _width <= 0 || _height <= 0 )
		return 0;

	DISPLAY_SIZE = point2i( _width, _height );
	WORLD_TO_WINDOW.x = (float)DISPLAY_SIZE.x/K_CAM->width;
	WORLD_TO_WINDOW.y = (float)DISPLAY_SIZE.y/K_CAM->height;
	WINDOW_TO_WORLD.x = (float)K_CAM->width/DISPLAY_SIZE.x;
	WINDOW_TO_WORLD.y = (float)K_CAM->height/DISPLAY_SIZE.y;
	return 1;
}

bool K_SetDepthColor( float _red, float _green, float _blue )
{
	// check color
	if( ( _red < 0.f || _red > 1.f ) || ( _green < 0.f || _green > 1.f ) || ( _blue < 0.f || _blue > 1.f ) )
	{
		DEPTH_RGB = .5f;
		DIST_RGB = .5f/K_CAM->range;
		return 1;
	}

	// set depth color
	DEPTH_RGB = color( _red, _green, _blue );
	DIST_RGB = DEPTH_RGB/K_CAM->range;
	return 0;
}

void K_DrawGray( const float *_gray )
{
	glDisable( GL_DEPTH_TEST );
		glPixelZoom( (float)WORLD_TO_WINDOW.x, (float)-WORLD_TO_WINDOW.y );
		glRasterPos2f( -1.f, 1.f );
		glPixelStorei( GL_PACK_ALIGNMENT, 4 );
		glDrawPixels( K_CAM->width, K_CAM->height, GL_LUMINANCE, GL_FLOAT, _gray );
	glEnable( GL_DEPTH_TEST );
}

void K_DrawWithAlpha( const color4f *_color4f )
{
	glDisable( GL_DEPTH_TEST );
		glPixelZoom( (float)WORLD_TO_WINDOW.x, (float)-WORLD_TO_WINDOW.y );
		glRasterPos2f( -1.f, 1.f );
		glPixelStorei( GL_PACK_ALIGNMENT, 4 );
		glDrawPixels( K_CAM->width, K_CAM->height, GL_RGBA, GL_FLOAT, _color4f );
	glEnable( GL_DEPTH_TEST );
}

void K_DrawBackground( const color *_color )
{
	glDisable( GL_DEPTH_TEST );
		glPixelZoom( (float)WORLD_TO_WINDOW.x, (float)-WORLD_TO_WINDOW.y );
		glRasterPos2f( -1.f, 1.f );
		glPixelStorei( GL_PACK_ALIGNMENT, 4 );
		glDrawPixels( K_CAM->width, K_CAM->height, GL_RGB, GL_FLOAT, _color );
	glEnable( GL_DEPTH_TEST );
}

void K_GetDepthBackground( const point *_point, color *_color )
{
	for( int i = 0; i < K_CAM->size; i++, _point++, _color++ )
	{
		if( _point[0].z )
			_color[0] = DIST_RGB*( K_CAM->range + _point[0].z );
		else
			memset( _color, 0, sizeof(color) );
	}
}

void K_GetAlphaBackground( const color *_color, color4f *_color4f, float _alpha )
{
	for( int i = 0; i < K_CAM->size; i++, _color++, _color4f++ )
		if( _color[0].r || _color[0].g || _color[0].b )
			_color4f[0] = color4f( _color[0].r, _color[0].g, _color[0].b, _alpha );
		else
			memset( _color4f, 0, sizeof(color4f) );
}

void K_DrawObjectColor( const point *_point, const color *_color, float _alpha )
{
	const color* ptr_c = _color + K_CAM->width;
	const point* ptr_p = _point + K_CAM->width;

	glBegin( GL_TRIANGLE_STRIP );
	{
		for( int i = 1, j = 1; i < K_CAM->size; i += JUMP, j += JUMP, _point += JUMP, ptr_p += JUMP, _color += JUMP, ptr_c += JUMP )
		{
			if( j >= K_CAM->width )
			{
				j = 0;
				glEnd();
				glBegin( GL_TRIANGLE_STRIP );
				continue;
			}

			if( _point[0].z )
			{
				glColor4f( _color[0].r, _color[0].g, _color[0].b, _alpha );
				glVertex3f( _point[0].x, _point[0].y, _point[0].z );
			}
			else
			{
				glEnd();
				glBegin( GL_TRIANGLE_STRIP );
				continue;
			}

			if( ptr_p[0].z )
			{
				glColor4f( ptr_c[0].r, ptr_c[0].g, ptr_c[0].b, _alpha );
				glVertex3f( ptr_p[0].x, ptr_p[0].y, ptr_p[0].z );
			}
			else
			{
				glEnd();
				glBegin( GL_TRIANGLE_STRIP );
				continue;
			}
		}
	}
	glEnd();
}

void K_DrawObjectDepth( const point *_point, float _alpha )
{
	float temp;

	const point *ptr_p = _point + K_CAM->width;
	glBegin( GL_TRIANGLE_STRIP );
	{
		for( int i = 1, j = 1; i < K_CAM->size; i += JUMP, j += JUMP, _point += JUMP, ptr_p += JUMP )
		{
			if( j >= K_CAM->width )
			{
				j = 0;
				glEnd();
				glBegin( GL_TRIANGLE_STRIP );
				continue;
			}

			if( _point[0].z )
			{
				temp = K_CAM->range + _point[0].z;
				glColor4f( DIST_RGB.b*temp, DIST_RGB.g*temp, DIST_RGB.r*temp, _alpha );
				glVertex3f( _point[0].x, _point[0].y, _point[0].z );
			}
			else
			{
				glEnd();
				glBegin( GL_TRIANGLE_STRIP );
				continue;
			}

			if( ptr_p[0].z )
			{
				temp = K_CAM->range + ptr_p[0].z;
				glColor4f( DIST_RGB.b*temp, DIST_RGB.g*temp, DIST_RGB.r*temp, _alpha );
				glVertex3f( ptr_p[0].x, ptr_p[0].y, ptr_p[0].z );
			}
			else
			{
				glEnd();
				glBegin( GL_TRIANGLE_STRIP );
				continue;
			}
		}
	}
	glEnd();
}

void K_DrawObjectPixel( const point *_point, const color *_color, float _alpha )
{
	glBegin( GL_POINTS );
		for( int i = 0; i < K_CAM->size; i++, _point++, _color++ )
			if( _point[0].z )
			{
				glColor4f( _color[0].r, _color[0].g, _color[0].b, _alpha );
				glVertex3f( _point[0].x, _point[0].y, _point[0].z );
			}
	glEnd();
}

bool K_InitGLView( void )
{
	if( !K_SetDevice( K_CAM ) )
		return false;

	double *a, *pa1, *pa2, *b, *pb1, *pb2, temp;
	void *mem;
	int i, j, k, n;

	// memory pool
	mem = malloc( 
		12*K_CAM->width*K_CAM->height*sizeof(double) + 
		6*6*sizeof(double) + 
		6*( sizeof(int) + sizeof(double) ) +
		K_CAM->size*sizeof(point)
	);
	a = (double*)( (int*)( (double*)mem + 6 ) + 6 );	// width * height * 12
	b = a + K_CAM->width*K_CAM->height*12;				// 6 * 6
	if( mem == NULL )
		return 1;

	point *map = (point*)( b + 6*6 );
	point *ptrm = map;

	const XnDepthPixel *XN_DEPTH = K_CAM->depthGener.GetDepthMap();
	for( i = 0; i < K_CAM->height; i++ )
	{
		for( j = 0; j < K_CAM->width; j++, ptrm++, XN_DEPTH++ )
		{
			ptrm[0].x = (float)j;
			ptrm[0].y = (float)i;
			if( XN_DEPTH[0] > K_CAM->range )
				ptrm[0].z = 0.f;
			else
				ptrm[0].z = (float)-XN_DEPTH[0];
		}
	}
	K_CAM->depthGener.ConvertProjectiveToRealWorld( K_CAM->size, (XnPoint3D*)map, (XnPoint3D*)map );

	// construct matrix A
	n = 0;								// count the meaningful data
	pa1 = a;	
	pa2 = a + 6;
	for( i = 0; i < K_CAM->height; i++ )
	{
		for( j = 0; j < K_CAM->width; j++, map++ )
		{
			if( -map[0].z < FLT_MIN )
				continue;

			// row 1
			pa1[0] = -map[0].x;
			pa1[1] = map[0].y;
			pa1[2] = -map[0].z;
			pa1[3] = pa1[4] = 0.0;
			pa1[5] = map[0].z*(double)j;
			// row 2
			pa2[0] = pa2[1] = pa2[2] = 0.0;
			pa2[3] = map[0].y;
			pa2[4] = -map[0].z;
			pa2[5] = map[0].z*(double)i;

			pa1 += 12;
			pa2 += 12;
			n += 2;
		}
	}

	// B = AT * A
	for( i = 0; i < 6; i++ )
	{
		pb1 = pb2 = b + i*6 + i;
		for( j = i; j < 6; j++, pb1++, pb2 += 6 )
		{
			pa1 = a + i;
			pa2 = a + j;
			temp = 0.0;
			for( k = 0; k < n; k++, pa1 += 6, pa2 += 6 )
				temp += pa1[0]*pa2[0];
			pb1[0] = pb2[0] = temp;
		}
	}
	
	// Find Minimal eigenvector
	K_InvIter( b, a, 0.0, 6, 30, (double*)mem );

	// Set projection
	for( i = 0; i < 5; i++ )
		a[i] *= 1.0/a[5];

	// Scale by width
	temp = (double)DISPLAY_SIZE.x/K_CAM->width;
	for( i = 0; i < 3; i++ )
		a[i] *= temp;

	// Scale by height
	temp = (double)DISPLAY_SIZE.y/K_CAM->height;
	for( ; i < 5; i++ )
		a[i] *= temp;
	        
	K_CAM->nearClip = 1.0;
	K_CAM->farClip = 10000.0;
	
	double Fw = K_CAM->nearClip*DISPLAY_SIZE.x/a[0];
	double Fh = K_CAM->nearClip*DISPLAY_SIZE.y/a[3];
	K_CAM->view_l = Fw*( -a[2]/DISPLAY_SIZE.x );
	K_CAM->view_r = Fw*( 1 - a[2]/DISPLAY_SIZE.x );
	K_CAM->view_t = Fh*( a[4]/DISPLAY_SIZE.y );
	K_CAM->view_b = Fh*( -1 + a[4]/DISPLAY_SIZE.y );

	// free memory
	free( mem );

	return 0;
}

int K_LUfact( double *A, int *index, int n, double *mem )
{
	double sum, big, *v, *a, *b;
	int i, j, k, row;
	void *memory;

	// memory pool
	if( mem == NULL ) memory = malloc( n * sizeof( double ) );
	else memory = mem;
	v = ( double * )memory;

	// calculate largest value per row
	a = A;
	for( i = 0; i < n; i++ ){
		big = 0.0;
		for( j = 0; j < n; j++, a++ ){
			if( ( sum = fabs( *a ) ) > big ) big = sum;
		}
		if( big == 0.0 ) return 1;
		v[i] = 1.0 / big;
	}
	for( j = 0; j < n; j++ ){
		// calculate U when i < j
		for( i = 0; i < j; i++ ){
			sum = A[ i * n + j ];
			for( k = 0; k < i; k++ ) sum -= A[ i * n + k ] * A[ k * n + j ];
			A[ i * n + j ] = sum;
		}
		big = 0.0;
		// calculate L wneh i >= j;
		for( i = j; i < n; i++ ){
			sum = A[ i * n + j ];
			for( k = 0; k < j; k++ ) sum -= A[ i * n + k ] * A[ k * n + j ];
			A[ i * n + j ] = sum;
			sum = v[i] * fabs( sum ) ;
			if( sum >= big ){
				big = sum;	row = i;
			}
		}
		if( j != row ){		// need to change pivot
			a = A + row * n;
			b = A + j * n;
			for( k = 0; k < n; k++, a++, b++ ){	// row interchange
				sum = *a;
				*a = *b;
				*b = sum;
			}
			v[ row ] = v[j];
		}
		index[j] = row;				// row index
		a = A + j * n + j;
		if( *a == 0 ) *a = FLT_MIN;	// sigular case
		if( j < n - 1 ){			// devided by pivot
			sum = 1.0 / *a;
			a += n;
			for( i = j + 1; i < n; i++, a += n ) *a *= sum;
		}
	}
	// free memory pool
	if( mem == NULL ) free( memory );
	return 0;
}

void K_SolveLU( double *A, int *ind, double *b, int n )
{
	double sum;
	int i, j, ip;

	for( i = 0; i < n; i++ ){		// solve Ly = b part
		ip = ind[i];
		sum = b[ip];
		b[ip] = b[i];
		if( i ){
			for( j = 0; j < i; j++ ) sum -= A[ i * n + j ] * b[j];
		}
		b[i] = sum;
	}
	for( i = n - 1; i >= 0; i-- ){	// Solve Ux = y part
		sum = b[i];
		for( j = i + 1; j < n; j++ ) sum -= A[ i * n + j ] * b[j];
		b[i] = sum / A[ i * n + i ];
	}
}

void K_LUcomb( double *A, int *ind, int n )
{
	double sum, *a, *b;
	int i, j, k;

	for( j = n - 1; j >= 0; j-- ){
		for( i = n - 1; i > j; i-- ){		// case i >= j 
			sum = 0.0;
			for( k = 0; k <= j; k++ ) sum += A[ i * n + k ] * A[ k * n + j ];
			A[ i * n + j ] = sum;
		}
		for( i = j; i >= 0; i-- ){			// acse i < j
			sum = 0.0;
			for( k = 0; k < i; k++ ) sum += A[ i * n + k ] * A[ k * n + j ];
			A[ i * n + j ] += sum;
		}
	}
	for( i = 0; i < n; i++ ){			// row interchange
		j = ind[i];
		if( j == i ) continue;
		// change!
		a = A + i * n;
		b = A + j * n;
		for( k = 0; k < n; k++, a++, b++ ){
			sum = *a;
			*a = *b;
			*b = sum;
		}
	}
}

double K_InvIter( double *A, double *b0, double shift, int n, int time, double *memory )
{
	double *b1, sum, g, prev;
	int *ind, i, k;
	void *mem;

	// memory pool
	if( memory == NULL ) mem = malloc( n * ( sizeof( double ) + sizeof( int ) ) );
	else mem = memory;
	b1 = ( double * )mem;
	ind = ( int* )( b1 + n );

	// initialize
	for( i = 0; i < n; i++ ) A[ i * n + i ] -= shift;
	K_LUfact( A, ind, n, b1 );

	// keep guess when grown factor too small
	do{
		for( i = 0; i < n; i++ ) b0[i] = b1[i] = ( double )rand() / ( double )RAND_MAX;
		// normalize for both
		sum = 0.0;
		for( i = 0; i < n; i++ ) sum += SQR( b0[i] );
		sum = sqrt( sum );
		for( i = 0; i < n; i++ ) b0[i] = b1[i] = b0[i] / sum;
		// solve ( A- aI )b1 = b0
		K_SolveLU( A, ind, b1, n );
		// normalized b1
		g = 0.0;					// grown factor: | b0 |
		for( i = 0; i < n; i++ ) g += SQR( b1[i] );
		g = sqrt( g );
	} while( g < ERR );	
	
	// normalize
	prev = 0.0;							// | b1 - b0 |
	for( i = 0; i < n; i++ ) {	
		b1[i] *= 1.0/g;						// nomralize
		prev += SQR( b0[i] - b1[i] );	
		b0[i] = b1[i];					// updata b0
	}

	for( k = 0; k < time; k++ ){
		// solve ( A- aI )y = b
		K_SolveLU( A, ind, b1, n );
		// normalized
		sum = 0.0;
		for( i = 0; i < n; i++ ) sum += SQR( b1[i] );
		sum = sqrt( sum );
		g = 0.0;					// | b1 - b0 |
		for( i = 0; i < n; i++ ) {
			b1[i] *= 1.0/sum;			// normalize 
			g += SQR( b0[i] - b1[i] );	// | b1 - b0 |
		}
		// break contition
		if( g < ERR ) break;
		// updata shift value when shift of | b1 - b0 | < error bound
		if( prev - g < ERR ){
			prev = g;
			g = 0.0;
			for( i = 0; i < n; i++ ) g += b0[i] * b1[i];
			g = 1.0 / ( g * sum );	// unnormalized y 
			K_LUcomb( A, ind, n );
			for( i = 0; i < n; i++ ) A[ i * n + i ] += shift;
			shift += g;
			for( i = 0; i < n; i++ ) A[ i * n + i ] -= shift;
			K_LUfact( A, ind, n, b0 );
		}
		else prev = g;
		// updata b0
		for( i = 0; i < n; i++ ) b0[i] = b1[i];
	}

	if( memory == NULL ) free( mem );
	return shift;
}


#endif