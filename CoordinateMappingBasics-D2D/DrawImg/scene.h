
#ifndef __SCENE__
#define __SCENE__


#ifdef __cplusplus
extern "C" {
#endif

/*	struct define
_______________________________________________________________*/

	typedef struct _point_xy_int			point2i;
	typedef struct _point_xy_float			point2f;
	typedef struct _point_xyz_int			point3i;
	typedef struct _point_xyz_float			point3f, vector3f;
	typedef struct _color_rgb_float			color3f, rgb;
	typedef struct _color_rgb_uchar			color3uc, rgbuc;
	typedef struct _color_rgba_float		color4f, rgba;
	typedef struct _color_rgba_uchar		color4uc, rgbauc;
	
	typedef struct _rectangle_float			rectanglef;
	typedef struct _rectangle_detail_int	rectangleDi;
	typedef struct _rectangle_detail_float	rectangleDf;
	typedef struct _ellipse_xy_int			ellipsei;
	typedef struct _ellipse_xy_float		ellipsef;
	typedef struct _triangle_xyz_float		trianglef;
	typedef struct _size_wh_int				size2i;
	typedef struct _size_wh_float			size2f;

/*	point :	int x, y
_______________________________________________________________*/

	struct _point_xy_int {
		int x, y;

		_point_xy_int( void ) { x = y = 0; }
		_point_xy_int( int _n ) { x = y = _n; }
		_point_xy_int( int _x, int _y ) { x = _x; y = _y; }

		
	};



/*	point :	float x, y
_______________________________________________________________*/

	struct _point_xy_float {
		float x, y;

		_point_xy_float( void ) { x = y = 0.f; }
		_point_xy_float( float _n ) { x = y = _n; }
		_point_xy_float( float _x, float _y ) { x = _x; y = _y; }

		
	};



/*	point :	int x, y, z
_______________________________________________________________*/

	struct _point_xyz_int {
		int x, y, z;

		_point_xyz_int( void ) { x = y = z = 0; }
		_point_xyz_int( int _n ) { x = y = z = _n; }
		_point_xyz_int( int _x, int _y, int _z ) { x = _x; y = _y; z = _z; }

		
	};



/*	point :	float x, y, z
_______________________________________________________________*/

	struct _point_xyz_float {
		float x, y, z;

		_point_xyz_float( void ) { x = y = z = 0.f; }
		_point_xyz_float( float _n ) { x = y = z = _n; }
		_point_xyz_float( float _x, float _y, float _z ) { x = _x; y = _y; z = _z; }

		
	};



/*	color :	float r, g, b
_______________________________________________________________*/

	struct _color_rgb_float {
		float r, g, b;

		_color_rgb_float( void ) { r = g = b = 0.f; }
		_color_rgb_float( float _n ) { r = g = b = _n; }
		_color_rgb_float( float _r, float _g, float _b ) { r = _r; g = _g; b = _b; }

		
	};



/*	color :	unsigned char r, g, b
_______________________________________________________________*/

	struct _color_rgb_uchar {
		unsigned char r, g, b;

		_color_rgb_uchar( const _color_rgb_float &_s ) { r = (int)( _s.r*255 ); g = (int)( _s.g*255 ); b = (int)( _s.b*255 ); }
		_color_rgb_uchar( void ) { r = g = b = 0; }
		_color_rgb_uchar( unsigned char _n ) { r = g = b = _n; }
		_color_rgb_uchar( unsigned char _r, unsigned char _g, unsigned char _b ) { r = _r; g = _g; b = _b; }

		
	};



/*	color :	float r, g, b, a
_______________________________________________________________*/

	struct _color_rgba_float {
		float r, g, b, a;

		_color_rgba_float( void ) { r = g = b = 0.f; a = 1.f; }
		_color_rgba_float( float _n ) { r = g = b = _n; a = 1.f; }
		_color_rgba_float( const color3f _c, float _a ) { r = _c.r; g = _c.g; b = _c.b; a = _a; }
		_color_rgba_float( float _r, float _g, float _b, float _a ) { r = _r; g = _g; b = _b; a = _a; }

		
	};



/*	color :	unsigned char r, g, b, a
_______________________________________________________________*/

	struct _color_rgba_uchar {
		unsigned char r, g, b, a;

		_color_rgba_uchar( void ) { r = g = b = 0; a = 255; }
		_color_rgba_uchar( unsigned char _n ) { r = g = b = _n; a = 255; }
		_color_rgba_uchar( const color3uc _c, unsigned char _a ) { r = _c.r; g = _c.g; b = _c.b; a = _a; }
		_color_rgba_uchar( const color4f _c ) { r = (unsigned char)( 255*_c.r ); g = (unsigned char)( 255*_c.g ); b = (unsigned char)( 255*_c.b ); a = (unsigned char)( 255*_c.a ); }
		_color_rgba_uchar( unsigned char _r, unsigned char _g, unsigned char _b, unsigned char _a ) { r = _r; g = _g; b = _b; a = _a; }

		
	};



/*	rectangle :	int location( x, y )
				int width, height
_______________________________________________________________*/

	struct _rectangle_int {
		point2i location;
		int width, height, area;

		_rectangle_int( void ) { width = height = area = 0; }
		_rectangle_int( int _x, int _y, int _w, int _h ) { location = point2i( _x, _y ); width = _w; height = _h; area = _w*_h; }
		_rectangle_int( point2i &_a, point2i &_b )
		{
			if( _a.x < _b.x ) { location.x = _a.x; width = _b.x - _a.x; } else { location.x = _b.x; width = _a.x - _b.x; }
			if( _a.y < _b.y ) { location.y = _a.y; height = _b.y - _a.y; } else { location.y = _b.y; height = _a.y - _b.y; }
			area = width*height;
		}
		
		
	};



/*	rectangle :	float location( x, y )
				float width, height
_______________________________________________________________*/

	struct _rectangle_float {
		point2f location;
		float width, height, area;

		_rectangle_float( void ) { width = height = area = 0.f; }
		_rectangle_float( float _x, float _y, float _w,float _h ) { location = point2f( _x, _y ); width = _w; height = _h; area = _w*_h; }
		_rectangle_float( point2f &_a, point2f &_b )
		{
			if( _a.x < _b.x ) { location.x = _a.x; width = _b.x - _a.x; } else { location.x = _b.x; width = _a.x - _b.x; }
			if( _a.y < _b.y ) { location.y = _a.y; height = _b.y - _a.y; } else { location.y = _b.y; height = _a.y - _b.y; }
			area = width*height;
		}
		
		
	};



/*	rectangle :	int location( x, y )
				int top_left( x, y ), top_right( x, y ), bottom_left( x, y ), bottom_right( x, y )
				int width, height
_______________________________________________________________*/



	struct _rectangle_detail_float {
		point2f location;
		float width, height, area;
		point2f top_left, top_right, bottom_left, bottom_right;
		
		void const SideCalculate( void )
		{
			area = width*height;
			top_left = location;
			top_right = point2f( location.x + width, location.y );
			bottom_left = point2f( location.x, location.y - height );
			bottom_right = point2f( location.x + width, location.y - height );
		}

		_rectangle_detail_float( void ) { width = height = area = 0; }
		_rectangle_detail_float( float _x, float _y, float _w, float _h )
		{
			location = point2f( _x, _y );
			width = _w;
			height = _h;
			SideCalculate();
		}		
		_rectangle_detail_float( point2i &_a, point2i &_b )
		{
			if( _a.x < _b.x ) { location.x = (float)_a.x; width = (float)( _b.x - _a.x ); } else { location.x = (float)_b.x; width = (float)( _a.x - _b.x ); }
			if( _a.y < _b.y ) { location.y = (float)_a.y; height = (float)( _b.y - _a.y ); } else { location.y = (float)_b.y; height = (float)( _a.y - _b.y ); }
			SideCalculate();
		}
		_rectangle_detail_float( point2f &_a, point2f &_b )
		{
			if( _a.x < _b.x ) { location.x = _a.x; width = _b.x - _a.x; } else { location.x = _b.x; width = _a.x - _b.x; }
			if( _a.y > _b.y ) { location.y = _a.y; height = _b.y - _a.y; } else { location.y = _b.y; height = _a.y - _b.y; }
			SideCalculate();
		}
		_rectangle_detail_float( _rectangle_float &_s )
		{
			location = _s.location;
			width = _s.width;
			height = _s.height;
			SideCalculate();
		}

		
	};



/*	ellipse :	int center( x, y )
				int alfa : length of long-axis
				int beta : length of short-axis
				int theta : angle between long-axis and x-axis
_______________________________________________________________*/

	struct _ellipse_xy_int {
		point2i center;
		int alfa, beta, theta;

		_ellipse_xy_int( void ) { alfa = beta = theta = 0; }
		_ellipse_xy_int( point2i &_s, int _a, int _b, int _t ) { center = _s; alfa = _a; beta = _b; theta = _t; }
		
		
	};



/*	ellipse :	float center( x, y )
				float alfa : length of long-axis
				float beta : length of short-axis
				float theta : angle between long-axis and x-axis
_______________________________________________________________*/

	struct _ellipse_xy_float {
		point2f center;
		float alfa, beta, theta;

		_ellipse_xy_float( void ) { alfa = beta = theta = 0; }
		_ellipse_xy_float( point2f &_s, float _a, float _b, float _t ) { center = _s; alfa = _a; beta = _b; theta = _t; }
		
		
	};



/*	triangle :	float p1( x, y, z ), p2( x, y, z ), p3( x, y, z )
_______________________________________________________________*/

	struct _triangle_xyz_float {
		point3f p1, p2, p3;

		_triangle_xyz_float( void ) { p1 = p2 = p3 = 0.f; }
		_triangle_xyz_float( point3f &_s1, point3f &_s2, point3f &_s3 ) { p1 = _s1; p2 = _s2; p3 = _s3; }

		
	};


	
/*	size :	int width, height
_______________________________________________________________*/
	
	struct _size_wh_int {
		int width, height;

		_size_wh_int( void ) { width = height = 0; }
		_size_wh_int( int _w, int _h ) { width = _w; height = _h; }
		
		
	};


	
/*	size :	float width, height
_______________________________________________________________*/
	
	struct _size_wh_float {
		float width, height;

		_size_wh_float( void ) { width = height = 0; }
		_size_wh_float( float _w, float _h ) { width = _w; height = _h; }
		
		
	};



#ifdef __cplusplus
}
#endif

#endif
