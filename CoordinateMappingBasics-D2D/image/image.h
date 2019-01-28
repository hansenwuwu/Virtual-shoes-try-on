#ifndef _BASIC_IMAGE_CLASS_H_
#define _BASIC_IMAGE_CLASS_H_


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <windows.h>

#include "include\jpeglib.h"
#include "include\version.h"

//
// Global variables
//

static char gstring[500];

/**/
//
// function template
//

template < typename T > BYTE CheckColor( T x ){
	if( ( float )x < 0.5f ) return 0x00;
	if( ( float )x > 254.5f ) return 0xff;
	else return ( BYTE )( x + 0.5 );
}

template < typename T > BYTE CheckColor2( T x ){
	if( ( float )x < 0.5f ) return 0x00;
	if( ( float )x > 65535.5f ) return 0xff;
	else return ( BYTE )( x + 0.5 );
}

//
// Basic ImageR Class
//

template < typename T >
class ImageR{
protected:
	T* _bits;
	int _w;
	int _h;
	int _stride;
	int _color;

public:
	ImageR( void );
	~ImageR( void );

	// Basic Functions
	virtual int Initialize( int w, int h, int color );	// allocate image buffer
	virtual int Height( void );
	virtual int Width( void );
	virtual int Stride( void );
	virtual int Color( void );							// 1, 3, 4
	virtual int Empty( void );
	virtual void SetEmpty( void );						// empty without freeing memory block
	virtual void Clear( void );							// free image buffer

	// File Processing
	virtual int ReadFile( const char * );					// check the file format, and read
	virtual int WriteBMP( const char * );
	virtual int ReadBMP( const char * );
	virtual int WriteJPEG( const char *, const int );
	virtual int ReadJPEG( const char * );
	virtual int WritePPM( const char *, const int );
	virtual int ReadPPM( const char * );

	// Editting Functions
	virtual T* GetImage( void );
	virtual int SetImage( T* bits, int w, int h, int color );	// assign the image only
	virtual T GetPixel( int x, int y, int );
	virtual int SetPixel( int x, int y, int, T );
	virtual int SetColor( int );
	virtual int Scale( int w, int h );
};

//
// Definition
//

template < typename T >
ImageR< T >::ImageR( void ){
	_bits = NULL;
	_w = _h = _stride = _color = 0;
}

template < typename T >
ImageR< T >::~ImageR( void ){
	Clear();
}

// 
// Basic Functions
//

template < typename T >
int ImageR< T >::Initialize( int w, int h, int color ){
	// check this
	if( _bits ) Clear();	
	// check image size
	if( w < 0 || h < 0 ) return 0;
	// stride setting
	_color = color;
	_stride = _color * w;
	if( _stride % 4 ) _stride += 4 - _stride % 4;
	_h = h;
	_w = w;
	// memory allocate
	_bits = ( T* )malloc( _stride * _h * sizeof( T ) );
	if( _bits == NULL ) return 0;
	return 1;
}

template < typename T >
int ImageR< T >::Height( void ){ return _h; }

template < typename T >
int ImageR< T >::Width( void ){ return _w; }

template < typename T >
int ImageR< T >::Stride( void ){ return _stride; }

template < typename T >
int ImageR< T >::Color( void ){ return _color; }

template < typename T >
int ImageR< T >::Empty( void ){ return ( _bits ) ? 0 : 1; }

template < typename T >
void ImageR< T >::SetEmpty( void ){
	_bits = NULL;
	_w = _h = _stride = _color = 0;
}

template < typename T >
void ImageR< T >::Clear( void ){ if( _bits ){ free( _bits ); _bits = NULL; } }

//
// File Processing
//

template < typename T >
int ImageR< T >::ReadFile( const char *path ){
	char *ptr = strrchr( ( char* )path, '.' ) + 1;

	if( ptr - 1 == NULL ) return 0;

	// JPEG format
	if( strcmp( ptr, "jpg" ) == 0 ) return ReadJPEG( path );
	if( strcmp( ptr, "jpeg" ) == 0 ) return ReadJPEG( path );
	if( strcmp( ptr, "jpe" ) == 0 ) return ReadJPEG( path );
	if( strcmp( ptr, "jfif" ) == 0 ) return ReadJPEG( path );
	// BMP format
	if( strcmp( ptr, "bmp" ) == 0 ) return ReadBMP( path );
	if( strcmp( ptr, "dib" ) == 0 ) return ReadBMP( path );
	// PPM format
	if( strcmp( ptr, "ppm" ) == 0 ) return ReadPPM( path );

	return 0;		// no feasible format
}

template < typename T >
int ImageR< T >::WriteBMP( const char *path ){
	FILE *outfile;
	T *imgbit;
	BYTE *bits, *bit;
	int i, j;
	DWORD buf;

	if( Empty() ){			// checking
		//OutputDebugString( "can't write bmp file: image is empty\n" );
		return 0;
	}

	// open output file
#if _NOT_VC6_
	fopen_s( &outfile, ( char * )path, "wb" );
#else
	outfile = fopen( ( char * )path, "wb" );
#endif
	if( outfile == NULL ){
#if _NOT_VC6_
		sprintf_s( gstring, 500, "can't open bmp file: %s\n", path );
#else
		sprintf( gstring, "can't open bmp file: %s\n", path );
#endif
		//OutputDebugString( gstring );
		return 0;
	}

	// write file header
	fwrite( &( buf = 0x4D42 ), 2, 1, outfile );
	fwrite( &( buf = _stride * _h + 54 ), 4, 1, outfile );
	fwrite( &( buf = 0x00000000 ), 4, 1, outfile );
	fwrite( &( buf = ( _color == 1 ) ? 1078 : 54 ), 4, 1, outfile );
	fwrite( &( buf = 40 ), 4, 1, outfile );
	fwrite( &_w, 4, 1, outfile );
	fwrite( &_h, 4, 1, outfile );
	fwrite( &( buf = 0x0001 ), 2, 1, outfile );
	fwrite( &( buf = _color * 8 ), 2, 1, outfile );
	fwrite( &( buf = 0x00 ), 4, 1, outfile );
	fwrite( &( buf = _w * _h * _color ), 4, 1, outfile );
	fwrite( &( buf = 0x00000000 ), 4, 2, outfile );
//	fwrite( &buf, 4, 1, outfile );
	fwrite( &( buf = ( _color == 1 ) ? 256 : 0x00000000 ), 4, 1, outfile );
	fwrite( &( buf = 0x00000000 ), 4, 1, outfile );

	// write gray-scale color platte for 256 color bitmap
	if( _color == 1 ){
		for( i = 0; i < 256; i++ ){
			j = i * 0x10000 + i * 0x100 + i;
			fwrite( &j, 1, 4, outfile );
		}
	}
	// generate image buffer
	bits = ( BYTE* )malloc( _stride * _h );

	// color transfer
	switch( _color ){
	case 1:
		imgbit = _bits;
		bit = bits;
		for( i = 0, j = _stride * _h; i < j; i++, imgbit++, bit++ )
			bit[0] = CheckColor( imgbit[0] );
		break;
	case 3:
		for( i = 0; i < _h; i++ ){
			imgbit = _bits + i * _stride;
			bit = bits + i * _stride;
			for( j = 0; j < _stride; j += 3, imgbit += 3, bit += 3 ){
				bit[0] = bit[1] = bit[2] = CheckColor( imgbit[2] );
				bit[1] = CheckColor( imgbit[1] );
				bit[2] = CheckColor( imgbit[0] );
			}
		}
		break;
	case 4:
		for( i = 0; i < _h; i++ ){
			imgbit = _bits + i * _stride;
			bit = bits + i * _stride;
			for( j = 0; j < _stride; j += 4, imgbit += 4, bit += 4 ){
				bit[0] = CheckColor( imgbit[2] );
				bit[1] = CheckColor( imgbit[1] );
				bit[2] = CheckColor( imgbit[0] );
				bit[3] = CheckColor( imgbit[3] );
			}
		}
		break;
	}
	// write bitmap
	fwrite( bits, 1, _stride * _h, outfile );

	// free memory
	free( bits );
	fclose( outfile );
	return 1;
}

template < typename T >
int ImageR< T >::ReadBMP( const char *path ){
	FILE *infile;
	T *imgbit;
	BYTE *bits, *bit, *platte, *color;
	BITMAPINFOHEADER bmpinfoheader;
	BITMAPFILEHEADER bmpfileheader;
	int i, j, stride;

	if( Empty() == 0 ) Clear();

	// open output file
#if _NOT_VC6_
	fopen_s( &infile, path, "rb" );
#else
	infile = fopen( path, "rb" );
#endif
	if( infile == NULL ){
#if _NOT_VC6_
		sprintf_s( gstring, 500, "can't open bmp file: %s\n", path );
#else
		sprintf( gstring, "can't open bmp file: %s\n", path );
#endif
		//OutputDebugString( gstring );
		return 0;
	}

	// read bitmap file header
	fread( &bmpfileheader, sizeof( BITMAPFILEHEADER ), 1, infile );
	if( bmpfileheader.bfType != 0x4D42 ) {// check file magic number
#if _NOT_VC6_
		sprintf_s( gstring, 500, "it's not a bmp file: %s\n", path );
#else
		sprintf( gstring, "it's not a bmp file: %s\n", path );
#endif
		//OutputDebugString( gstring );
		fclose( infile );
		return 0;
	}

	// read bitmap infomation header
	fread( &bmpinfoheader, sizeof( BITMAPINFOHEADER ), 1, infile );
	_w = bmpinfoheader.biWidth;
	_h = bmpinfoheader.biHeight;
	_color = ( bmpinfoheader.biBitCount == 24 ) ? 3 : 4;			// save as RGB or RGBA

	// color platte setting
	if( bmpinfoheader.biBitCount < 24 && bmpinfoheader.biClrUsed == 0 ){
		// set ClrUsed
		bmpinfoheader.biClrUsed = 0x00000001 << bmpinfoheader.biBitCount;
	}
	if( bmpinfoheader.biClrUsed ){
		platte = ( BYTE* )malloc( bmpinfoheader.biClrUsed * 4 );
		fread( platte, 1, bmpinfoheader.biClrUsed * 4, infile );
	}

	// memory setting
	_stride = _w * _color;
	if( _stride % 4 ) _stride += 4 - _stride % 4;
	_bits = ( T* )malloc( _h * _stride * sizeof( T ) );				// image
	bits = ( BYTE* )malloc( bmpinfoheader.biSizeImage );			// image file array
	fread( bits, 1, bmpinfoheader.biSizeImage, infile );			// read bitmap

#define SetColor \
	imgbit[0] = ( T ) color[2]; imgbit[1] = ( T ) color[1]; imgbit[2] = ( T ) color[0];\
	imgbit[3] = ( T ) color[3]; imgbit += 4;

	// color transfer
	stride = bmpinfoheader.biSizeImage / _h;
	switch( bmpinfoheader.biBitCount ){
	case 1:
		for( i = 0; i < _h; i++ ){
			bit = bits + i * stride;
			imgbit = _bits + i * _stride;
			for( j = 0; j < stride; j++, bit++ ){
				color = platte + 4 * ( ( bit[0] & 0x80 ) ? 1 : 0 ); SetColor;
				color = platte + 4 * ( ( bit[0] & 0x40 ) ? 1 : 0 ); SetColor;
				color = platte + 4 * ( ( bit[0] & 0x20 ) ? 1 : 0 ); SetColor;
				color = platte + 4 * ( ( bit[0] & 0x10 ) ? 1 : 0 ); SetColor;
				color = platte + 4 * ( ( bit[0] & 0x08 ) ? 1 : 0 ); SetColor;
				color = platte + 4 * ( ( bit[0] & 0x04 ) ? 1 : 0 ); SetColor;
				color = platte + 4 * ( ( bit[0] & 0x02 ) ? 1 : 0 ); SetColor;
				color = platte + 4 * ( ( bit[0] & 0x01 ) ? 1 : 0 ); SetColor;
			}
		}
		break;
	case 4:
		for( i = 0; i < _h; i++ ){
			bit = bits + i * stride;
			imgbit = _bits + i * _stride;
			for( j = 0; j < stride; j++, bit++ ){
				color = platte + 4 * ( ( bit[0] & 0xf0 ) >> 4 ); SetColor;
				color = platte + 4 * ( ( bit[0] & 0x0f ) ); SetColor;
			}
		}
		break;
	case 8:
		for( i = 0; i < _h; i++ ){
			bit = bits + i * stride;
			imgbit = _bits + i * _stride;
			for( j = 0; j < stride; j++, bit++ ){
				color = platte + 4 * ( bit[0] ); SetColor;
			}
		}
		break;
	case 16:
		for( i = 0; i < _h; i++ ){
			bit = bits + i * stride;
			imgbit = _bits + i * _stride;
			for( j = 0; j < stride; j++, bit += 2 ){
				color = platte + 4 * ( bit[0] * 0xff + bit[1] ); SetColor;
			}
		}
		break;
	case 24:
		for( i = 0; i < _h; i++ ){
			imgbit = _bits + i * _stride;
			bit = bits + i * _stride;
			for( j = 0; j < _stride; j += 3, imgbit += 3, bit += 3 ){
				imgbit[0] = ( T )bit[2];
				imgbit[1] = ( T )bit[1];
				imgbit[2] = ( T )bit[0];
			}
		}
		break;
	case 32:
		for( i = 0; i < _h; i++ ){
			imgbit = _bits + i * _stride;
			bit = bits + i * _stride;
			for( j = 0; j < _stride; j += 4, imgbit += 4, bit += 4 ){
				imgbit[0] = ( T )bit[2];
				imgbit[1] = ( T )bit[1];
				imgbit[2] = ( T )bit[0];
				imgbit[3] = ( T )bit[3];
			}
		}
		break;
	}
#undef SetColor

	// free memory
	free( bits );
	if( bmpinfoheader.biClrUsed ) free( platte );
	fclose( infile );

	return 1;
}

template < typename T >
int ImageR< T >::ReadJPEG( const char *path ){	
	struct jpeg_decompress_struct jpg;
	struct jpeg_error_mgr err;
	FILE *infile;
	BYTE *bits, *bit;
	T *imgbit;
	int i, j;
	JSAMPROW row[1];
	
	if( Empty() == 0 ) Clear();

	// initialize
	jpg.err = jpeg_std_error( &err );
	jpeg_create_decompress( &jpg );

	// file process
#if _NOT_VC6_
	fopen_s( &infile, path, "rb" );
#else
	infile = fopen( path, "rb" );
#endif

	if( infile == NULL ){
#if _NOT_VC6_
		sprintf_s( gstring, 500, "can't open jpeg file: %s\n", path );
#else
		sprintf( gstring, "can't open jpeg file: %s\n", path );
#endif
		//OutputDebugString( gstring );
		return 0;
	}

	// hpeg header
	jpeg_stdio_src( &jpg, infile );
	jpeg_read_header( &jpg, true );
	jpeg_start_decompress( &jpg );

	// create bitmap buffer
	_w = jpg.image_width;
	_h = jpg.image_height;
	_color = jpg.output_components;
	_stride = _w * _color;
	if( _stride % 4 ) _stride += 4 - _stride % 4;
	bits = ( BYTE* )malloc( _stride * _h );
	_bits = ( T* )malloc( _stride * _h * sizeof( T ) );

	// fail to mamory allocate
	if( _bits == NULL || bits == NULL ){
		//OutputDebugString( "can't allocate memory\n" );
		free(bits);
		fclose( infile );
		return 0;
	}

	// decompress
	while( jpg.output_scanline < jpg.output_height ){
		row[0] = bits + ( jpg.output_height - jpg.output_scanline - 1 ) * _stride;
		jpeg_read_scanlines( &jpg, row, 1 );
	}
	jpeg_finish_decompress( &jpg );

	// data transfer
	for( i = 0; i < _h; i++ ){
		bit = bits + i * _stride;
		imgbit = _bits + i * _stride;
		for( j = 0; j < _stride; j++, bit++, imgbit++ ) 
			imgbit[0] = ( T )bit[0];
	}

	// free
	jpeg_destroy_decompress( &jpg );
	free( bits );
	fclose( infile );

	return 1;
}

template < typename T >
int ImageR< T >::WriteJPEG( const char *path, const int quality ){	
	struct jpeg_compress_struct jpg;
	struct jpeg_error_mgr err;
	FILE *outfile;
	BYTE *bits, *bit;
	T *imgbit;
	JSAMPROW row[1];
	int i, j;

	if( Empty() ){			// checking
		//OutputDebugString( "can't write jpeg file: image is empty\n" );
		return 0;
	}

	// initialize
	jpg.err = jpeg_std_error( &err );
	jpeg_create_compress( &jpg );

	// file process
#if _NOT_VC6_
	fopen_s( &outfile, path, "wb" );
#else
	outfile = fopen( path, "wb" );
#endif
	if( outfile == NULL ){
#if _NOT_VC6_
		sprintf_s( gstring, 500, "can't write jpeg file: %s\n", path );
#else
		sprintf( gstring, "can't write jpeg file: %s\n", path );
#endif
		//OutputDebugString( gstring );
		return 0;
	}
	jpeg_stdio_dest( &jpg, outfile );

	// parameter setting
	jpg.image_width = _w;
	jpg.image_height = _h;
	jpg.input_components = ( _color == 1 ) ? 1 : 3;
	jpg.in_color_space = ( _color == 1 ) ? JCS_GRAYSCALE : JCS_RGB;
	jpeg_set_defaults( &jpg );
	jpeg_set_quality( &jpg, quality, true );
/**/
	// type transfer
	bits = ( BYTE* )malloc( _h * _stride );
	switch( _color ){
	case 1:
	case 3:
		bit = bits;
		imgbit = _bits;
		for( i = 0, j = _stride * _h; i < j; i++, bit++, imgbit++ ) 
			bit[0] = CheckColor( imgbit[0] );
		break;
	case 4:
		// eliminate the alpha chanel
		for( i = 0; i < _h; i++ ){
			bit = bits + i * _stride;
			imgbit = _bits + i * _stride;
			for( j = 0; j < _w; j++, bit += 3, imgbit += 4 ){
				bit[0] = CheckColor( imgbit[0] );
				bit[1] = CheckColor( imgbit[1] );
				bit[2] = CheckColor( imgbit[2] );
			}
		}
		break;
	}
/**/
	// compress
	jpeg_start_compress( &jpg, true );
	while( jpg.next_scanline < jpg.image_height ){
		row[0] = bits + ( _h - jpg.next_scanline ) * _stride;
		jpeg_write_scanlines( &jpg, row, 1 );
	}
	jpeg_finish_compress( &jpg );

	// free memory
	jpeg_destroy_compress( &jpg );
	free( bits );
	fclose( outfile );

	return 1;
}

template < typename T >
int ImageR< T >::ReadPPM( const char *path ){
	FILE *infile;
	T *bit;
	void *buf;
	unsigned char *p1;
	unsigned short int *p2;
	int i, j, bytes, binary, w3;

	// file open
	if( Empty() == 0 ) Clear();			// clear itself
#if _NOT_VC6_
	fopen_s( &infile, path, "r" );
#else
	infile = fopen( path, "r" );
#endif

	if( infile == NULL ){
#if _NOT_VC6_
		sprintf_s( gstring, 500, "can't read PPM file: %s\n", path );
#else
		sprintf( gstring, "can't read PPM file: %s\n", path );
#endif
		//OutputDebugString( gstring );
		return 0;
	}

	// check magic number
#if _NOT_VC6_
	fscanf_s( infile, "%s\n", gstring );
#else
	fscanf( infile, "%s\n", gstring );
#endif
	if( gstring[0] != 'P' ){
		//OutputDebugString( "It's not a PPM file\n" );
		fclose( infile );
		return 0;
	}
	switch( gstring[1] ){
	case '3': binary = 0; break;
	case '6': binary = 1; break;
	default: 
		//OutputDebugString( "It's not a PPM file\n" ); 
		fclose( infile ); return 0;
	}

	// read header
#if _NOT_VC6_
	fscanf_s( infile, "%d %d %d\n", &_w, &_h, &bytes );
#else
	fscanf( infile, "%d %d %d\n", &_w, &_h, &bytes );
#endif
	// setting image
	_color = 3;							// must RGB
	_stride = w3 = _w * 3;
	if( _stride % 4 ) _stride += 4 - _stride % 4;
	// allocate image memory buffer
	_bits = ( T* )malloc( _stride * _h * sizeof( T ) );
	if( _bits == NULL ){
		//OutputDebugString( "fail to allocate memory\n" );
		fclose( infile );
		return 0;
	}

	// read image array
#define Loading( x ) \
	for( i = _h - 1; i >= 0; i-- ){ bit = _bits + i * _stride; 	for( j = 0; j < w3; j++, bit++, p##x++ ) bit[0] = ( T )p##x[0]; }

	bytes = ( bytes >> 8 ) ? 2 : 1;			// if max > 255 => 2 bytes
	if( binary ){							// binary mode
		// re-open file in binary mode
		i = ftell( infile );				// get current position
#if _NOT_VC6_
		freopen_s( &infile, path, "rb", infile );
#else
		infile = freopen( path, "rb", infile );
#endif
		fseek( infile, i - 3, SEEK_SET );
		buf = malloc( w3 * _h * bytes );
		if( buf == NULL ){
			//OutputDebugString( "fail to allocate memory\n" );
			fclose( infile );
			return 0;
		}
		fread( buf, bytes, w3 * _h, infile );
		if( bytes == 1 ) { p1 = ( unsigned char * )buf; Loading( 1 ); }
		else if( bytes == 2 ) { p2 = ( unsigned short int * )buf; Loading( 2 ); }
		free( buf );
	}
	else{
		for( i = _h - 1; i >= 0; i-- ){
			bit = _bits + i * _stride;
			for( j = 0; j < w3; j++, bit++ ) 
#if _NOT_VC6_
				bit[0] = fscanf_s( infile, "%d", bit );
#else
				bit[0] = fscanf( infile, "%d", bit ); 
#endif
		}
	}
#undef Loading
	fclose( infile );
	return 1;
}

template < typename T >
int ImageR< T >::WritePPM( const char *path, const int max ){
	FILE *outfile;
	T* bit;
	void *buf;
	unsigned char *p1;
	unsigned short int *p2;
	const int bytes = ( max >> 8 ) ? 2 : 1;	// if max > 255 => 2 bytes
	const int w3 = _w * 3;
	int i, j;

	if( Empty() ){			// checking
		//OutputDebugString( "can't write jpeg file: image is empty\n" );
		return 0;
	}

	// check file open

#if _NOT_VC6_
	fopen_s( &outfile, path, "w" );
#else
	outfile = fopen( path, "w" );
#endif

	if( outfile == NULL ){
#if _NOT_VC6_
		sprintf_s( gstring, 500, "can't write PPM file: %s\n", path );
#else
		sprintf( gstring, "can't write PPM file: %s\n", path );
#endif
		//OutputDebugString( gstring );
		return 0;
	}

	// generate output buffer
	buf = malloc( _w * _h * 3 * bytes );	// must RGB only
	if( buf == NULL ){
		//OutputDebugString( "fail to allocate memory\n" );
		fclose( outfile );
		return 0;
	}

	// setting buffer
#define CheckColor1	CheckColor
#define SetGray( x ) \
	for( i = _h - 1; i >= 0; i-- ){ bit = _bits + i * _stride; 	for( j = 0; j < _w; j++, bit++, p##x += 3 )	p##x[0] = p##x[1] = p##x[2] = CheckColor##x( bit[0] ); }
#define SetRGB( x ) \
	for( i = _h - 1; i >= 0; i-- ){ bit = _bits + i * _stride; 	for( j = 0; j < w3; j++, bit++, p##x++ ) p##x[0] = CheckColor##x( bit[0] ); }
#define SetRGBA( x ) \
	for( i = _h - 1; i >= 0; i-- ){ bit = _bits + i * _stride; 	for( j = 0; j < _w; j++, bit += 4, p##x += 3 ) \
	p##x[0] = CheckColor##x( bit[0] ); p##x[1] = CheckColor##x( bit[1] ); p##x[2] = CheckColor##x( bit[2] ); }

	switch( _color ){
	case 1:
		switch( bytes ){
		case 1:								// gray-scale, 1 bytes
			p1 = ( unsigned char * )buf;
			SetGray( 1 ); break;
		case 2:								// gray-scale, 2 bytes
			p2 = ( unsigned short int * )buf;
			SetGray( 2 ); break;
		}
		break;
	case 3:
		switch( bytes ){
		case 1:								// RGB-scale, 1 bytes
			p1 = ( unsigned char * )buf;
			SetRGB( 1 ); break;
		case 2:				// gray-scale, 2 bytes
			p2 = ( unsigned short int * )buf;
			SetRGB( 2 ); break;
		}
		break;
	case 4:
		switch( bytes ){
		case 1:								// RGB-scale, 1 bytes
			p1 = ( unsigned char * )buf;
			SetRGBA( 1 ); break;
		case 2:				// gray-scale, 2 bytes
			p2 = ( unsigned short int * )buf;
			SetRGBA( 2 ); break;
		}
		break;
	}
#undef CheckColor1
#undef SetGray
#undef SetRGB
#undef SetRGBA

	// write magic number P6, width, height, maximum number
	fprintf( outfile, "P6\n%d %d\n%d\n", _w, _h, max );

	// write binary image array
#if _NOT_VC6_
	freopen_s( &outfile, path, "ab", outfile );
#else
	outfile = freopen( path, "ab", outfile );
#endif
	fwrite( buf, bytes, w3 * _h, outfile );
	// free memory
	free( buf );
	fclose( outfile );
	return 1;
}

//
// Editing Functions
//

template < typename T >
T* ImageR< T >::GetImage( void ){ return _bits; }

template < typename T >
int ImageR< T >::SetImage( T* bits, int w, int h, int color ){
	if( _bits ) free( _bits );
	_bits = bits;
	_w = w;
	_h = h;
	// stride setting
	_color = color;
	_stride = _color * _w;
	if( _stride % 4 ) _stride += 4 - _stride % 4;
	return 1;
}

template < typename T >
T ImageR< T >::GetPixel( int x, int y, int c ){
	T *bit = _bits + y * _stride + x * _color;
	return bit[c];
}

template < typename T >
int ImageR< T >::SetPixel( int x, int y, int c, T v ){
	if( x < 0 || x >= _w || y < 0 || y >= _h ) return 0;
	T *bit = _bits + y * _stride + x * _color;
	bit[c] = v;
	return 1;
}

template < typename T >
int ImageR< T >::SetColor( int c ){
	T *bit1, *bit2, *bits;
	int stride, i, j;

	if( c == _color ) return 1;		// do need to change

	switch( c ){
	case 1:		// gray-scale image
		switch( _color ){
		case 3:
		case 4:
			stride = _w;
			if( stride % 4 ) stride += 4 - stride % 4;
			bits = ( T* )malloc( _h * stride * sizeof( T ) );
			if( bits == NULL ){
				//OutputDebugString( "fail to allocate memory\n" );
				return 0;
			}
			for( i = 0; i < _h; i++ ){
				bit1 = _bits + i * _stride;
				bit2 = bits + i * stride;
				for( j = 0; j < _w; j++, bit1 += _color, bit2++ )
					bit2[0] = ( T )( ( ( float )bit1[0] + ( float )bit1[1] + ( float )bit1[2] ) / 3.0f );
			}
			if( _bits == NULL ) free( _bits ); 
			_bits = bits; _stride = stride; _color = 1;
			break;
		default: return 0;
		}
		return 1;
	case 3:		// RGB
		stride = _w * 3;
		if( stride % 4 ) stride += 4 - stride % 4;
		bits = ( T* )malloc( _h * stride * sizeof( T ) );
		if( bits == NULL ){
			//OutputDebugString( "fail to allocate memory\n" );
			return 0;
		}
		switch( _color ){
		case 1:					// from gray-style
			for( i = 0; i < _h; i++ ){
				bit1 = _bits + i * _stride;
				bit2 = bits + i * stride;
				for( j = 0; j < _w; j++, bit1++, bit2 += 3 ) 
					bit2[0] = bit2[1] = bit2[2] = bit1[0];
			}
			break;
		case 4:					// from RGBA
			for( i = 0; i < _h; i++ ){
				bit1 = _bits + i * _stride;
				bit2 = bits + i * stride;
				for( j = 0; j < _w; j++, bit1 += 4, bit2 += 3 )
					memcpy( bit2, bit1, 3 * sizeof( T ) );	// copy RGB only	
			}
			break;
		default: return 0;
		}
		if( _bits == NULL ) free( _bits ); 
		_bits = bits; _stride = stride; _color = 3;
		return 1;
	case 4:		// RGBA
		stride = _w * 4;
		bits = ( T* )calloc( _h * stride, sizeof( T ) );	// initialize to 0
		if( bits == NULL ){
			//OutputDebugString( "fail to allocate memory\n" );
			return 0;
		}
		switch( _color ){
		case 1:					// from gray-style
			for( i = 0; i < _h; i++ ){
				bit1 = _bits + i * _stride;
				bit2 = bits + i * stride;
				for( j = 0; j < _w; j++, bit1++, bit2 += 4 ) 
					bit2[0] = bit2[1] = bit2[2] = bit1[0];
			}
			break;
		case 3:					// from RGB
			for( i = 0; i < _h; i++ ){
				bit1 = _bits + i * _stride;
				bit2 = bits + i * stride;
				for( j = 0; j < _w; j++, bit1 += 3, bit2 += 4 )
					memcpy( bit2, bit1, 3 * sizeof( T ) );	// copy RGB only	
			}
			break;
		default: return 0;
		}
		if( _bits == NULL ) free( _bits ); 
		_bits = bits; _stride = stride; _color = 4;
		return 1;
	default: return 0;
	}
	return 0;
}

template < typename T >
int ImageR< T >::Scale( int w2, int h2 ){
	const int w1 = _w, h1 = _h,
		stride1 = _stride,
		stride2 = ( ( w2 * _color ) % 4 ) ? ( w2 * _color / 4 * 4 + 4 ) : ( w2 * _color );
	int y1, y2, x1, x2, i;
	float v1, v2, u1, u2, temp[4];
	T *img1, *img2, *bit1, *bit2;

	// Set scaled image
	img1 = _bits;
	img2 = ( T* )malloc( stride2 * h2 * sizeof( T ) );
	if( img2 == NULL ) return 0;

	// inverse warping
	for( y2 = 0; y2 < h2; y2++ ){
		bit2 = img2 + y2 * stride2;
		v2 = ( float )y2 / ( float )h2 * ( float )h1;
		y1 = ( int )v2;					// the old pixel
		v2 = v2 - ( float )y1;			// the ratio of interpolation
		v1 = 1.0f - v2;
		for( x2 = 0; x2 < w2; x2++, bit2 += _color ){
			u2 = ( float )x2 / ( float )w2 * ( float )w1;
			x1 = ( int )u2;				// the old pixel
			u2 = u2 - ( float )x1;		// the ratio of interpolation
			u1 = 1.0f - u2;
			bit1 = img1 + y1 * stride1 + x1 * _color;
			// bilinear
			for( i = 0; i < _color; i++ ) 
				temp[i] = u1 * v1 * ( float )bit1[i];
			if( y1 + 1 < h1 ){
				for( i = 0; i < _color; i++ ) 
					temp[i] += u1 * v2 * ( float )bit1[ stride1 + i ];
			}
			if( x1 + 1 < w1 ){
				for( i = 0; i < _color; i++ ) 
					temp[i] += u2 * v1 * ( float )bit1[ _color + i ];
			}
			if( x1 + 1 < w1 && y1 + 1 < h1 ){
				for( i = 0; i < _color; i++ ) 
					temp[i] += u2 * v2 * ( float )bit1[ _color + stride1 + i ] ;
			}
			for( i = 0; i < _color; i++ ) bit2[i] = ( T )temp[i];
		}
	}
	// setting
	free( _bits );	_bits = img2;
	_h = h2;	_w = w2;	_stride = stride2;
	return 1;
}



#endif