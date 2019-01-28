#ifndef _IMAGE_H_
#define _IMAGE_H_

#include <windows.h>

//
// Macro
//

#define IMAGE_GRAY	8
#define IMAGE_RGB	24	
#define IMAGE_RGBA	32

// Color Spaces

enum BGR{ B, G, R};
enum LUV{ L, U, V};
enum XYZ{ X, Y, Z};

class IMAGE{
protected:
	BITMAPINFO _inf;
	BITMAPFILEHEADER _fileheader;
	BITMAP _bmp;
	HBITMAP _hbmp;

public:
	IMAGE(void);
	IMAGE(char* path);
	IMAGE(IMAGE& img);
	IMAGE(BITMAP& bmp);
	IMAGE( int w, int h );
	~IMAGE(void);

	// Basic Information
	BITMAPINFO Info(void);
	BITMAPINFO* pInfo(void);
	BITMAPFILEHEADER FileHeader(void);
	BITMAPFILEHEADER* pFileHeader(void);
	BITMAP BMP(void);
	BITMAP* pBMP(void);
	HBITMAP& Handle(void);
	bool Empty( void );
	void Resize( int w, int h );

	// Image Information
	BYTE& Pixel(int x, int y, BGR c);
	double Pixel(int x, int y, XYZ c);
	double Pixel(int x, int y, LUV c);
	int GrayPixel( int x, int y );
	bool SetPixel( int x, int y, BGR c, BYTE value );
	bool SetPixel( int x, int y, int bgr );
	int Width(void);
	int Height(void);

	// Image Process
	void Scale( int w, int h );

	// WIN32-API function
	bool LoadBMPFromFile( char* filename );
	bool WriteBMPFile( char *filename );
	int Display( HDC hdc, int x, int y, int width, int height, DWORD rop);
};

typedef IMAGE Image;

// File flow //
bool ReadBMP(char* path, BITMAP& img, BITMAPINFO& info, BITMAPFILEHEADER& fh);
bool ReadBMP(char* path, BITMAP& img);
bool WriteBMP(char* path, BITMAP& img, BITMAPINFO& info, BITMAPFILEHEADER& fh);
bool WriteBMP(char* path, BITMAP& img);

// Data Process //
HBITMAP CreateBMP( BITMAP *, int w, int h );
int DeleteBMP( BITMAP *, HBITMAP );

#endif