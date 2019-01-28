#include "bitmap.h"
#include "imgproc.h"
#include <iostream>
#include <math.h>
#include <stdexcept>
using namespace std;

//
// class IMAGE
//

// constructor  

IMAGE::IMAGE(void){
	_bmp.bmBits = NULL;
	_hbmp = NULL;
}

IMAGE::IMAGE(char* filename){
	_bmp.bmBits = NULL;
	_hbmp = NULL;
	LoadBMPFromFile( filename );
}

IMAGE::IMAGE(IMAGE& img){
	// copy all the data
	_inf = img.Info();
	_fileheader = img.FileHeader();
	_bmp = img.BMP();

	// Create handle of bitmap
	_hbmp = CreateDIBSection( NULL, &_inf, DIB_RGB_COLORS, &_bmp.bmBits, NULL, 0 );
	SetDIBits( NULL, _hbmp, 0, _bmp.bmHeight, img.pBMP()->bmBits, &_inf, DIB_RGB_COLORS );
}

IMAGE::IMAGE(BITMAP& bmp){
	// construct the bitmapfileheader
	_fileheader.bfType = 0x4D42;
	_fileheader.bfOffBits = sizeof( BITMAPINFOHEADER );
	_fileheader.bfSize = 54 + bmp.bmHeight * bmp.bmWidthBytes;

	// construct the bitmapinfo
	_inf.bmiHeader.biSize = sizeof( BITMAPINFOHEADER );
	_inf.bmiHeader.biWidth = bmp.bmWidth;
	_inf.bmiHeader.biHeight = bmp.bmHeight;
	_inf.bmiHeader.biPlanes = bmp.bmPlanes;
	_inf.bmiHeader.biBitCount = bmp.bmBitsPixel;
	_inf.bmiHeader.biCompression = BI_RGB;
	_inf.bmiHeader.biSizeImage = 3 * bmp.bmHeight * bmp.bmWidth;
	_inf.bmiHeader.biXPelsPerMeter = 0;
	_inf.bmiHeader.biYPelsPerMeter = 0;
	_inf.bmiHeader.biClrUsed = 0;
	_inf.bmiHeader.biClrImportant = 0;
	_bmp = bmp;
	_hbmp = CreateDIBSection( NULL, &_inf, DIB_RGB_COLORS, &_bmp.bmBits, NULL, 0 );
	if( _hbmp != NULL ){
		SetDIBits( NULL, _hbmp, 0, bmp.bmHeight, bmp.bmBits, &_inf, DIB_RGB_COLORS );
	}
}

IMAGE::IMAGE( int w, int h ){
	// construct the bitmap
	_bmp.bmType = 0;
	_bmp.bmWidth = w;
	_bmp.bmHeight = h;
	_bmp.bmWidthBytes = w * 3;
	if( _bmp.bmWidthBytes % 4 != 0 ) _bmp.bmWidthBytes += ( 4 - ( _bmp.bmWidthBytes % 4 ) );
	_bmp.bmPlanes = 1;
	_bmp.bmBitsPixel = 24;

	// construct the bitmapfileheader
	_fileheader.bfType = 0x4D42;
	_fileheader.bfOffBits = sizeof( BITMAPINFOHEADER );
	_fileheader.bfSize = 54 + w * _bmp.bmWidthBytes;

	// construct the bitmapinfo
	_inf.bmiHeader.biSize = sizeof( BITMAPINFOHEADER );
	_inf.bmiHeader.biWidth = w;
	_inf.bmiHeader.biHeight = h;
	_inf.bmiHeader.biPlanes = _bmp.bmPlanes;
	_inf.bmiHeader.biBitCount = _bmp.bmBitsPixel;
	_inf.bmiHeader.biCompression = BI_RGB;
	_inf.bmiHeader.biSizeImage = 3 * h * w;
	_inf.bmiHeader.biXPelsPerMeter = 0;
	_inf.bmiHeader.biYPelsPerMeter = 0;
	_inf.bmiHeader.biClrUsed = 0;
	_inf.bmiHeader.biClrImportant = 0;
	_hbmp = CreateDIBSection( NULL, &_inf, DIB_RGB_COLORS, &_bmp.bmBits, NULL, 0 );
}

// destructor

IMAGE::~IMAGE(void){
	if( _hbmp != NULL ) DeleteObject( _hbmp );
}

// Basic Information

BITMAPINFO IMAGE::Info(void){
	return _inf;
}

BITMAPINFO* IMAGE::pInfo(void){
	return &_inf;
}

BITMAPFILEHEADER IMAGE::FileHeader(void){
	return _fileheader;
}

BITMAPFILEHEADER* IMAGE::pFileHeader(void){
	return &_fileheader;
}

BITMAP IMAGE::BMP(void){
	return _bmp;
}

BITMAP* IMAGE::pBMP(void){
	return &_bmp;
}

HBITMAP& IMAGE::Handle(void){
	return _hbmp;
}

void IMAGE::Resize( int w, int h ){
	HBITMAP hbmp = _hbmp;
	int h_old, w_old, stride_old;

	if( _bmp.bmBits == NULL ) {
		IMAGE( w, h );
		return;
	}

	h_old = _bmp.bmHeight;
	w_old = _bmp.bmWidth;
	stride_old = _bmp.bmWidthBytes;

	// update bitmap
	_bmp.bmHeight = h;
	_bmp.bmWidth = w;
	_bmp.bmWidthBytes = _bmp.bmWidth * 3;
	if( _bmp.bmWidthBytes % 4 != 0 ) _bmp.bmWidthBytes+= ( 4 - ( _bmp.bmWidthBytes % 4 ) );

	// update header information
	_inf.bmiHeader.biHeight = h;
	_inf.bmiHeader.biWidth = w;
	_inf.bmiHeader.biSizeImage = 3 * h * w;

	// update fileheader
	_fileheader.bfSize = 54 + _bmp.bmHeight * _bmp.bmWidthBytes;

	DeleteObject( _hbmp );
	_hbmp = CreateDIBSection( NULL, &_inf, DIB_RGB_COLORS, &_bmp.bmBits, NULL, 0 );

	if( _hbmp == NULL ) {
		_hbmp = hbmp;

		// bitmap
		_bmp.bmHeight = h_old;
		_bmp.bmWidth = w_old;
		_bmp.bmWidthBytes = _bmp.bmWidth * 3;
		if( _bmp.bmWidthBytes % 4 != 0 ) _bmp.bmWidthBytes+= ( 4 - ( _bmp.bmWidthBytes % 4 ) );

		// header information
		_inf.bmiHeader.biHeight = h_old;
		_inf.bmiHeader.biWidth = w_old;
		_inf.bmiHeader.biSizeImage = 3 * h_old * w_old;

		// fileheader
		_fileheader.bfSize = 54 + _bmp.bmHeight * _bmp.bmWidthBytes;
		return;
	}
}

bool IMAGE::Empty( void ){
	if( _bmp.bmBits == NULL ) return true;
	else return false;
}

// Image Information

BYTE& IMAGE::Pixel(int x, int y, BGR c){
	BYTE* output = (BYTE*)_bmp.bmBits;
	if( x < 0 || x >= _bmp.bmWidth || y < 0 || y >= _bmp.bmHeight )
		return output[0];
	return output[y * _bmp.bmWidthBytes + x * 3 + c];
}

double IMAGE::Pixel( int x, int y, XYZ c ){
	if( x < 0 || x >= _bmp.bmWidth || y < 0 || y >= _bmp.bmHeight )
		return false;
	double b, g, r;
	b = (double)Pixel( x, y, B ) / 255.0;
	g = (double)Pixel( x, y, G ) / 255.0;
	r = (double)Pixel( x, y, R ) / 255.0;
	switch( c ){
	case X: return ( 0.49 * r + 0.31 * g + 0.2 * b ) / 0.17697;
	case Y: return ( 0.17697 * r + 0.81240 * g + 0.01063 * b ) / 0.17697;
	case Z: return ( 0.00 * r + 0.01 * g + 0.88 * b ) / 0.17697;
	}
	return 0.0;
}

double IMAGE::Pixel( int x, int y, LUV c ){
	if( x < 0 || x >= _bmp.bmWidth || y < 0 || y >= _bmp.bmHeight )
		return false;
#define u (4*X/(Xc+15*Yc+3*Zc))
#define v (9*Y/(Xc+15*Yc+3*Zc))
#define un 0.2009
#define vn 0.4610
	double Xc, Yc, Zc, l;
	Xc = (double)Pixel( x, y, X );
	Yc = (double)Pixel( x, y, Y );
	Zc = (double)Pixel( x, y, Z );
	l = ( Yc > 0.00885645 ) ? ( 116.0 * pow( Yc, 0.333 ) - 16.0 ) : 903.2962963 * Yc;
	switch( c ){
	case L: return l;
	case U: return 13.0 * l * ( u - un );
	case V: return 13.0 * l * ( v - vn );
	}
	return 0.0;
#undef u
#undef v
#undef vn
#undef un
}

int IMAGE::GrayPixel( int x, int y ){
	if( x < 0 || x >= _bmp.bmWidth || y < 0 || y >= _bmp.bmHeight ) return false;
	BYTE* output = (BYTE*)_bmp.bmBits;
	return (
		(int)output[y * _bmp.bmWidthBytes + x * 3 + 0] +
		(int)output[y * _bmp.bmWidthBytes + x * 3 + 1] +
		(int)output[y * _bmp.bmWidthBytes + x * 3 + 2] ) / 3;
}

bool IMAGE::SetPixel( int x, int y, BGR c, BYTE v ){
	if( x < 0 || x >= _bmp.bmWidth || y < 0 || y >= _bmp.bmHeight )
		return false;
	else{
		BYTE *bits = ( BYTE* )_bmp.bmBits;
		bits[ y * _bmp.bmWidthBytes + 3 * x + c ] = v;
		return true;
	}
}

bool IMAGE::SetPixel( int x, int y, int bgr ){
	if( x < 0 || x >= _bmp.bmWidth || y < 0 || y >= _bmp.bmHeight )
		return false;
	else{
		BYTE *bits = ( BYTE* )_bmp.bmBits;
		bits[ y * _bmp.bmWidthBytes + 3 * x + 2 ] = bgr % 0x100 ;
		bgr /= 0x100;
		bits[ y * _bmp.bmWidthBytes + 3 * x + 1 ] = bgr % 0x100 ;
		bgr /= 0x100;
		bits[ y * _bmp.bmWidthBytes + 3 * x + 0 ] = bgr % 0x100 ;
		return true;
	}
}

int IMAGE::Width(void){
	return _bmp.bmWidth;
}

int IMAGE::Height(void){
	return _bmp.bmHeight;
}

// Image Process

void IMAGE::Scale( int w, int h ){
	HBITMAP hbmp = _hbmp;
	BYTE* pixel_new, *pixel_old;
	int h_old, w_old, stride_new, stride_old;

	if( _bmp.bmBits == NULL ) return;
	else pixel_old = (BYTE*)_bmp.bmBits;
	h_old = _bmp.bmHeight;
	w_old = _bmp.bmWidth;
	stride_old = _bmp.bmWidthBytes;

	// update bitmap
	_bmp.bmHeight = h;
	_bmp.bmWidth = w;
	_bmp.bmWidthBytes = _bmp.bmWidth * 3;
	if( _bmp.bmWidthBytes % 4 != 0 ) _bmp.bmWidthBytes+= ( 4 - ( _bmp.bmWidthBytes % 4 ) );

	// update header information
	_inf.bmiHeader.biHeight = h;
	_inf.bmiHeader.biWidth = w;
	_inf.bmiHeader.biSizeImage = 3 * h * w;
	_inf.bmiHeader.biXPelsPerMeter = 0;
	_inf.bmiHeader.biYPelsPerMeter = 0;

	// update fileheader
	_fileheader.bfSize = sizeof( BITMAPINFOHEADER ) + sizeof( BITMAPFILEHEADER ) + _bmp.bmWidthBytes * _bmp.bmHeight;

	_hbmp = CreateDIBSection( NULL, &_inf, DIB_RGB_COLORS, &_bmp.bmBits, NULL, 0 );
	if( _hbmp == NULL ) {
		_hbmp = hbmp;
		_bmp.bmBits = pixel_old;

		// bitmap
		_bmp.bmHeight = h_old;
		_bmp.bmWidth = w_old;
		_bmp.bmWidthBytes = _bmp.bmWidth * 3;
		if( _bmp.bmWidthBytes % 4 != 0 ) _bmp.bmWidthBytes+= ( 4 - ( _bmp.bmWidthBytes % 4 ) );

		// header information
		_inf.bmiHeader.biHeight = h_old;
		_inf.bmiHeader.biWidth = w_old;
		_inf.bmiHeader.biSizeImage = 3 * h_old * w_old;

		// fileheader
		_fileheader.bfSize = sizeof( BITMAPINFOHEADER ) + sizeof( BITMAPFILEHEADER ) + _bmp.bmWidthBytes * _bmp.bmHeight;
		return;
	}

	stride_new = _bmp.bmWidthBytes;
	pixel_new = (BYTE*) _bmp.bmBits;

	::Scale( pixel_old, w_old, h_old, stride_old, pixel_new, w, h, stride_new );
	DeleteObject( hbmp );
}
	

// API Finction

bool IMAGE::LoadBMPFromFile( char *filename ){
	LPVOID bits;
	if( _bmp.bmBits != NULL ) {
		DeleteObject( _hbmp );
		_bmp.bmBits = NULL;
	}
	if( !ReadBMP( filename, _bmp, _inf, _fileheader ) ){	// if file dosn't exist
		_bmp.bmBits = NULL;									// set as NULL
		return false;
	}
	bits = _bmp.bmBits;
	_hbmp = CreateDIBSection( NULL, &_inf, DIB_RGB_COLORS, &_bmp.bmBits, NULL, 0 );	// using bmBits to Create DIB
	if( _hbmp == NULL ) return false;
	SetDIBits( NULL, _hbmp, 0, _bmp.bmHeight, bits, &_inf, DIB_RGB_COLORS );
	free( bits );
	return true;
}

bool IMAGE::WriteBMPFile( char *file ){
	return WriteBMP( file, _bmp );
}

int IMAGE::Display( HDC hdc, int x, int y, int w, int h, DWORD rop ){
	int re;
	HDC hdc_mem = CreateCompatibleDC( hdc );
	SelectObject( hdc_mem, _hbmp );
	re = StretchBlt( hdc, x, y, w, h, hdc_mem, 0, 0, _bmp.bmWidth, _bmp.bmHeight, rop );
	SelectObject( hdc_mem, NULL );
	DeleteDC( hdc_mem );
	return re;
}

// 
// File Flow
//

bool ReadBMP( char* path, BITMAP& img, BITMAPINFO& info, BITMAPFILEHEADER& fh ){
	FILE* infile;
	WORD buf;
	fopen_s(&infile, path, "rb");
	if(infile==NULL){
		cerr<< "can't open "<< path<< endl;
		return false;
	}
	fread(&fh.bfType, 2, 1, infile);
	if(fh.bfType!=0x4D42){
		cerr<< "it's not a bmp file.\n";
		return false;
	}
	img.bmType=0;
	fread(&fh.bfSize, 4, 1, infile);	//file size
	info.bmiHeader.biSize = sizeof( BITMAPINFOHEADER );
	fread(&fh.bfReserved1, 2, 1, infile);	//reserved
	fread(&fh.bfReserved2, 2, 1, infile);	//reserved
	fread(&fh.bfOffBits, 4, 1, infile);	//image offset
	fread(&buf, 4, 1, infile);	//header size
	fread(&info.bmiHeader.biWidth, 4, 1, infile);	//image DX
	fread(&info.bmiHeader.biHeight, 4, 1, infile);	//imgae DY
	img.bmHeight=info.bmiHeader.biHeight;
	img.bmWidth=info.bmiHeader.biWidth;
	img.bmWidthBytes=img.bmWidth*3;
	if( img.bmWidthBytes % 4 != 0 ) img.bmWidthBytes+= ( 4 - ( img.bmWidthBytes % 4 ) );
	fread(&info.bmiHeader.biPlanes, 2, 1, infile);	//bit plane
	img.bmPlanes=info.bmiHeader.biPlanes;
	fread(&info.bmiHeader.biBitCount, 2, 1, infile);	//bit count
	img.bmBitsPixel=info.bmiHeader.biBitCount;
	fread(&info.bmiHeader.biCompression, 4, 1, infile);	//compression, 0=RGB, 1=RLE8, 2=RLE4
	fread(&info.bmiHeader.biSizeImage, 4, 1, infile);	//image size
	fread(&info.bmiHeader.biXPelsPerMeter, 4, 1, infile);	//XDPM
	fread(&info.bmiHeader.biYPelsPerMeter, 4, 1, infile);	//YDPM
	fread(&info.bmiHeader.biClrUsed, 4, 1, infile);	//Color Used
	fread(&info.bmiHeader.biClrImportant, 4, 1, infile);	//color important
	info.bmiHeader.biSizeImage = 3 * info.bmiHeader.biHeight * info.bmiHeader.biWidth;

	if( img.bmBits != NULL ) free( img.bmBits );
	img.bmBits=malloc(img.bmHeight*img.bmWidthBytes*sizeof(BYTE));
	fread(img.bmBits, 1, img.bmHeight*img.bmWidthBytes, infile);
	fclose(infile);
	return true;
}
bool ReadBMP(char* path, BITMAP& img){
	FILE* infile;
	WORD buf;
	unsigned int i;
	fopen_s(&infile, path, "rb");
	if(infile==NULL){
		cerr<< "can't open "<< path<< endl;
		return false;
	}
	fread(&buf, 2, 1, infile);
	if(buf!=0x4D42){
		cerr<< "it's not a bmp file.\n";
		return false;
	}
	img.bmType=0;
	for( i=0; i<4; i++) fread(&buf, 4, 1, infile);	//buffer
	fread(&img.bmWidth, 4, 1, infile);	//image DX
	fread(&img.bmHeight, 4, 1, infile);	//imgae DY
	img.bmWidthBytes=img.bmWidth*3;
	if( img.bmWidthBytes % 4 != 0 ) img.bmWidthBytes += ( 4 - ( img.bmWidthBytes % 4 ) );
	fread(&img.bmPlanes, 2, 1, infile);	//bit plane
	fread(&img.bmBitsPixel, 2, 1, infile);	//bit count
	for(i=0; i<6; i++) fread(&buf, 4, 1, infile);	//buffer

	if( img.bmBits != NULL ) free( img.bmBits );
	img.bmBits=malloc(3*img.bmHeight*img.bmWidthBytes*sizeof(BYTE));
	fread(img.bmBits, 1, 3*img.bmHeight*img.bmWidthBytes, infile);
	fclose(infile);
	return true;
}
bool WriteBMP(char* path, BITMAP& img, BITMAPINFO& info, BITMAPFILEHEADER& fh){
	FILE* outfile;
	WORD buf=40;
	fopen_s(&outfile, path, "wb");
	if(outfile==NULL){
		cerr<< "can't open "<< path<< endl;
		return false;
	}
	fwrite(&fh.bfType, 2, 1, outfile);	//BM
	fwrite(&fh.bfSize, 4, 1, outfile);	//file size
	fwrite(&fh.bfReserved1, 2, 1, outfile);	//reserved
	fwrite(&fh.bfReserved2, 2, 1, outfile);	//reserved
	fwrite(&fh.bfOffBits, 4, 1, outfile);	//image offset
	fwrite(&buf, 4, 1, outfile);	//header size
	fwrite(&info.bmiHeader.biWidth, 4, 1, outfile);	//image DX
	fwrite(&info.bmiHeader.biHeight, 4, 1, outfile);	//imgae DY
	fwrite(&info.bmiHeader.biPlanes, 2, 1, outfile);	//bit plane
	fwrite(&info.bmiHeader.biBitCount, 2, 1, outfile);	//bit count
	fwrite(&info.bmiHeader.biCompression, 4, 1, outfile);	//compression, 0=RGB, 1=RLE8, 2=RLE4
	fwrite(&info.bmiHeader.biSizeImage, 4, 1, outfile);	//image size
	fwrite(&info.bmiHeader.biXPelsPerMeter, 4, 1, outfile);	//XDPM
	fwrite(&info.bmiHeader.biYPelsPerMeter, 4, 1, outfile);	//YDPM
	fwrite(&info.bmiHeader.biClrUsed, 4, 1, outfile);	//Color Used
	fwrite(&info.bmiHeader.biClrImportant, 4, 1, outfile);	//color important
	fwrite(img.bmBits, 1, img.bmHeight*img.bmWidthBytes, outfile);
	fclose(outfile);
	return true;
}
bool WriteBMP(char* path, BITMAP& img){
	FILE* outfile;
	DWORD buf;
	fopen_s(&outfile, path, "wb");
	if(outfile==NULL){
		cerr<< "can't open "<< path<< endl;
		return false;
	}
	fwrite(&(buf=0x4D42), 2, 1, outfile);
	fwrite(&(buf=img.bmWidthBytes*img.bmHeight+54), 4, 1, outfile);	//file size
	fwrite(&(buf=0x00000000), 4, 1, outfile);	//reserved
	fwrite(&(buf=54), 4, 1, outfile);	//image offset
	fwrite(&(buf=40), 4, 1, outfile);	//header size
	fwrite(&img.bmWidth, 4, 1, outfile);	//image DX
	fwrite(&img.bmHeight, 4, 1, outfile);	//imgae DY
	fwrite(&img.bmPlanes, 2, 1, outfile);	//bit plane
	fwrite(&img.bmBitsPixel, 2, 1, outfile);	//bit count
	fwrite(&(buf=0), 4, 1, outfile);	//compression, 0=RGB, 1=RLE8, 2=RLE4
	fwrite(&(buf=img.bmWidth*img.bmHeight*3), 4, 1, outfile);	//image size
	fwrite(&(buf=0x00000000), 4, 1, outfile);	//XDPM
	fwrite(&(buf=0x00000000), 4, 1, outfile);	//YDPM
	fwrite(&(buf=0x00000000), 4, 1, outfile);	//Color Used
	fwrite(&(buf=0x00000000), 4, 1, outfile);	//color important
	fwrite(img.bmBits, 1, img.bmHeight*img.bmWidthBytes, outfile);
	fclose(outfile);
	return true;
}

//
// Data Process //
//
BITMAP CopyBITMAP(BITMAP& img){
	BITMAP output=img;
	BYTE* pixel=(BYTE*)malloc(img.bmHeight*img.bmWidthBytes*sizeof(BYTE));
	BYTE* o_pix=(BYTE*)img.bmBits;
	for(int i=0; i<img.bmHeight*img.bmWidthBytes; i++) pixel[i]=o_pix[i];
	output.bmBits=(BYTE*)pixel;
	return output;
}

BYTE* CopyPixel(LPVOID img, int h, int stride){
	BYTE* output=(BYTE*)malloc(h*stride*sizeof(BYTE));
	BYTE* pixel=(BYTE*)img;
	for(int i=0; i<h*stride; i++) output[i]=pixel[i];
	return output;
}
/*
void InitBMP24Bits( BITMAP *bmp, int w, int h ){
	bmp->bmBitsPixel = 24;
	bmp->bmHeight = h;
	bmp->bmWidth = w;
	bmp->bmPlanes = 1;
	bmp->bmType = 0;
	bmp->bmWidthBytes = ( (bmp->bmWidth*3) % 4 == 0 ) ? (bmp->bmWidth*3) : ( (bmp->bmWidth*3) + 4 - ( (bmp->bmWidth*3) % 4 ) );
	bmp->bmBits = malloc( 3 * bmp->bmWidthBytes * bmp->bmHeight );
}
/**/
HBITMAP CreateBMP( BITMAP *bmp, int w, int h ){
	HBITMAP hbmp;
	BITMAPINFO bmpinfo;

	// construct the bitmapinfo
	bmpinfo.bmiHeader.biSize = sizeof( BITMAPINFOHEADER );
	bmpinfo.bmiHeader.biWidth = w;
	bmpinfo.bmiHeader.biHeight = h;
	bmpinfo.bmiHeader.biPlanes = 1;
	bmpinfo.bmiHeader.biBitCount = 24;
	bmpinfo.bmiHeader.biCompression = BI_RGB;
	bmpinfo.bmiHeader.biSizeImage = 3 * h * w;
	bmpinfo.bmiHeader.biXPelsPerMeter = 0;
	bmpinfo.bmiHeader.biYPelsPerMeter = 0;
	bmpinfo.bmiHeader.biClrUsed = 0;
	bmpinfo.bmiHeader.biClrImportant = 0;

	// bitmap structure
	bmp->bmWidth = w;
	bmp->bmHeight = h;
	bmp->bmWidthBytes = w * 3;
	if( bmp->bmWidthBytes % 4 ) bmp->bmWidthBytes += 4 - bmp->bmWidthBytes % 4;
	bmp->bmPlanes = 1;
	bmp->bmBitsPixel = 24;

	// Create
	hbmp = CreateDIBSection( NULL, &bmpinfo, DIB_RGB_COLORS, &(*bmp).bmBits, NULL, 0 );
	if( hbmp != NULL ){
		SetDIBits( NULL, hbmp, 0, (*bmp).bmHeight, (*bmp).bmBits, &bmpinfo, DIB_RGB_COLORS );
	}
	else {
		OutputDebugString( L"can't create BITMAP image\n" );
		(*bmp).bmBits = NULL;
	}
	return hbmp;
}

int DeleteBMP( BITMAP *bmp, HBITMAP hbmp ){
	if( hbmp ) {
		DeleteObject( hbmp );
		return 1;
	}
	else return 0;
}