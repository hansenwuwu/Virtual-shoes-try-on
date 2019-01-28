#include "ppm.h"
#include <stdlib.h>
#include <stdio.h>

int ReadPPM( BYTE **img, int *w, int *h, char *path ){
	char string[1000];
	FILE *infile;
	int i, n;

	// open file
	fopen_s( &infile, path, "rb" );
	if( infile == NULL ) {
		OutputDebugString( L"can't open PPM file" );
		return 1;
	}

	// check file format
	fscanf_s( infile, "%s", string );
	if( string[0] != 'P' || ( string[1] != '3' && string[1] != '6' ) ){
		OutputDebugString( L"it's not a PPM file" );
		fclose( infile );
		return 1;
	}

	fscanf_s( infile, "%d", w );
	fscanf_s( infile, "%d", h );
	*img = ( BYTE* )malloc( 3 * ( *w ) * ( *h ) * sizeof( BYTE ) );
	if( string[1] == '6' )
		fread( *img, 3 * ( *w ) * ( *h ), sizeof( BYTE ), infile );
	else
		for( i = 0, n = 3 * ( *w ) * ( *h ); i < n; i++ ) fscanf_s( infile, "%c", img[i] );
	
	fclose( infile );
	return 0;
}

int ReadPPM( double **img, int *w, int *h, char *path ){
	BYTE *bits, *bit;
	double *pix;
	int i, j;

	i = ReadPPM( &bits, w, h, path );
	if( i ) return i;


	pix = *img;
	for( i = *h - 1; i >= 0; i-- ){
		bit = bits + i * w[0] * 3;
		for( j = 0; j < *w; j++, bit += 3, pix += 3 ) {
			pix[0] = bit[0];
			pix[1] = bit[1];
			pix[2] = bit[2];
		}
	}
	return 0;
}

int ReadPPM( BITMAP *img, char *path ){
	int w, h, h2, i, j;
	BYTE *bit1, *bit2, temp;

	i = ReadPPM( ( BYTE** )( &img->bmBits ), &w, &h, path );
	if( i ) return i;

	img->bmWidth = w;
	img->bmHeight = h;
	img->bmWidthBytes = img->bmWidth * 3;
	img->bmBitsPixel = 24;
	img->bmPlanes = 1;
	img->bmType = 0;

	// variefy the order in space
	h2 = h / 2;
	bit1 = ( BYTE* )img->bmBits;
	for( i = 0; i < h2; i++ ){
		bit2 = ( BYTE* )img->bmBits + ( h - i - 1 ) * w * 3;
		for( j = 0; j < w; j++, bit1 += 3, bit2 += 3 ) {
			temp = bit1[0];	bit1[0] = bit2[0]; bit2[0] = temp;
			temp = bit1[1];	bit1[1] = bit2[1]; bit2[1] = temp;
			temp = bit1[2];	bit1[2] = bit2[2]; bit2[2] = temp;
		}
	}
/*
	// change RGB to BGR
	bit1 = ( BYTE* )img->bmBits;
	h2 = h * w;
	for( i = 0; i < h2; i++, bit1 += 3 ){
		temp = bit1[0];
		bit1[0] = bit1[2];
		bit1[2] = temp;
	}
/**/
	return 0;
}