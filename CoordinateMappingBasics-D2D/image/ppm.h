#ifndef _PPM_H_
#define _PPM_H_

#include <windows.h>

int ReadPPM( BYTE **img, int *w, int *h, char *path );
int ReadPPM( double **img, int *w, int *h, char *path );
int ReadPPM( BITMAP *bmp, char *path );

#endif