#include "../stdafx.h"
#include "scene.h"
#include "imgmap.h"
#include "../include/GL/glut.h"



void DrawBackgroundUC( color3uc *_color)
{
	glDisable( GL_DEPTH_TEST );
	//glPixelZoom( (float)2160/1920.0, (float)1215/1080.0 );
	//glPixelZoom((float)2160 / 1388.0, (float)1680 / 1080.0);
	glPixelZoom(1.0,1.0);
	//glRasterPos2f( -1.f, 1.f );
	glPixelStorei( GL_PACK_ALIGNMENT, 4 );
	glDrawPixels(1920, 1080, GL_RGB, GL_UNSIGNED_BYTE, _color);
	glEnable( GL_DEPTH_TEST );	
}

void DrawDepthPointcloud( color3uc *_color, point3f *_depth)
{
	glBegin( GL_POINTS );
	int w, h;
	w = 1920;
	h = 1080;
	
	for (int i = 0; i < w * h; i++){
		glColor3f( (float)_color[i].r/255, (float)_color[i].g/255, (float)_color[i].b/255 );	
		glVertex3f( _depth[i].x/1000, _depth[i].y/1000, _depth[i].z/1000 );
	}
	glEnd();
}
