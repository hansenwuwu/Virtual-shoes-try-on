#include <GL/glut.h>
#include "scene.h"



// Initialze OpenGL perspective matrix
void GL_Setup()
{
	glShadeModel(GL_SMOOTH);							
	glClearColor(0.0f, 0.0f, 0.0f, 0.5f);			// Black Background
	glClearDepth(1.0f);					// Depth Buffer Setup
	glEnable(GL_DEPTH_TEST);				// Enables Depth Testing
	glDepthFunc(GL_LEQUAL);					// The Type Of Depth Testing To Do
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);	// Really Nice Perspective Calculations
}

void Light()
{
	float specular[4] = { .8f, .8f, .8f, 1.0f };
	float diffuse[4] = { .8f, .8f, .8f, 1.0f };
	float ambient[4] = { .2f, .2f, .2f, 1.0f };
	float position[4] = { 10.0f, 1000.0f, 0.0f, 1.0f };

	glLightfv( GL_LIGHT0, GL_SPECULAR, specular );
	glLightfv( GL_LIGHT0, GL_DIFFUSE, diffuse );
	glLightfv( GL_LIGHT0, GL_AMBIENT, ambient );
	glLightfv( GL_LIGHT0, GL_POSITION, position );

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
}
