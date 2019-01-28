
#ifndef __AR_TOOLKIT_EX__
#define __AR_TOOLKIT_EX__

#include "ar_toolkit.hpp"
#include <scene.h>

#define		ATE_ROUND( X )		(int)( (float)X + .5f )

point2i		ATE_SIDE[4]			;

void		ATE_DrawCoordinate( void );
void		ATE_SetSidePoint( ARMarkerInfo & /* marker info */ );
bool		ATE_CheckPointInSide( int /* x */, int /* y */ );


/*	Do
________________________________________________________________*/

void ATE_DrawCoordinate( void )
{
	glBegin( GL_LINES );
		glColor4f( 1.f, 0.f, 0.f, .5f );
		glVertex3f( 0.f, 0.f, 0.f );
		glVertex3f( 1.f, 0.f, 0.f );
		glColor4f( 0.f, 1.f, 0.f, .5f );
		glVertex3f( 0.f, 0.f, 0.f );
		glVertex3f( 0.f, 1.f, 0.f );
		glColor4f( 0.f, 0.f, 1.f, .5f );
		glVertex3f( 0.f, 0.f, 0.f );
		glVertex3f( 0.f, 0.f, 1.f );
	glEnd();
}

void ATE_SetSidePoint( ARMarkerInfo &_info )
{
	static int i;
	for( i = 0; i < 4; i++ )
		ATE_SIDE[i] = point2i( ATE_ROUND( _info.vertex[i][0] ), ATE_ROUND( _info.vertex[i][1] ) );
}

bool ATE_CheckPointInSide( int _x, int _y )
{
	// fast algorithm
	static bool cross;
	static int i, j;

	cross = false;
	for( i = 0, j = 3; i < 4; j = i++ )
	{
		if( ( ATE_SIDE[i].y > _y ) != ( ATE_SIDE[j].y > _y ) && 
			_x < ( ATE_SIDE[j].x - ATE_SIDE[i].x )*( _y - ATE_SIDE[i].y )/( ATE_SIDE[j].y - ATE_SIDE[i].y ) + ATE_SIDE[i].x )
		{
			cross = !cross;
		}
	}
	return cross;

	/*
	static int cross, i, x;
	static point2i p1, p2;

	cross = 0;

	for( i = 0; i < 4; i++ )
	{
		p1 = ATE_SIDE[i];
		p2 = ATE_SIDE[ ( i + 1 )%4 ];

		// solve intersection point between ( y = _y ) and ( p1p2 )
			if( p1.y == p2.y )
				// parallel
				continue;
			if( _y < min( p1.y, p2.y ) )
				// intersection point on the extend line of ( p1p2 ) 
				continue;
			if( _y >= max( p1.y, p2.y ) )
				// intersection point on the extend line of ( p1p2 )
				continue;

			// find the x of intersection point
			x = ( _y - p1.y )*( p2.x - p1.x )/( p2.y - p1.y ) + p1.x;
			if( x > _x )
				cross++;
	}
	return ( cross%2 == 1 );
	*/
}



#endif