#ifndef _FOOT_TRACKING_CONFIG_H_
#define _FOOT_TRACKING_CONFIG_H_

#include "../image/image.h"
#include "../lib3D/depthimage.h"
#include "../memory/mem.h"

//
// inline function
//

inline bool IsGreen( int *color ){
	return ( color[2] < -40 && abs( color[1] ) + abs( color[1] ) < abs( color[2] ) );
}


inline bool IsRed( int *color ){
	return ( color[1] > 80 && abs( color[2] ) + abs( color[2] ) < abs( color[1] ) );
}

inline bool IsSkin( int *color){
	//return ( color[1] <- 40 &&  abs( color[1] ) < abs( color[2] ) );
	//return (color[1]<0 && color[2]>0 );
	//return 1;
	//return ( color[1] > 50 &&  abs( color[2] ) < abs( color[1] ) );
	return ( color[1] < 30  &&  color[1] > 20 &&  color[2] < 15 &&  color[2] > -50 );
//	return ( color[2] < -40 &&  abs( color[1] ) < abs( color[2] ) && color[1]<0);

}

#endif
