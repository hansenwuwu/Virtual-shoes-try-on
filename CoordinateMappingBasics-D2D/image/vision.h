#ifndef _RICE_COMPUTER_VISION_H_
#define _RICE_COMPUTER_VISION_H_

//
// structures
//

typedef struct _sequence_buffer{
	void **seq;
	int w;
	int h;
	int stride;
	int color;
	int n;
} SequenceBuffer;

typedef struct _optical_flow{
	float *u;
	float *v;
	float *buf;
	float alpha;
	float iter;
} OpticalFlow;

//
// Macro
//

#define OF_HORN_ALGO	1
#define OF_LK_ALGO		2
typedef int				OptFlowAlgo;
//typedef int( *OptFlowAlgo )( OpticalFlow *, const SequenceBuffer );

//
// Functions
//

// Memory allocate
int NewOpticalFlow( OpticalFlow *, const SequenceBuffer, OptFlowAlgo );
int DeleteOpticalFlow( OpticalFlow * );

// Motion analysis

int HornOptFlowi( OpticalFlow *, const SequenceBuffer );
int HornOptFlowb( OpticalFlow *, const SequenceBuffer );
int HornOptFlowf( OpticalFlow *, const SequenceBuffer );
int HornOptFlowd( OpticalFlow *, const SequenceBuffer );

int LKOptFlowi( OpticalFlow *, const SequenceBuffer );
int LKOptFlowb( OpticalFlow *, const SequenceBuffer );
int LKOptFlowf( OpticalFlow *, const SequenceBuffer );
int LKOptFlowd( OpticalFlow *, const SequenceBuffer );

#endif