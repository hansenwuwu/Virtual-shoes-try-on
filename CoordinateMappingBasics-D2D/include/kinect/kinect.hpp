
#ifndef __KINECT__
#define __KINECT__

#include <XnCppWrapper.h>
#include <iostream>
#include <string>
#include <math.h>
#include <gl/glew.h>
#include <scene.h>

#define		NONE		0
#define		ANY			1
#define		ALL			2


class KINECT
{
	public:
		/* create object
		_________________________________________________________*/

			// use sensor to play
			KINECT( int width, int height, int fps );
			
			// load *.oni file to replay record file
			KINECT( int width, int height, int fps, const std::string file_name );
			
			// default width=640, height=480, fps=30
			KINECT( void );
			~KINECT( void );
			int width;			// camera width
			int height;			// camera height
			int size;			// width * height

		/* control device
		_________________________________________________________*/

			// start device ( must use this before get data )
			bool Start( void );

			// stop device
			bool Stop( void );

			// update data by NONE, ANY, or ALL
			bool Update( int wait_type );
			
			// update data by specific node
			bool Update( xn::ProductionNode &type );

		/* display
		_________________________________________________________*/
			
			// set depth detection distance
			bool SetViewRange( int view_distance_in_mm );
			int range;			// depth distance

		/* data
		_________________________________________________________*/

			// store depth data to map
			void GetDepth( point *depth );
			void GetRawDepth( point *depth );
			void GetBothDepth( point *project_depth, point *real_world_depth );

			// convert point between projective and realworld
			void ConvertProjectToRealWorld( int point_number, const point *input, point *output );
			void ConvertRealWorldToProject( int point_number, const point *input, point *output );

			// store color data to map with arrange type
			void GetImage( color *image );

		/* record
		_________________________________________________________*/

			// start record frame with *.oni file to save
			bool StartRecord( const std::string sFile );

			// stop record
			void StopRecord( void );

		/* gesture
		_________________________________________________________*/

		/* openGL
		_________________________________________________________*/

			// set openGL viewing ( must use this when first Start device )
			friend static bool K_InitGLView( void );
			double view_l;			// left, right, top, and bottom of display
			double view_r;
			double view_t;
			double view_b;
			double nearClip;		// near & far clipping, near is also the focal length
			double farClip;

	protected:
		
		/* error
		_________________________________________________________*/

			// check openNI errors
			//bool K_CheckError( XnStatus result, const std::string state );
			bool kinect_error;				// store error

		/* settings
		_________________________________________________________*/

			// set device width, height, and fps
			bool Setting( int width, int height, int fps );
			XnMapOutputMode mapMode;		// display mode

			// initialize device
			bool Init( void );
			xn::Context context;			// device context
			bool InitDepth( void );
			xn::DepthGenerator depthGener;	// depth generator
			bool InitImage( void );
			xn::ImageGenerator imageGener;	// image generator

			xn::Recorder recGener;			// record generator

		/* data
		_________________________________________________________*/

			// convert source point to use with arrange type
			const XnDepthPixel *XN_DEPTH;	// pointer to depth
			const XnRGB24Pixel *XN_RGB;		// pointer to image color

		/* temperate
		_________________________________________________________*/

			int run_i, run_j;				// runner for GetDepth and GetImage

		/* specific member
		_________________________________________________________*/

			float rgb_trans;
};


/*===============================================================
	Do
===============================================================*/

bool K_CheckError( XnStatus _r, const std::string _s )
{
	if( _r != XN_STATUS_OK )
	{
		std::cout << "[ERROR] " << _s << "\n\t" << xnGetStatusString( _r ) << std::endl;
		std::cout << "Press any key to continue...";
		system( "pause" );
		exit( 1 );
		return 1;
	}
	return 0;
}


KINECT::KINECT( void )
{
	if( !Setting( 640, 480, 30 ) )
		if( Init() )
			width = height = -1;
}

KINECT::KINECT( int _w, int _h, int _fps )
{
	if( !Setting( _w, _h, _fps ) )
		if( Init() )
			width = height = -1;
}

KINECT::KINECT( int _w, int _h, int _fps, const std::string _fname )
{
	if( !Setting( _w, _h, _fps ) )
		if( !K_CheckError( context.OpenFileRecording( _fname.c_str() ), "open *.oni file" ) )
			if( Init() )
				width = height = -1;
}

KINECT::~KINECT( void )
{
}

bool KINECT::Setting( int _w, int _h, int _fps )
{
	rgb_trans = 1.f/255.f;
	if( K_CheckError( context.Init(), "initialize device" ) )
		return 1;

	if( _w <= 0 || _h <= 0 )
		return 1;

	width = _w;
	height = _h;
	size = width*height;

	mapMode.nXRes = width;
	mapMode.nYRes = height;
	mapMode.nFPS = _fps;

	SetViewRange( 5000 );

	return 0;
}

bool KINECT::Init( void )
{
	kinect_error = 0;
	kinect_error |= InitImage();
	kinect_error |= InitDepth();
	kinect_error |= K_CheckError( depthGener.GetAlternativeViewPointCap().SetViewPoint( imageGener ), "alternative view points" );
	return kinect_error;
}

bool KINECT::InitDepth( void )
{
	kinect_error = 0;
	kinect_error |= K_CheckError( depthGener.Create( context ), "initialize depth sensor" );

	/* the following functions may only work on Kinect for Windows rather than for Xbox */
	//kinect_error |= K_CheckError( depthGener.SetIntProperty( "nearMode", 1 ), "turn near mode on" );
	//kinect_error |= K_CheckError( depthGener.SetIntProperty( "distinctOverflowDepthValues", 1 ), "turn distinct overflow depth mode on" );
	
	kinect_error |= K_CheckError( depthGener.SetMapOutputMode( mapMode ), "set depth mode" );
	return kinect_error;
}

bool KINECT::InitImage( void )
{
	kinect_error = 0;
	kinect_error |= K_CheckError( imageGener.Create( context ), "initialize image sensor" );
	kinect_error |= K_CheckError( imageGener.SetMapOutputMode( mapMode ), "set image mode" );
	return kinect_error;
}

bool KINECT::StartRecord( const std::string _fname )
{
	kinect_error = 0;
	kinect_error |= K_CheckError( recGener.Create( context ), "initialize record device" );
	kinect_error |= K_CheckError( recGener.SetDestination( XN_RECORD_MEDIUM_FILE, _fname.c_str() ), "set store file" );

	// compression type: XN_CODEC_NULL, XN_CODEC_UNCOMPRESSED, XN_CODEC_JPEG
	kinect_error |= K_CheckError( recGener.AddNodeToRecording( imageGener, XN_CODEC_UNCOMPRESSED ), "set image store type" );

	// compression type: XN_CODEC_NULL, XN_CODEC_UNCOMPRESSED, XN_CODEC_16Z, XN_CODEC_16Z_EMB_TABLES, or XN_CODEC_8Z
	kinect_error |= K_CheckError( recGener.AddNodeToRecording( depthGener, XN_CODEC_UNCOMPRESSED ), "set depth store type" );

	return kinect_error;
}

void KINECT::StopRecord( void )
{
	recGener.Release();
}

bool KINECT::Start( void )
{
	kinect_error = 0;
	kinect_error |= K_CheckError( context.StartGeneratingAll(), "start device" );
	kinect_error |= Update( ALL );
	return kinect_error;
}

bool KINECT::Stop( void )
{
	kinect_error = 0;
	kinect_error |= K_CheckError( context.StopGeneratingAll(), "stop device" );
	context.Release();
	return kinect_error;
}

bool KINECT::Update( int _type )
{
	switch( _type )
	{
		default:
		case NONE:
			context.WaitNoneUpdateAll();
			return 0;

		case ANY:
			context.WaitAnyUpdateAll();
			return 0;

		case ALL:
			context.WaitAndUpdateAll();
			return 0;
	}
}

bool KINECT::Update( xn::ProductionNode &_type )
{
	if( K_CheckError( context.WaitOneUpdateAll( _type ), "wait one update" ) )
		return 1;
	return 0;
}

void KINECT::ConvertProjectToRealWorld( int _num, const point *_in, point *_out )
{
	depthGener.ConvertProjectiveToRealWorld( _num, (XnPoint3D*)_in, (XnPoint3D*)_out );
}

void KINECT::ConvertRealWorldToProject( int _num, const point *_in, point *_out )
{
	depthGener.ConvertRealWorldToProjective( _num, (XnPoint3D*)_in, (XnPoint3D*)_out );
}

void KINECT::GetRawDepth( point *_map )
{
	XN_DEPTH = depthGener.GetDepthMap();
	for( run_i = height; run_i > 0; run_i-- )
	{
		for( run_j = width; run_j > 0; run_j--, _map++, XN_DEPTH++ )
		{
			_map[0].x = (float)run_j;
			_map[0].y = (float)run_i;
			if( XN_DEPTH[0] > range )
				_map[0].z = 0.f;
			else
				_map[0].z = (float)-XN_DEPTH[0];
		}
	}
}

void KINECT::GetDepth( point *_map )
{
	GetRawDepth( _map );
	depthGener.ConvertProjectiveToRealWorld( size, (XnPoint3D*)_map, (XnPoint3D*)_map );
}

void KINECT::GetBothDepth( point *_pj, point *_rw )
{
	GetRawDepth( _pj );
	depthGener.ConvertProjectiveToRealWorld( size, (XnPoint3D*)_pj, (XnPoint3D*)_rw );
}

void KINECT::GetImage( color *_map )
{
	XN_RGB = imageGener.GetRGB24ImageMap();
	for( run_i = 0; run_i < height; run_i++ )
	{
		for( run_j = 0; run_j < width; run_j++, XN_RGB++, _map++ )
		{
			_map[0].r = XN_RGB[0].nRed*rgb_trans;
			_map[0].g = XN_RGB[0].nGreen*rgb_trans;
			_map[0].b = XN_RGB[0].nBlue*rgb_trans;
		}
	}
}

bool KINECT::SetViewRange( int _r )
{
	if( _r > 0 )
	{
		range = _r;
		return 0;
	}
	range = -_r;
	return 1;
}


#endif