
#ifndef __AR_TOOLKIT__
#define __AR_TOOLKIT__

#include <ar/gsub.h>
#include <ar/ar.h>
#include <string>

const char		*AT_CAM_PARA_FILE				= "patt/camera_para.dat";
const char		*AT_PATT_CALIB_FILE				= "patt/patt.calib";
const char		*AT_PATT_HIRO_FILE				= "patt/patt.hiro";
const char		*AT_PATT_KANJI_FILE				= "patt/patt.kanji";
const char		*AT_PATT_SAMPLE1_FILE			= "patt/patt.sample1";
const char		*AT_PATT_SAMPLE2_FILE			= "patt/patt.sample2";

#define			AT_PATT_CALIB					AT_PATT_CALIB_FILE
#define			AT_PATT_HIRO					AT_PATT_HIRO_FILE
#define			AT_PATT_KANJI					AT_PATT_KANJI_FILE
#define			AT_PATT_SAMPLE1					AT_PATT_SAMPLE1_FILE
#define			AT_PATT_SAMPLE2					AT_PATT_SAMPLE2_FILE

class AT_PATTERN
{
	public:
		AT_PATTERN( const char * /* marker file */ );
		int id;
		double marker_width;
		double center[2];
		double opengl_trans[16];
		double trans[3][4];
		std::string error;
};

AT_PATTERN::AT_PATTERN( const char *_pattfile )
{
	if( ( id = arLoadPatt( _pattfile ) ) < 0 )
		error = "patten file error";
	center[0] = center[1] = 0.0;
	marker_width = 80.0;
}

ARParam			AT_WP							;
ARParam			AT_CP							;
int				AT_MARKER_NUM					;
int				AT_THRESH						= 100;
ARMarkerInfo	*AT_MARKER_INFO					;
std::string		AT_ERROR						= "";

bool	AT_Init				( void );
bool	AT_Init				( int /* camera width */, int /* camera height */, const char * /* camera parameter file */ cam_para = AT_CAM_PARA_FILE );
int		AT_SetThresh		( int /* thresh value */ );
int		AT_Detect			( ARUint8 * /* RGBA */ );
int		AT_Detect			( ARUint8 * /* RGBA */, int /* thresh */ );
int		AT_Detect			( ARUint8 * /* RGBA */, int /* thresh */, ARMarkerInfo ** /* marker info */, int & /* number of detected markers */ );
int		AT_Detect			( ARUint8 * /* RGBA */, ARMarkerInfo ** /* marker info */, int & /* number of detected markers */ );
int		AT_Detect			( ARUint8 * /* RGBA */, AT_PATTERN & /* pattern */ );
int		AT_Detect			( ARUint8 * /* RGBA */, AT_PATTERN & /* pattern */, int /* thresh */ );
int		AT_Detect			( ARUint8 * /* RGBA */, AT_PATTERN & /* pattern */, int /* thresh */, ARMarkerInfo ** /* marker info */, int & /* number of detected markers */ );
int		AT_Detect			( ARUint8 * /* RGBA */, AT_PATTERN & /* pattern */, ARMarkerInfo* * /* marker info */, int & /* number of detected markers */ );
int		AT_CheckDetected	( int /* marker pattern id */ );
int		AT_CheckDetected	( int /* marker pattern id */, ARMarkerInfo ** /* marker info */, int & /* number of detected markers */ );
bool	AT_GetTransMatrix	( AT_PATTERN & /* pattern */ );
bool	AT_GetTransMatrix	( AT_PATTERN & /* pattern */, double * /* transformation matrix */ );
bool	AT_GetTransMatrix	( AT_PATTERN & /* pattern */, ARMarkerInfo ** /* marker info */ );
bool	AT_GetTransMatrix	( AT_PATTERN & /* pattern */, ARMarkerInfo ** /* marker info */, double * /* transformation matrix */ );
void	AT_MultiMatrix		( AT_PATTERN & /* pattern */ );
void	AT_SwapBuffers		( void );
void	AT_Destory			( void );


/*	Do
________________________________________________________________*/

bool AT_Init( void )
{
	return AT_Init( 640, 480 );
}

bool AT_Init( int _w, int _h, const char *_campara )
{
	if( arParamLoad( _campara, 1, &AT_WP ) < 0 )
	{
		AT_ERROR = "camera parameter file error";
		return 1;
	}
	arParamChangeSize( &AT_WP, _w, _h, &AT_CP );
	arInitCparam( &AT_CP );
	arParamDisp( &AT_CP );
	argInit( &AT_CP, 1.0, 0, 0, 0, 1 );
	return 0;
}

int AT_SetThresh( int _t )
{
	if( _t > 0 )
		return ( AT_THRESH = _t );
	return -1;
}

int AT_Detect( ARUint8 *_map )
{
	return arDetectMarker( _map, AT_THRESH, &AT_MARKER_INFO, &AT_MARKER_NUM );
}

int AT_Detect( ARUint8 *_map, int _thresh )
{
	return arDetectMarker( _map, _thresh, &AT_MARKER_INFO, &AT_MARKER_NUM );
}

int AT_Detect( ARUint8 *_map, int _thresh, ARMarkerInfo **_info, int &_num )
{
	return arDetectMarker( _map, _thresh, _info, &_num );
}

int AT_Detect( ARUint8 *_map, ARMarkerInfo **_info, int &_num )
{
	return arDetectMarker( _map, AT_THRESH, _info, &_num );
}

int AT_Detect( ARUint8 *_map, AT_PATTERN &_patt )
{
	arDetectMarker( _map, AT_THRESH, &AT_MARKER_INFO, &AT_MARKER_NUM );
	return AT_CheckDetected( _patt.id );
}

int AT_Detect( ARUint8 *_map, AT_PATTERN &_patt, int _thresh )
{
	arDetectMarker( _map, _thresh, &AT_MARKER_INFO, &AT_MARKER_NUM );
	return AT_CheckDetected( _patt.id );
}

int AT_Detect( ARUint8 *_map, AT_PATTERN &_patt, int _thresh, ARMarkerInfo **_info, int &_num )
{
	arDetectMarker( _map, _thresh, _info, &_num );
	return AT_CheckDetected( _patt.id, _info, _num );
}

int AT_Detect( ARUint8 *_map, AT_PATTERN &_patt, ARMarkerInfo **_info, int &_num )
{
	arDetectMarker( _map, AT_THRESH, _info, &_num );
	return AT_CheckDetected( _patt.id, _info, _num );
}

int AT_CheckDetected( int _id )
{
	static int j, k;
	k = -1;
	for( j = 0; j < AT_MARKER_NUM; j++ )
	{
		if( _id == AT_MARKER_INFO[j].id )
		{
			if( k == -1 )
				k = j;
			else if( AT_MARKER_INFO[k].cf < AT_MARKER_INFO[j].cf )
				k = j;
		}
	}
	return k;
}

int AT_CheckDetected( int _id, ARMarkerInfo **_info, int &_num )
{
	static int j, k;
	k = -1;
	for( j = 0; j < _num; j++ )
	{
		if( _id == (*_info)[j].id )
		{
			if( k == -1 )
				k = j;
			else if( (*_info)[k].cf < (*_info)[j].cf )
				k = j;
		}
	}
	return k;
}

bool AT_GetTransMatrix( AT_PATTERN &_patt )
{	
	if( arGetTransMat( &AT_MARKER_INFO[_patt.id], _patt.center, _patt.marker_width, _patt.trans ) < 0 )
		return 1;
	argConvGlpara( _patt.trans, _patt.opengl_trans );
	return 0;
}

bool AT_GetTransMatrix( AT_PATTERN &_patt, double *_gltrans )
{	
	if( arGetTransMat( &AT_MARKER_INFO[_patt.id], _patt.center, _patt.marker_width, _patt.trans ) < 0 )
		return 1;
	argConvGlpara( _patt.trans, _gltrans );
	return 0;
}

bool AT_GetTransMatrix( AT_PATTERN &_patt, ARMarkerInfo **_info )
{	
	if( arGetTransMat( &(*_info)[_patt.id], _patt.center, _patt.marker_width, _patt.trans ) < 0 )
		return 1;
	argConvGlpara( _patt.trans, _patt.opengl_trans );
	return 0;
}

bool AT_GetTransMatrix( AT_PATTERN &_patt, ARMarkerInfo **_info, double *_gltrans )
{	
	if( arGetTransMat( &(*_info)[_patt.id], _patt.center, _patt.marker_width, _patt.trans ) < 0 )
		return 1;
	argConvGlpara( _patt.trans, _gltrans );
	return 0;
}

void AT_SwapBuffers( void )
{
	argSwapBuffers();
}

void AT_Destory( void )
{
	argCleanup();
}



#endif