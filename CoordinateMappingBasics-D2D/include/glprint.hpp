
#ifndef __GL_PRINT__
#define __GL_PRINT__

#include <stdio.h>
#include <string>
#include <windows.h>
#include <time.h>
#include <math.h>
#include <gl/glut.h>


#define			GP_MAX_CHAR			128

#ifndef			GP_PI
	#define		GP_PI				3.1415926535f
#endif

#define			GP_ENG				ANSI_CHARSET
#define			GP_CHS				GB2312_CHARSET
#define			GP_CHT				DEFAULT_CHARSET


/*	font
_______________________________________________________________*/

	// set font type
	void	glFont		( int size, int char_set, const char *font_family );

	// print english string
	void	glPrint		( const char *string );

	// print integer number
	void	glPrint		( int integer );

	// print chinese string
	void	glPrintCT	( const char *string );

/*	draw
_______________________________________________________________*/

	void	glCircle	( float radius, int edge_number );
	void	glEllipse	( float x_axis, float y_axis, int edge_number );

/*	fps
_______________________________________________________________*/

class	glFPS
{
	public:
		glFPS( void );
		~glFPS( void );
		int Get( void );

		int Size( void ) { return size; }
		bool Size( int font_size );

		int Lang( void ) { return lang; }
		void Lang( int font_lang );

		char *Font( void ) { return font; }
		void Font( const std::string string );

	protected:
		int fps, frame;
		clock_t lastTime, currentTime;

		int size;
		std::string font_s;
		char *font;
		int lang;
};


/*===============================================================
	Do
===============================================================*/

/*
/*	font
_______________________________________________________________*/

void glFont( int _size, int _char, const char *_font )
{
	HFONT font = CreateFontA( _size, 0, 0, 0, FW_MEDIUM, 0, 0, 0, _char, 
		OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, DEFAULT_QUALITY, DEFAULT_PITCH | FF_SWISS, _font );
	HFONT lastFont = (HFONT)SelectObject( wglGetCurrentDC(), font );
	DeleteObject( lastFont );
}

void glPrint( const char *_string )
{
	bool firstCall = true;
	static GLuint lists;

	if( firstCall )
	{
		firstCall = false;
		lists = glGenLists( GP_MAX_CHAR );
		wglUseFontBitmaps( wglGetCurrentDC(), 0, GP_MAX_CHAR, lists );
	}

	for( ; _string[0] != '\0'; _string++ )
		glCallList( lists + _string[0] );
}

void glPrint( int _integer )
{
	char c[1];
	sprintf_s( c, sizeof(int), "%d", _integer );
	glPrint( c );
}

void glPrintCT( const char *_string )
{
	int leng, i;
	const char *ptrc = _string;
	wchar_t *wstring, *ptrw;
	HDC hDC = wglGetCurrentDC();
	GLuint list = glGenLists( 1 );

	leng = 0;
	for( ; ptrc[0] != '\0'; ptrc++, leng++ )
		if( IsDBCSLeadByte( ptrc[0] ) )
			ptrc++;

	wstring = (wchar_t*)malloc( ( leng + 1 )*sizeof(wchar_t) );
	MultiByteToWideChar( CP_ACP, MB_PRECOMPOSED, _string, -1, wstring, leng );
	wstring[leng] = L'\0';

	ptrw = wstring;
	for( i = 0; i < leng; i++, ptrw++ )
	{
		wglUseFontBitmapsW( hDC, ptrw[0], 1, list );
		glCallList( list );
	}

	free( wstring );
	glDeleteLists( list, 1 );
}


/*
/*	draw
_______________________________________________________________*/

void glCircle( float _radius, int _edge_number )
{
	float temp, mult = 2*GP_PI/_edge_number;
	glBegin( GL_LINE_LOOP );
		for( int i = 0; i < _edge_number; i++ )
		{
			temp = i*mult;
			glVertex3f( _radius*cosf( temp ), _radius*sinf( temp ), 0.f );
		}
	glEnd();
}

void glEllipse( float _x_length, float _y_length, int _edge_number )
{
	float temp, mult = 2*GP_PI/_edge_number;
	glBegin( GL_LINE_LOOP );
	for( int i = 0; i < _edge_number; i++ )
	{
		temp = i*mult;
		glVertex3f( _x_length*cosf( temp ), _y_length*sinf( temp ), 0.f );
	}
	glEnd();
}


/*
/*	fps
_______________________________________________________________*/

#pragma region :: CLASS : FPS ::

glFPS::glFPS( void )
{
	size = 30;
	lang = GP_CHT;
	font_s = "Fixedsys";
	font = &font_s[0];

	fps = 0;
	frame = 0;
	lastTime = clock();
}

glFPS::~glFPS( void )
{
}

int glFPS::Get( void )
{
	frame++;
	currentTime = clock();
	if( currentTime - lastTime > CLOCKS_PER_SEC )
	{
		fps = frame;
		lastTime = currentTime;
		frame = 0;
	}
	return fps;
}

bool glFPS::Size( int _size )
{
	if( _size <= 0 )
		return false;
	size = _size;
	return true;
}

void glFPS::Lang( int _lang )
{
	lang = _lang;
}

void glFPS::Font( const std::string _string )
{
	font_s = _string;
}

#pragma endregion


#endif