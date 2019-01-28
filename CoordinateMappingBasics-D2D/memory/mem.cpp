
#include "mem.h"
#include <stdlib.h>

//
// Global Variables
// 

static	char *gptr;

//
// Definitions
//

Memory::Memory( void ){
	_mem = NULL; 
	_m_size = _used = 0;
}

Memory::~Memory( void ){
	if( _m_size ) Clear();
}

void * Memory::Allocate( int x ){
	if( x > ( int )_m_size ) {
		gptr = ( char * )realloc( _mem, x );
		if( gptr != NULL ) {
			_mem = gptr;
			_m_size = _used = x;
			return gptr;
		}
		else return NULL;
	}
	_used = x;
	return _mem;
}

int Memory::Assign( int n, ... ){
	int i, mem = 0;
	va_list va;

	// count total memory
	va_start( va, n );
	for( i = 0; i < n; i++ ){
		va_arg( va, void * );			// pointer
		mem += va_arg( va, int );		// allocated size
	}
	if( mem > ( int )_m_size ) return 0;
	// allocate to pointer
	mem = 0;
	va_start( va, n );
	for( i = 0; i < n; i++ ){
		*( va_arg( va, void ** ) ) = Offset( mem );
		mem += va_arg( va, int );
	}
	va_end( va );
	// update used memory block
	_used = mem;
	return mem;
}

void * Memory::PushAssign( int x ){
	if( _used + x > _m_size ) return NULL;
	_used += x;
	return ( char * )_mem + _used - x;
}

void Memory::PopAssign( int x ){
	_used = ( _used - x < 0 ) ? 0 : _used - x;
}

void * Memory::Offset( int x ){
	gptr = ( char * )_mem;
	return gptr + x;
}

void Memory::Clear( void ){
	if( _m_size ){	
		free( _mem );
		_mem = NULL;
		_m_size = 0;
	}
}

int Memory::Size( void ){
	return _m_size;
}
