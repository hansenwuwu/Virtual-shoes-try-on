#ifndef _RICE_MEMORY_POOL_H_
#define _RICE_MEMORY_POOL_H_

#include <stdarg.h>

//
// Memory pool: a stack structure
//

class Memory{
protected:
	void *_mem;
	unsigned int _used;
	unsigned int _m_size;

public:
	Memory( void );
	~Memory( void );

	virtual void * Allocate( int );
	virtual int Assign( int, ... );
	virtual void * PushAssign( int );
	virtual void PopAssign( int );
	virtual void * Offset( int );
	virtual void Clear( void );
	virtual int Size( void );
};

#endif