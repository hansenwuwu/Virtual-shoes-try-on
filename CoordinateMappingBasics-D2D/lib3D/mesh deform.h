#pragma once

#include <stdio.h>
#include <vector>
#include <list>
#include "memory\mem.h"
#include "lib3D\obj.h"

using namespace std;

typedef unsigned int	size_t;

namespace MeshDeform{

// datatypes
typedef float *					Parameters;
typedef list< vector< int > >	GroupList;

	class BaseDeform: virtual public Memory{
	protected:
		size_t _nv;
	public:
		virtual ~BaseDeform( void ){};
		
		// Modifiers
		virtual void Clear( void ) = 0;

		// Function
		virtual int Training( vector< Mesh* >&, Parameters ) = 0;
		virtual int Training( vector< Mesh >&, vector< Parameters >& ) = 0;
		virtual int Deform( Mesh &, Parameters ) = 0;
		virtual float Residual( Mesh &, Parameters );		// error between Mesh(S) and Mesh'

		// File operation
		virtual bool Import( wchar_t [] );
		virtual bool Export( wchar_t [] );
		virtual bool Import( FILE * ) = 0;
		virtual bool Export( FILE * ) = 0;

		// friend functions
		friend float ShapeSSE( BaseDeform &, GroupList &, Mesh &, Parameters );
		friend float DrawError( BaseDeform &, Mesh &, Parameters, float max );
		friend float DrawShapeError( BaseDeform &, GroupList &, Mesh &, Parameters );
	};

	class LinearDeform : public BaseDeform{
	protected:
		float *_coef;					// coefficient;
		size_t _dim;

	public:
		LinearDeform( size_t vertex = 0, size_t dim = 0 );	// default constructor
		~LinearDeform( void );								// destructor
		virtual LinearDeform & operator= ( const LinearDeform & );	// copy constructor

		// Modifiers
		virtual void Clear( void );

		// Functions
		virtual int Training( vector< Mesh* >&, Parameters );
		virtual int Training( vector< Mesh >&, vector< Parameters >& );
		virtual int Deform( Mesh &, Parameters );

		// File operation
		virtual bool Import( wchar_t p[] ){ return BaseDeform::Import(p); }
		virtual bool Export( wchar_t p[] ){ return BaseDeform::Export(p); }
		virtual bool Import( FILE * );
		virtual bool Export( FILE * );
	};

	class ConstGeometryDeform: public  LinearDeform{
	protected:
		typedef struct _cgd_group{
			vector< int > vertex;		// verteces
			float *w;					// rotation angles
			float *coef;				// linear regression of translation
		} CGD_Group;

	protected:
		Point3f *_v;					// mesh
		list< CGD_Group > _group;

	protected:
		void CheckGroup( void );

	public:
		ConstGeometryDeform( size_t vertex = 0, size_t dim = 0 );		// default constructor
		~ConstGeometryDeform( void );									// destructor
		virtual ConstGeometryDeform & operator= ( const ConstGeometryDeform & );// copy constructor

		// Modifiers 
		virtual int SetMesh( Mesh & );
		virtual void PushGroup( vector< int >& );
		virtual int PopGroup( vector< int >& );
		virtual void Clear( void );
		virtual GroupList GetGroupList( void );

		// Functions
		virtual int Training( vector< Mesh* >&, Parameters );
		virtual int Deform( Mesh &, Parameters );

		// File operation
		virtual bool Import( wchar_t p[] ){ return BaseDeform::Import(p); }
		virtual bool Export( wchar_t p[] ){ return BaseDeform::Export(p); }
		virtual bool Import( FILE * );
		virtual bool Export( FILE * );
	};

	// functions

	float ShapeSSE( BaseDeform &, GroupList &, Mesh &, Parameters );
	float ShapeMSE( BaseDeform &, GroupList &, Mesh &, Parameters );

	// drawing
	float DrawError( BaseDeform &, Mesh &, Parameters, float max );
	float DrawShapeError( BaseDeform &, GroupList &, Mesh &, Parameters, float max );
};