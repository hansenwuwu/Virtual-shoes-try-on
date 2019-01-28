#ifndef _NON_LINEAR_H_
#define _NON_LINEAR_H_

#ifndef _LINEAR_2_H_
#error Neet to include <linear2.h>
#else

#ifdef __cplusplus
extern "C"{
#endif

typedef Dim		Dim_in, Dim_out;
typedef Vecf	Vecf_unknown, Vecf_value;
typedef Vecd	Vecd_unknown, Vecd_value;
typedef Mat2f	Mat2f_value;
typedef Mat2d	Mat2d_value;

//
// Leverberg-Marqardt
//

typedef void ( *LevMarFuncf )( Vecf_in, Vecf_out );
typedef void ( *LevMarFuncd )( Vecd_in, Vecd_out );
typedef void ( *LevMarFunc2f )( Vecf_in, Mat2f_out );
typedef void ( *LevMarFunc2d )( Vecd_in, Mat2d_out );

Errf LevMarf( LevMarFuncf, Vecf_unknown, Vecf_value, Dim_in, Dim_out, 
			  IterTime, Paramf, Errf, Errf, Errf );

Errd LevMard( LevMarFuncd, Vecd_unknown, Vecd_value, Dim_in, Dim_out, 
			  IterTime, Paramd, Errd, Errd, Errd );

Errf LevMar2f( LevMarFunc2f, Vecf_unknown, Mat2f_value, Dim_in, Row_out, Col_out, 
			  IterTime, Paramf, Errf, Errf, Errf );

Errd LevMar2d( LevMarFunc2d, Vecd_unknown, Mat2d_value, Dim_in, Row_out, Col_out, 
			  IterTime, Paramd, Errd, Errd, Errd );

#ifdef __cplusplus
}
#endif

#endif

#endif