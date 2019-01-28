
#ifndef __MATHBOX_NUMBER__
#define __MATHBOX_NUMBER__


#define		nRound( X )				(int)( (X) - .5 )
#define		nMax( A, B )			( (A) > (B) ) ? (A) : (B)
#define		nMin( A, B )			( (A) < (B) ) ? (A) : (B)
#define		n3Max( A, B, C )		( (A) > (B) ) ? ( ( (A) > (C) ) ? (A) : (C) ) : ( ( (B) > (C) ) ? (B) : (C) )
#define		n3Min( A, B, C )		( (A) < (B) ) ? ( ( (A) < (C) ) ? (A) : (C) ) : ( ( (B) < (C) ) ? (B) : (C) )
#define		nLUcmp( X, L, U )		( (X) < (L) ) ? (L) : ( ( (X) < (U) ) ? (X) : (U) )


#endif