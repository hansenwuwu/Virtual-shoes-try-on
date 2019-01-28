
#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

//
// Macro
//

#ifndef _NOT_VC6_
#define _NOT_VC6_	( WINVER > 0x0400 )
#endif

#ifndef _IS_VC6_
#define _IS_VC6_	!_NOT_VC6_
#endif

