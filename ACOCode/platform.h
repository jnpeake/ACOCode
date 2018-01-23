#ifndef _PLATFORM_INC_
#define _PLATFORM_INC_

#define EMULATE
//#define USE_VROULETTE

#ifndef EMULATE
#define ALLOC( _x ) _mm_malloc( (_x), 64 )
#define FREE( _x ) _mm_free( _x )
#define USE_OMP
#else
#define ALLOC( _x ) malloc( (_x) )
#define FREE( _x ) free( _x )
#endif


#ifdef USE_OMP
#include <omp.h>
#endif

#endif
