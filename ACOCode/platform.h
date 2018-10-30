#ifndef _PLATFORM_INC_
#define _PLATFORM_INC_

//#define EMULATE
//#define USE_VROULETTE
//#define AVX2
#define SISD
#define BITS 32
#define _VECSIZE 8
#ifndef EMULATE
#define ALLOC( _x ) _mm_malloc( (_x), 64 )
#define FREE( _x ) _mm_free( _x )
#define USE_OMP
#else
#define ALLOC( _x ) malloc( (_x) )
#define FREE( _x ) free( _x )
#endif

#ifdef __GNUC__
#define ALIGN(a,b) b __attribute__ ((aligned (a)))
#else
#define ALIGN(a,b) __declspec(align(a)) b
#endif



#ifdef USE_OMP
#include <omp.h>
#endif

#endif
