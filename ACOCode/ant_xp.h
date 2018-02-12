#ifndef _ANT_INC_
#define _ANT_INC_

#include "platform.h"
#include "tsp.h"
#include "list"
#ifndef EMULATE
#include <immintrin.h>
#endif

//#define USE_VROULETTE

class AntSystem; // forward declaration

class Ant
{
public:
	float tourDist;
	float *rouletteVals;
	int *tour;
	int *remaining;
	int nTour;
	AntSystem *m_as;

	int *rIndices;

	// random number generator data
#ifndef EMULATE
	__m512i rSeed;
	__m512i rC0;
	__m512i rC1;
	__m512 ones; // constant
#else
	unsigned int rSeed[16];
#endif
	int nVert16;
	float *index; // array of indices (0..nVerts), for iRoulette

	int *tabu; // tabu list

	void seedAvxRandom( int *seeds ); 
#ifndef EMULATE
	inline __m512 avxRandom( void );
#else
	void avxRandom( float *r );
#endif

	int iRoulette( float *weights, int *tabu, int nWeights /*std::list<nearestNeighbour> nnList*/);
	int vRoulette( float *weights, int *tabu, int nWeights );
	int csRoulette(float *weights, int *tabu, int nWeights, nearestNeighbour *nnList, int numNN);

	void ConstructTour( void ); 

	void Init( AntSystem *as, int *seeds );
	void PrintTour( void );
	void TestRandom( void );
};

#endif
