#ifndef _ANT_INC_
#define _ANT_INC_

#include "platform.h"
#include "tsp.h"
#include "list"
#include "vector.h"


//#define USE_VROULETTE

class AntSystem; // forward declaration

class Ant
{
public:
	int randSeeds[16];
	float tourDist;
	float *rouletteVals;
	int *tour;
	int *remaining;
	int nTour;
	AntSystem *m_as;

	int *rIndices;

	// random number generator data

	unsigned int rSeed[16];
	Vector rC0;
	Vector rC1;
	Vector factor;
	Vector ones; // constant

	int nVert16;
	float *index; // array of indices (0..nVerts), for iRoulette

	int *tabu; // tabu list

	void seedAvxRandom( int *seeds ); 

	int iRoulette( float *weights, int *tabu, int nWeights);
	int csRoulette(float *weights, int *tabu, int nWeights, nearestNeighbour *nnList, int numNN);



	void ConstructTour( void ); 

	void Init( AntSystem *as, int * seeds );
};

#endif
