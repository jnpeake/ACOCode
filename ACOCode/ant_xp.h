#pragma once

#ifndef _ANT_INC_
#define _ANT_INC_

#include "vector.h"
#include "platform.h"
#include "localsearch.h"
#include "tsp.h"
#include "list"


class Vector; // forward declaration

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
	LocalSearch *localSearch;

	int *rIndices;
	RanluxGenerator rlgen;


	// random number generator data

	Vector rSeed;
	Vector rC0;
	Vector rC1;
	Vector factor;
	Vector ones; // constant

	int nVertPadded;
	float *index; // array of indices (0..nVerts), for iRoulette

	int *tabu; // tabu list

	int fallback( float *weights, int *tabu, int currentIndex, TSP *tsp);
	int csRoulette(float *weights, int *tabu, int nWeights, nearestNeighbour *nnList, int numNN);



	void ConstructTour( void ); 

	void Init( AntSystem *as, int * seeds );
};

#endif
