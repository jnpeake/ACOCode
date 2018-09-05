#include "ant_xp.h"
#include "antsystem_xp.h"
#include "tsp.h"
#include <cstdio> 
#include <cstdlib>
#include <iostream>

void Ant::Init( AntSystem *as, int *seeds )
{

	m_as = as;
	remaining = (int*)ALLOC( (m_as->GetTSP()->numVerts * sizeof(int)) );
	tour = (int*)ALLOC( (m_as->GetTSP()->numVerts * sizeof(int)) );
	nTour = 0;

	rIndices = (int*)ALLOC( m_as->GetTSP()->numVerts * sizeof(int) );

	// allocate the index and tabu arrays
	nVert16 = m_as->GetTSP()->numVerts;
	if ( nVert16%16 )
		nVert16 = 16 * ( nVert16/16 +1 ); // pad to multiple of 16 floats
	index = (float*)ALLOC( nVert16 * sizeof(float) ); // 64-byte aligned

	tabu = (int*)ALLOC( nVert16 * sizeof(int) ); // 64-byte aligned

	// fill the index array
	for ( int i = 0; i < nVert16; i++ )
		index[i] = (float)i;

	// intialize the RNG

	seedVecRandom(rC0,rC1,factor,seeds,rSeed);
	// handy constant
	float f1 = 1.0f;
	ones.set1(f1);
}

int Ant::iRoulette( float *weights, int *tabu, int nWeights )
{

	__declspec(align(64)) float indexSeed[16] = {0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 
												 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f};
	__declspec(align(64)) float indexStep[16] = { 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f,
												  16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f };
	__declspec(align(64)) float minusOnes[16] = { -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f,
												  -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f };

	Vector minusOne;
	minusOne.load(minusOnes);
	Vector runningIndex;
	runningIndex.load(indexSeed);
	Vector delta16;
	delta16.load( indexStep );
	Vector curIndices = minusOne;
	Vector curWeights = minusOne;
 	Vector tabuMask = int2mask( tabu[0] );
	for ( int i = 0; i < nWeights/16; i++ )
	{
		Vector nextWeights;
		nextWeights.load(weights + i * 16);
		Vector nextIndices =  runningIndex;
		Vector randoms = vecRandom(rC0,rC1,factor,rSeed);
		
		
		nextWeights = nextWeights * randoms;
		nextWeights = mask_mov( nextWeights, tabuMask, minusOne );

		tabuMask = int2mask( tabu[i+1] );
		maxLocStep( curWeights, curIndices, nextWeights, nextIndices );
		runningIndex = runningIndex + delta16 ;
	}
	// now reduce the elements of curWeights
	int reduced = reduceMax( curWeights, curIndices );
	return reduced;
}

int Ant::csRoulette(float *weights, int *tabu, int nVerts, nearestNeighbour *nnList, int numNN)
{
	__declspec(align(64)) float indexSeed[16] = { 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f,
		8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f};
	__declspec(align(64)) float minusOnes[16] = { -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f,
		-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f };

	Vector minusOne;
	minusOne.load(minusOnes);
	Vector baseIndex;
	baseIndex.load(indexSeed);
	Vector curIndices = minusOne;
	Vector curWeights = minusOne;

	for(int i = 0; i < numNN; i++)
	{	

		if(nnList[i].vectIndex != -1)
		{

			Vector randoms = vecRandom(rC0,rC1,factor,rSeed);
			Vector nextIndices;		
			Vector tabuMask = int2mask(tabu[nnList[i].vectIndex]);					
			Vector nnMask = int2mask(nnList[i].nnMask);
		
			Vector nextWeights;
			nextWeights.load(weights + nnList[i].vectIndex *16);
			//nextWeights.maskLoad(minusOne, nnMask, weights + nnList[i].vectIndex *16);
			nextWeights = nextWeights * randoms;
			nextWeights = mask_mov(minusOne, nnMask, nextWeights);
			nextWeights = mask_mov(nextWeights, tabuMask, minusOne);
		
			float offset = 16*nnList[i].vectIndex;
			nextIndices.set1(offset);	
			nextIndices =  nextIndices + baseIndex;
			maxLocStep(curWeights, curIndices, nextWeights, nextIndices);

		}
		else
		{

			break;
		}
	}
	
	// now reduce the elements of curWeights
	int reduced = reduceMax(curWeights, curIndices);
	if (reduced < 0) {
		reduced = -1;
	}
	
	return reduced;
}

#define USE_IMCI 0


void Ant::ConstructTour( void )
{
	//printf("construct tour");
	int i, j;
	tourDist = 0.0f;
	nTour = 0;
	TSP *tsp = m_as->GetTSP();
	int numNN = tsp->numNN;
	// zero the tabu list
	memset( tabu, 0, nVert16*sizeof(int));
	// pick a random start city
	Vector randoms = vecRandom(rC0,rC1,factor,rSeed);
	__declspec(align(64)) float r[16];
	store(r, randoms);
	int iLast = tsp->numVerts * r[0]; 
	if ( iLast == tsp->numVerts )
		iLast = 0;
	tour[0] = iLast;

	int iTabu = (iLast/16);
	int jTabu = iLast%16;
	tabu[iTabu] |= (1<<jTabu);
	if(nVert16 != tsp->numVerts)
	{
		for(int i = tsp->numVerts; i < nVert16; i++)
		{
			int jTabu = i%16;
			tabu[i/16] |= (1<<jTabu);
		}
	}
	for ( i = 1; i < tsp->numVerts; i++ )
	{
#pragma noinline

		tour[i] = csRoulette(m_as->m_weights[tour[i - 1]], tabu, nVert16, tsp->neighbourVectors[tour[i - 1]], numNN);
		
		if (tour[i] == -1)
		{
			tour[i] = iRoulette(m_as->m_weights[tour[i - 1]], tabu, nVert16);
		}
		
		
		int iTabu = (tour[i]/16);
		int jTabu = tour[i]%16;
		//printf("\n%d",iTabu);fflush(stdout);
		tabu[iTabu] |= (1<<jTabu);
		//printf("\n%d, %d, %d",i,tour[i],tour[i-1]);fflush(stdout);
		tourDist += tsp->edgeDist[tour[i]][tour[i-1]];
	}
	tourDist += tsp->edgeDist[tour[0]][tour[tsp->numVerts-1]];
}


