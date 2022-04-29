#include "ant.h"
#include "antsystem.h"
#include "ranluxgen.h"
#include "timers.h"

#include "tsp.h"
#include <cstdio> 
#include <cstdlib>
#include <iostream>

void Ant::Init( AntSystem *as, int *seeds )
{

	m_as = as;
	TSP *tsp = m_as->GetTSP();
	remaining = (int*)ALLOC( (m_as->GetTSP()->numVerts * sizeof(int)) );
	tour = (int*)ALLOC( ((m_as->GetTSP()->numVerts+1) * sizeof(int)) ); //+1 for local search code
	nTour = 0;

	rIndices = (int*)ALLOC( m_as->GetTSP()->numVerts * sizeof(int) );
	rlgen.init(seeds[0] );
	localSearch = new LocalSearch(tsp, rlgen, m_as->GetTSP()->numNN);
	timers = new Timers(4);


	// allocate the index and tabu arrays
	nVertPadded = m_as->GetTSP()->numVerts;
	if (nVertPadded %_VECSIZE )
		nVertPadded = _VECSIZE * (nVertPadded /_VECSIZE +1 ); // pad to multiple of 16 floats
	//index = (float*)ALLOC( nVertPadded * sizeof(float) ); // 64-byte aligned

	tabu = (int*)ALLOC(nVertPadded * sizeof(int)); // 64-byte aligned

	/*// fill the index array
	for ( int i = 0; i < nVertPadded; i++ )
		index[i] = (float)i;*/

	// intialize the RNG

	seedVecRandom(rC0,rC1,factor,seeds,rSeed);
	// handy constant
	float f1 = 1.0f;
	ones.set1(f1);
}


int Ant::fallback( float *weights, int *tabu, int currentIndex, TSP *tsp)
{
	float* vertX = tsp->vertX;
	float* vertY = tsp->vertY;
	int currentX = vertX[currentIndex];
	int currentY = vertY[currentIndex];
	bool nextFound = false;
	int tryCount = 10;
	float shortestDist = 3.4e+38;
	int nextPoint;


	for(int i = 0; i < tsp-> numVerts; i++)
	{
		float dist = tsp->CalcEdgeDistSquared(vertX[i],vertY[i],currentX, currentY);
		if( dist < shortestDist && !((tabu[i/_VECSIZE] >> i % _VECSIZE) & 1))
		{
			shortestDist = dist;
			nextPoint = i;
		}
	}

	return nextPoint;
}





int Ant::csRoulette(float *weights, int *tabu, int nVerts, nearestNeighbour *nnList, int numNN)
{

#ifdef AVX2
	ALIGN(float indexSeed[_VECSIZE]) = { 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f};
	ALIGN(float indexStep[_VECSIZE]) = { 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f, 8.0f};
	ALIGN(float minusOnes[_VECSIZE]) = { -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f};
	ALIGN(float nextIndicesArray[_VECSIZE]) = {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};
	ALIGN(float ones[_VECSIZE]) = { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};

#elif defined AVX512
	ALIGN(float indexSeed[_VECSIZE]) = { 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f};
	ALIGN(float indexStep[_VECSIZE]) = { 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f};
	ALIGN(float minusOnes[_VECSIZE]) = { -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f ,-1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f};
	ALIGN(float nextIndicesArray[_VECSIZE]) = {0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};
#endif
	Vector minusOne;
	minusOne.load(minusOnes);
	Vector runningIndex;
	runningIndex.load(indexSeed);
	Vector delta16;
	delta16.load( indexStep );
	Vector curIndices = minusOne;
	Vector curWeights = minusOne;
	Vector one;
	one.load(ones);

	for(int i = 0; i < numNN; i++)
	{	

		if(nnList[i].vectIndex != -1)
		{

			Vector randoms = vecRandom(rC0,rC1,factor,rSeed);
			Vector nextIndices = runningIndex;	
			Vector tabuMask = int2mask(tabu[nnList[i].vectIndex]);					
			Vector nnMask = int2mask(nnList[i].nnMask);
		
			Vector nextWeights;
			//printf("\n i:%d, i*_VECSIZE: %d, numVerts: %d",i,i *_VECSIZE,nVerts);
			nextWeights.load(weights + (i *_VECSIZE));
			//nextWeights.maskLoad(minusOne, nnMask, weights + nnList[i].vectIndex *16);

			/*if(m_as->stagnated == false)
				nextWeights = nextWeights * randoms;*/

			if(m_as->stagnated == false)
				nextWeights = nextWeights * randoms;

			else
				nextWeights = (one-nextWeights) * randoms;
			nextWeights = mask_mov(minusOne, nnMask, nextWeights);
			nextWeights = mask_mov(nextWeights, tabuMask, minusOne);

			
			maxLocStep(curWeights, curIndices, nextWeights, nextIndices);
			runningIndex = runningIndex + delta16 ;

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

	timers->Clear();
	timers->StartTimer(0);
	int fallbackTotal = 0;
	int i, j;
	tourDist = 0.0f;
	nTour = 0;
	TSP *tsp = m_as->GetTSP();
	int numNN = tsp->numNN;

	// zero the tabu list
	memset( tabu, 0, nVertPadded*sizeof(int));

	// pick a random start city
	Vector randoms = vecRandom(rC0,rC1,factor,rSeed);
	ALIGN(float r[_VECSIZE]);
	store(r, randoms);
	int iLast = tsp->numVerts * r[0]; 
	if ( iLast == tsp->numVerts )
		iLast = 0;
	tour[0] = iLast;

	int iTabu = (iLast/_VECSIZE);
	int jTabu = iLast%_VECSIZE;
	tabu[iTabu] |= (1<<jTabu);
	if(nVertPadded != tsp->numVerts)
	{
		for(int i = tsp->numVerts; i < nVertPadded; i++)
		{
			int jTabu = i%_VECSIZE;
			tabu[i/_VECSIZE] |= (1<<jTabu);
		}
	}
	int fallbackCount = 0;

	for ( i = 1; i < tsp->numVerts; i++ )
	{
#pragma noinline


		tour[i] = csRoulette(m_as->m_weights[tour[i - 1]], tabu, nVertPadded, tsp->neighbourVectors[tour[i - 1]], numNN);
		
		//std::cout << tour[i] << "\n";
		if (tour[i] == -1)
		{
			timers->StartTimer(2);
			tour[i] = fallback(m_as->m_weights[tour[i - 1]], tabu, tour[i-1], tsp);
			tourDist += tsp->CalcEdgeDist(tour[i],tour[i-1]);
			fallbackTotal += tsp->CalcEdgeDist(tour[i],tour[i-1]);
			timers->StopTimer(2);
			fallbackCount++;
			
			//printf("\n%d: %d | DISTANCE: %f | DIFFERENCE: %f | FALLBACK %d",i,tour[i],tourDist,tsp->CalcEdgeDist(tour[i],tour[i-1]),++fallbackCount);
		}

		else
		{
			tourDist += tsp->edgeDist[tour[i-1]][tour[i]];
			tour[i] = tsp->nnList[tour[i - 1]][tour[i]];
			//printf("\n%d: %d | DISTANCE: %f",i,tour[i],tourDist);
			
		}
		

		int iTabu = (tour[i]/_VECSIZE);
		int jTabu = tour[i]%_VECSIZE;
		//printf("\n%d",iTabu);fflush(stdout);
		tabu[iTabu] |= (1<<jTabu);
		//printf("\n%d, %d, %d",i,tour[i],tour[i-1]);fflush(stdout);
		

		
	}
	tourDist += tsp->CalcEdgeDist(tour[0],tour[tsp->numVerts-1]);
	//timers->StopTimer(0);
	//timers->StartTimer(1);
	//printf("\nFINAL DIST: %f",tourDist);
	//printf("\nTOTAL FALLBACK WEIGHT: %d",fallbackTotal);
	tour[tsp->numVerts] = tour[0];
	
	timers->StartTimer(1);
	localSearch->ThreeOpt(tour);
	timers->StopTimer(1);
	int newDist = 0;

	for(int i = 1; i < tsp->numVerts; i++)
	{
		newDist+= tsp->CalcEdgeDist(tour[i],tour[i-1]);
	}

	newDist+=tsp->CalcEdgeDist(tour[0],tour[tsp->numVerts-1]);
	tourDist = newDist;
	timers->StopTimer(0);
	float fallbackPercentage = ((float)fallbackCount/(float)tsp->numVerts)*100;


	//printf("\n Thread Times | Overall %f | Fallback %f | Local Search %f | Weight Calculation %f | Fallback Percentage %f%% | Total Fallback Weight %d", timers->GetTimer(0),timers->GetTimer(2),timers->GetTimer(1),timers->GetTimer(3),fallbackPercentage,fallbackTotal);
	//timers->StopTimer(1);
	//printf("\nTour Time:%f | LS Time:%f",timers->GetTimer(0), timers->GetTimer(1));

}




