#include "antsystem.h"
#include "ranluxgen.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <cstdio> 
#include <iostream>

#ifndef EMULATE
#include <immintrin.h>
#endif

#ifdef _WIN32
#define fmax( _a, _b ) ( _a > _b ? _a : _b )
#define fmin( _a, _b ) ( _a < _b ? _a : _b )
#endif

#define _dbgp( _x ) {printf(_x); printf("\n"); fflush(stdout);}

void AntSystem::Init( int nAnts, TSP *tsp, int seed )
{
	
	m_nAnts = nAnts;
	m_pTSP = tsp;
	#ifdef mapFB
	weightMap.clear();
	#endif
	m_shortestDist = 1e20f;

	timers = new Timers(3);

	fallbackCount = 0;
	usingNNCount = 0;

	int i;

	// intiailise pheromone matrix 
	//ALLOC varies depending on whether EMULATE is defined or not - if yes it's the same as malloc, if not it's _mm_malloc

	//m_pher, m_weights, m_iDistSq and m_fNN are all pointer-to-pointer-to-float matrices of numVerts * 8
	m_pher = (float**)ALLOC( (m_pTSP->numVerts * sizeof( float* )) );
	m_weights = (float**)ALLOC( (m_pTSP->numVerts * sizeof( float* )) );
	m_iDistSq = (float**)ALLOC( (m_pTSP->numVerts * sizeof( float* )) );
	m_fNN = (float **)ALLOC( (m_pTSP->numVerts * sizeof( float* )) );

	// for Xeon Phi, each row of the array needs to be padded to a whole number * _VECSIZE float and 64-byte aligned.
	int nRowAlloc = m_pTSP->matrixSize;

	if ( nRowAlloc%_VECSIZE )
		nRowAlloc = _VECSIZE * (nRowAlloc/_VECSIZE + 1 );

	for ( i = 0; i < m_pTSP->numVerts; i++ )
	{
		//each entry in the matrix is an array of nRowAlloc * 8

	    m_pher[i] = (float*)ALLOC( (nRowAlloc * sizeof( float )) ); 
	    m_weights[i] = (float*)ALLOC(( nRowAlloc * sizeof( float ) )); 
	    m_iDistSq[i] = (float*)ALLOC(( nRowAlloc * sizeof( float ) )); 
	    m_fNN[i] = (float*)ALLOC(( nRowAlloc * sizeof( float ) )); 
		// zero the memory
		memset( m_pher[i], 0, nRowAlloc * sizeof( float ) );
		memset( m_weights[i], 0, nRowAlloc * sizeof( float ) );
		memset( m_iDistSq[i], 0, nRowAlloc * sizeof( float ) );
		memset( m_fNN[i], 0, nRowAlloc * sizeof( float ) );

		
	}

	// create ants
	//m_pAnts is an array of ants, m_nAnts * 144
	m_pAnts = (Ant*)ALLOC( (m_nAnts * sizeof( Ant )) );
	// create a ranlux generator to generate seeds for the
	// ants (which use a numerical recipes quick and dirty generator).
	
	rlgen.init( seed );
	localSearch = new LocalSearch(tsp, rlgen, tsp->numNN);
	//the ants are initialised in this loop
	for ( i = 0; i < m_nAnts; i++ )
	{
		int ALIGN(seeds[_VECSIZE]);
		for ( int j = 0; j < _VECSIZE; j++ )
		{
			seeds[j] = rlgen.irand( 0xFFFF );
		}
		m_pAnts[i].Init( this, seeds );
	}
	m_shortestTour = (int*)ALLOC( (m_pTSP->numVerts * sizeof(int)) );

	rho = 0.3f;

	Clear();
}



void AntSystem::Clear( void )
{
	timers->Clear();
	m_shortestDist = 1e20f;

	float aDist = 0.0f;
	int i, j;
	float val;
	
	val = 1.0f/(m_pTSP->nnDist*rho);
	pher0 = val;
	
	//printf("\nnnTour: %f Initial pheromone: %f\n",m_pTSP->nnDist,val);fflush(stdout);


	//from 0 to numVerts
	for ( i = 0; i < m_pTSP->numVerts; i++ )
	{
		//if the specified number of nearest neighbours is not 0

		nearestNeighbour* nn = m_pTSP->neighbourVectors[i];
	//from 0 to numverts
		for ( j = 0; j < m_pTSP->matrixSize; j++ )
		{
			if(m_pTSP->edgeDist[i][j] != 0)
			{
			//pheromone value is set
			val = 1.0f/(m_pTSP->nnDist*rho);
			SetPheromoneValue(i,j,val, false);

			//inverse square of edge distances
			m_iDistSq[i][j] = 1.0f/(m_pTSP->edgeDist[i][j]*m_pTSP->edgeDist[i][j]);
			//pheromone / edgeDist^2
			m_weights[i][j] = (GetPheromoneValue(i,j, false)/(m_pTSP->edgeDist[i][j]*m_pTSP->edgeDist[i][j]));
			//printf("\n%f",m_weights[i][j]);
			
			//printf("%d, %d: %f | %f | %f \n",i,j,m_weights[i][j], GetPheromoneValue(i,j, false), m_pTSP->edgeDist[i][j]);
			}

			//else
			//{
			//	break;
			//}

		}
	}


	// e^(-1.30102999566) / numVerts)
	float a = exp( log(0.05) / (float)m_pTSP->numVerts );

	int n = m_pTSP->numVerts;

	//if there are not 0 nearest neighbours
	if ( m_pTSP->numNN != 0 )
		n = m_pTSP->numNN;

	//mmasConst is (1-a) / (n+1) * a * 0.5);
	mmasConst = (1.0f - a ) / ( ((float)n + 1.0f) * a * 0.5f);

}

void AntSystem::Evaporate( void )
{
	int nVerts = m_pTSP->numVerts;
#ifdef USE_OMP
#pragma omp parallel for
#endif
	for ( int i = 0; i < nVerts; i++ )
	{
	  for ( int j = 0; j < m_pTSP -> matrixSize; j++ )
	    {
			SetPheromoneValue(i,j,GetPheromoneValue(i,j, false) * (1.0f - rho), false);
	    }
	}
}

void AntSystem::DoTours( void )
{
	int iNewShortest = -1;
	float iterationShortest = 1e20f;
	int iterationShortestAnt = m_nAnts;
#ifdef USE_OMP // run an ant per thread in the tour construction phase
#pragma omp parallel for 
#endif

	//each ant constructs a tour
	for ( int i = 0; i < m_nAnts; i++ )
	{
		m_pAnts[i].ConstructTour();
		//printf("tour constructed");
	}

	
	
	//checks the ants for the shortest tour distance
	for ( int i = 0; i < m_nAnts; i++ )
	{
		if(m_pAnts[i].tourDist < iterationShortest)
		{
			iterationShortestAnt = i;
		}

		if ( m_pAnts[i].tourDist < m_shortestDist )
		{
			m_shortestDist = m_pAnts[i].tourDist;
			iNewShortest = i;
			// check the tour...
		}		
	}

	m_pAnts[iterationShortestAnt].tour[m_pTSP->numVerts] = m_pAnts[iterationShortestAnt].tour[0];
	//localSearch->ThreeOpt (m_pAnts[iterationShortestAnt].tour);

	int newDist = 0;

	for(int i = 1; i < m_pTSP->numVerts; i++)
	{
		newDist+= m_pTSP->CalcEdgeDist(m_pAnts[iterationShortestAnt].tour[i], m_pAnts[iterationShortestAnt].tour[i-1]);
	}

	newDist+=m_pTSP->CalcEdgeDist(m_pAnts[iterationShortestAnt].tour[0], m_pAnts[iterationShortestAnt].tour[m_pTSP->numVerts-1]);

	if ( newDist < m_shortestDist )
		{
			m_shortestDist = newDist;
			iNewShortest = iterationShortestAnt;
		}	

	

	// copy new shortest tour, if any
	if ( iNewShortest != -1 )
	{
		//m_pAnts[iNewShortest].tour[m_pTSP->numVerts] = m_pAnts[iNewShortest].tour[0];
		//localSearch->TwoOpt (m_pAnts[iNewShortest].tour);
//		printf("new shortest ant %d\n",iNewShortest);
		memcpy( m_shortestTour, m_pAnts[iNewShortest].tour, m_pTSP->numVerts * sizeof(int) );
	}
}

void AntSystem::DepositFromTour( int *tour, float tourLength )
{
	float deltaPher = 1.0f / tourLength;
	
	//adds pheromone to the specified edges
	for ( int i = 0; i < m_pTSP->numVerts; i++ )
	{
		SetPheromoneValue(tour[i],tour[(i+1)%m_pTSP->numVerts],deltaPher, true);
		//m_pher[tour[i]][tour[(i+1)%m_pTSP->numVerts]] += deltaPher;
		//m_pher[tour[(i+1)%m_pTSP->numVerts]][tour[i]] += deltaPher;
	}

}
void AntSystem::Deposit( void )
{
	// deposit from the best ant only
	int iBest = -1;
	static int count = 0;

	float dBest = 1e20f;

	//determines the best ant
	for ( int i = 0; i < m_nAnts; i++ )
	{
		if ( m_pAnts[i].tourDist < dBest )
		{
			dBest = m_pAnts[i].tourDist;
			iBest = i;
		}
	}

	DepositFromTour( m_pAnts[iBest].tour, dBest );
	ALIGN(float pherMin);
	ALIGN(float pherMax);
	ALIGN(float evapFac) = 1.0f - rho;
	ALIGN(float one) = 1.0f;

	pherMax = 1.0f/(rho * m_shortestDist );
	pherMin = pherMax * mmasConst;

	Vector vEvap;
	vEvap.set1(evapFac);
	Vector vOnes;
	vOnes.set1(one);
	#ifdef mapFB
	MapEvap(evapFac, pherMax, pherMin);
	#endif


	// precompute the probability into, and enforce min and max pheromone levels for MMAS
#ifdef USE_OMP
#pragma omp parallel for
#endif
	for ( int i = 0; i < m_pTSP->numVerts; i++ )
	{
 		for ( int j = 0; j < m_pTSP->matrixSize; j+=_VECSIZE )
		{
			
			Vector pher;
			pher.load(m_pher[i] + j);
			

			Vector weights;
			weights.load(m_weights[i] + j);

			Vector iDist;
			iDist.load(m_iDistSq[i] + j);

			Vector nnFac;
			nnFac.load(m_fNN[i] + j);

			// evaporate
			pher = pher * vEvap;
			// MMAS - clamp pheromone between limits
			
			pher.vecMax(pherMax);
			pher.vecMin(pherMin);
				
			weights = pher * iDist;

			store(m_pher[i] + j, pher);
		    store(m_weights[i] + j, weights);
		}
	}
}
void AntSystem::CalcStagnationMetrics( void )
{
	meanEntropy = 0.0f;
	meanLambda = 0.0f;
	int i, j;

	for ( i = 0; i < m_pTSP->numVerts; i++ )
	{
		float entropy = 0.0f;
		float lambda = 0.0f;
		float pMin, pMax;
		pMin = 1e20f;
		pMax = 0.0f;
		for ( j = 0; j < m_pTSP->numVerts; j++ )
		{
			float p;
			p = m_pher[i][j] / (m_pTSP->edgeDist[i][j] * m_pTSP->edgeDist[i][j] );

			if ( m_pher[i][j] < pMin )
				pMin = m_pher[i][j];
			if ( m_pher[i][j] > pMax )
				pMax = m_pher[i][j];
			if ( p > 0.0f )
				entropy -= p*log(p);
		}
		float pThresh = pMin + 0.05f*(pMax - pMin);
		for ( j = 0; j < m_pTSP->numVerts; j++ )
		{
			float p = m_pher[i][j];
			if ( p > pThresh &&  i != j)
				lambda += 1.0f;
			if ( p > pThresh )
			{
				lambda += 1.0f;			
			}
		}
		meanEntropy += entropy;
	meanLambda += lambda;
	}
	meanEntropy /= (float)m_pTSP->numVerts;
	meanLambda /= (float)m_pTSP->numVerts;
}

void AntSystem::Iterate( void )
{
	timers->StartTimer(0);

  	timers->StartTimer(1);
  	DoTours(); 
  	timers->StopTimer(1);
  	timers->StartTimer(2);
  	Deposit();
  	timers->StopTimer(2);
	timers->StopTimer(0);

}

void AntSystem::Solve( int maxIterations, int maxStagnantIterations, bool continueStagnant )
{
	
	//performs an initial clear of the any system
	Clear();
	float shortestSoFar = 1e20f;
	int sinceChange = 0;
	int i;
	stagnated = false;
	timers->Clear();
	iStagnation = -1; // sentinel value indicates no stagnation
	


	//loop will continue until the max number of iterations are reached
	for ( i = 0; i < maxIterations && !(stagnated && !continueStagnant); i++ )
	{
		
		Iterate();
		//std::cout << i << "\n";
		if (i % 10 == 0)
		{
			#ifdef mapFB
			printf("\n%d %f %f %f %f %d", i, m_shortestDist, timers->GetTimer(1), timers->GetTimer(2), timers->GetTimer(0), weightMap.size());fflush(stdout);
			#else
			//printf("\n%d %f %f %f %f %d", i, m_shortestDist, timers->GetTimer(1), timers->GetTimer(2), timers->GetTimer(0));fflush(stdout);
			#endif

		}
		
		
		if ( m_shortestDist < shortestSoFar )
		{
			shortestSoFar = m_shortestDist;
			sinceChange = 0;
			stagnated = false;
		}
		else
			sinceChange++;
		if ( !stagnated && (sinceChange >= maxStagnantIterations) )
		{
			iStagnation = i;
			stagnated = true;
		}



		
			
	}
#ifdef EMULATE
	printf("\nFINAL: Iteration: %d, Shortest Distance: %f, Shortest Tour: %d, ", i, m_shortestDist, m_shortestTour, timers->GetTimer(0), timers->GetTimer(1));
	printf("\nTIMERS: Tour %f Pheromone %f\n",timers->GetTimer(0),timers->GetTimer(1));
#endif
//	CalcStagnationMetrics();
}
/*
bool compareFunc( TourSort &t1, TourSort &t2 )
{
	return t1.length < t2.length;
}
struct TourSort
{
	int index;
	float length;
};
*/
//std::vector<TourSort> stuff;
//std::sort( stuff.begin(), stuff.end(), compareFunc );


float AntSystem::GetPheromoneValue(int pointA, int pointB, bool afterTour)
{

	if(afterTour == false)
	{
		return m_pher[pointA][pointB];
	}
	#ifdef mapFB

	else
	{
		if(pointB > pointA)
		{
			int tempA = pointA;
			pointA = pointB;
			pointB = tempA;
		}

		uint64_t hash = pointA*m_pTSP->numVerts + pointB; 
		
		if(weightMap.find(hash) != weightMap.end())
		{
			return weightMap[hash];
		}

		else
		{
			return pher0;
		}
	}

	#endif
	
}

void AntSystem::SetPheromoneValue(int pointA, int pointB, float deltaPher, bool afterTour)
{
	bool aFound = false;
	bool bFound = false;
	bool nearestNeighbour;
	int aPosition;
	int bPosition;
	int** nnList = m_pTSP->nnList;
	int _numNN = m_pTSP->numNN;
	

	if(afterTour == true)
	{
		for(int i = 0; i < _numNN; i++)
		{
			if(nnList[pointA][i] == pointB)
			{
				bFound = true;
				bPosition = i;
			}

			if(nnList[pointB][i] == pointA)
			{
				aFound = true;
				aPosition = i;
			}

			if(aFound && bFound)
			{
				break;
			}

		}

		if(bFound)
		{
			m_pher[pointA][bPosition] += deltaPher;
		}

		if(aFound)
		{
			m_pher[pointB][aPosition] += deltaPher;
		}

		#ifdef mapFB
		if(!aFound && !bFound)
		{

			/*int vectorIndexA = pointA / _VECSIZE;
			int vectorIndexB = pointB / _VECSIZE;
			ALIGN(int indexStep[_VECSIZE]) = {0,1,2,3,4,5,6,7};
			Vector indexStepVec;
			Vector indexVec;

			indexVec.set1(vectorIndex);

			indexStepVec.load(indexStep);*/


			if(pointB > pointA)
			{
				int tempA = pointA;
				pointA = pointB;
				pointB = tempA;
			}


			uint64_t hash = pointA*m_pTSP->numVerts + pointB; 

			if(weightMap.find(hash) != weightMap.end())
			{
				weightMap[hash] += deltaPher;
			}

			else
			{
				weightMap[hash] = pher0 + deltaPher;
			}
		}
		#endif 
	}

	else
	{
		m_pher[pointA][pointB] += deltaPher;
	}
}



#ifdef mapFB
void AntSystem::MapEvap(float evapFac, float pherMax, float pherMin){
	
	std::map<uint64_t, float>::iterator it;
	for ( it = weightMap.begin(); it != weightMap.end(); it++ )
	{
		it->second = it->second * evapFac;
		if(it->second > pherMax)
		{
			it->second = pherMax;
		}

		else if (it->second < pherMin)
		{
			it->second = pherMin;
		}
	}

}
#endif
