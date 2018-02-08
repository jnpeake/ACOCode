#include "antsystem_xp.h"
#include "ranluxgen.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

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
	m_shortestDist = 1e20f;

	timers = new Timers(5);

	int i;

	// intiailise pheromone matrix 
	//ALLOC varies depending on whether EMULATE is defined or not - if yes it's the same as malloc, if not it's _mm_malloc
	printf("\n Size of float*: %d", sizeof(float*));

	//m_pher, m_weights, m_iDistSq and m_fNN are all pointer-to-pointer-to-float matrices of numVerts * 8
	m_pher = (float**)ALLOC( (m_pTSP->numVerts * sizeof( float* )) );
	m_weights = (float**)ALLOC( (m_pTSP->numVerts * sizeof( float* )) );
	m_iDistSq = (float**)ALLOC( (m_pTSP->numVerts * sizeof( float* )) );
	m_fNN = (float **)ALLOC( (m_pTSP->numVerts * sizeof( float* )) );

	// for Xeon Phi, each row of the array needs to be padded to a whole number * 16 float and 64-byte aligned.
	int nRowAlloc = m_pTSP->numVerts;

	if ( nRowAlloc%16 )
		nRowAlloc = 16 * (nRowAlloc/16 + 1 );

	printf("\n nRowAlloc: %d", nRowAlloc);
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
	RanluxGenerator rlgen;
	rlgen.init( seed );

	//the ants are initialised in this loop
	for ( i = 0; i < m_nAnts; i++ )
	{
		int seeds[16];
		for ( int j = 0; j < 16; j++ )
		{
			seeds[j] = rlgen.irand( 0xFFFF );
		}
		m_pAnts[i].Init( this, seeds );
	}
	m_shortestTour = (int*)ALLOC( (m_pTSP->numVerts * sizeof(int)) );

	rho = 0.02f;

	Clear();
}

void AntSystem::Clear( void )
{
	//timers are reset
	timers->Clear();
	
	//shortest dist is initially set as a very large number
	m_shortestDist = 1e20f;

	printf("\n Shortest Distance: %d", m_shortestDist);

	float aDist = 0.0f;
	int i, j;
	float val;

	int *tour = (int*)malloc(m_pTSP->numVerts * sizeof(int));

	//sets each step of the tour to i
	for ( i = 0; i < m_pTSP->numVerts; i++ )
	{
		tour[i] = i;
	}

	//calculates the nnTour distance
	for ( i = 0; i < m_pTSP->numVerts-1; i++ )
	{
		float nearD = 1e20f;
		int nearI;
		for ( j = i+1; j < m_pTSP->numVerts; j++ )
		{
			//if the distance between edge i and edge j is less than near D, nearD is set to the distance and nearI becomes j
			if ( m_pTSP->edgeDist[tour[i]][tour[j]] < nearD )
			{
				nearD = m_pTSP->edgeDist[tour[i]][tour[j]];
				nearI = j;
			}
		}
		// add distance and swap new index into tour
		aDist += nearD;

	    int swap = tour[i+1];
		tour[i+1] = tour[nearI];
		tour[nearI] = swap;
	}

	//the distance between the last two vertices is added to aDist, as well as the distance between the last and first vertex
	aDist += m_pTSP->edgeDist[tour[m_pTSP->numVerts-2]][tour[m_pTSP->numVerts-1]];
	aDist += m_pTSP->edgeDist[tour[0]][tour[m_pTSP->numVerts-1]];

	float sanityCheck = 0.0f;
	for ( i = 0; i < m_pTSP->numVerts; i++ )
	{
		int i0 = tour[i];

		int i1 = tour[(i+1)%m_pTSP->numVerts];
		//printf("\n SANITY CHECK i1: %d", i1);
		sanityCheck += m_pTSP->edgeDist[i0][i1];
	}
	//memcpy( m_shortestTour, tour, m_pTSP->numVerts * sizeof( int ) );
	//m_shortestDist = aDist;

	//tour is de-allocated
	free(tour);

	//initial pheremone level is determined by nnTour

	val = 1.0f/((float)aDist*rho);
	printf("\n nnTour: %f Initial pheromone: %f\n",aDist,val);


	//from 0 to numVerts
	for ( i = 0; i < m_pTSP->numVerts; i++ )
	{
		//if the specified number of nearest neighbours is not 0
		if ( m_pTSP->numNN != 0 )
		{
			//first numVerts values of m_fNN set to 0
			memset( m_fNN[i], 0, m_pTSP->numVerts*sizeof(float));
			for ( j = 0; j < m_pTSP->numNN; j++ )
			{
				//sets the first numNN nearest neighbours to 1
				//m_fNN[i][m_pTSP->nnList[i][j]] = 1.0f;
			}
		}

		//if no nearest neighbours specified
		else
		{
			//everything is set as nearest neighbour
			for ( j = 0; j < m_pTSP->numVerts; j++ )
				m_fNN[i][j] = 1.0f;
		}

		//from 0 to numverts
		for ( j = 0; j < m_pTSP->numVerts; j++ )
		{
			//pheromone value is set
			m_pher[i][j] = val;

			//inverse square of edge distances
			m_iDistSq[i][j] = 1.0f/(m_pTSP->edgeDist[i][j]*m_pTSP->edgeDist[i][j]);

			//pheromone / edgeDist^2
			m_weights[i][j] = m_pher[i][j]/(m_pTSP->edgeDist[i][j]*m_pTSP->edgeDist[i][j]);
#ifndef VANILLA

			//weight = (weight) * (1 + (1000*m_fNN)) - makes nearest neighbour edges 1000x more likely to be picked
			m_weights[i][j] *= ( 1.0f + 1000.0f);
#endif
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
	  for ( int j = 0; j < nVerts; j++ )
	    {
			m_pher[i][j] *= (1.0f - rho);
	    }
	}
}

void AntSystem::DoTours( void )
{
	int iNewShortest = -1;
#ifdef USE_OMP // run an ant per thread in the tour construction phase
#pragma omp parallel for 
#endif

	//each ant constructs a tour
	for ( int i = 0; i < m_nAnts; i++ )
	{
		m_pAnts[i].ConstructTour();
	}


	//checks the ants for the shortest tour distance
	for ( int i = 0; i < m_nAnts; i++ )
	{
		if ( m_pAnts[i].tourDist < m_shortestDist )
		{
			m_shortestDist = m_pAnts[i].tourDist;
			iNewShortest = i;
			// check the tour...
		}
	}

	// copy new shortest tour, if any
	if ( iNewShortest != -1 )
	{
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
		m_pher[tour[i]][tour[(i+1)%m_pTSP->numVerts]] += deltaPher;
		m_pher[tour[(i+1)%m_pTSP->numVerts]][tour[i]] += deltaPher;
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
//	if ( (++count)%25 == 0 )
//	{
//		count = 0;
		DepositFromTour( m_pAnts[iBest].tour, dBest );
//	}
//	else
//	{
//		DepositFromTour( m_shortestTour, m_shortestDist );
//	}
	// calculate the max and min pheromone values here

#ifndef EMULATE
	__declspec(align(64)) float pherMin;
	__declspec(align(64)) float pherMax;
	__declspec(align(64)) float evapFac = (1.0f - rho);
	__declspec(align(64)) float one = 1.0f;
	__declspec(align(64)) float nnMult = 999.0f; // nearest-neighbour boost factor minus 1

	pherMax = 1.0f/(rho * m_shortestDist );
	pherMin = pherMax * mmasConst;

	__m512 vPherMin = _mm512_extload_ps( &pherMin, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, 0 );
	__m512 vPherMax = _mm512_extload_ps( &pherMax, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, 0 );
	__m512 vEvap =  _mm512_extload_ps( &evapFac, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, 0 );
	__m512 vOnes = _mm512_extload_ps( &one, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, 0 );
	__m512 vNNMult = _mm512_extload_ps( &nnMult, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, 0 );


	// precompute the probability into, and enforce min and max pheromone levels for MMAS
#ifdef USE_OMP
#pragma omp parallel for
#endif
	for ( int i = 0; i < m_pTSP->numVerts; i++ )
	{

 		for ( int j = 0; j < m_pTSP->numVerts; j+=16 )
		{
			__m512 pher = _mm512_load_ps( m_pher[i]+j );
			__m512 weights = _mm512_load_ps( m_weights[i]+j );
			__m512 iDist = _mm512_load_ps( m_iDistSq[i]+j );
			__m512 nnFac = _mm512_load_ps( m_fNN[i]+j );

			// evaporate
			pher = _mm512_mul_ps( pher, vEvap );
			// MMAS - clamp pheromone between limits
			pher = _mm512_max_ps( pher, vPherMin );
			pher = _mm512_min_ps( pher, vPherMax );

			weights = _mm512_mul_ps( pher, iDist );
			// apply nearest neighbour correction to weights
			nnFac = _mm512_mul_ps( nnFac, vNNMult );
			nnFac = _mm512_add_ps( nnFac, vOnes );
			weights = _mm512_mul_ps( weights, nnFac );

			_mm512_store_ps( m_pher[i]+j, pher );
			_mm512_store_ps( m_weights[i]+j, weights );
		}
	}
#else
	float pherMin;
	float pherMax;
	float evapFac = (1.0f - rho);
	float one = 1.0f;
	float nnMult = 999.0f; // nearest-neighbour boost factor minus 1

	pherMax = 1.0f/(rho * m_shortestDist );
	pherMin = pherMax * mmasConst;

	float pher[16];
	float weights[16];
	float iDist[16];
	float nnFac[16];

#ifdef USE_OMP
#pragma omp parallel for
#endif

	//PHEROMONE EVAPORATION

	//loops through the vertices
	for ( int i = 0; i < m_pTSP->numVerts; i++ )
	{
		//loops through every 16th vertex
 		for ( int j = 0; j < m_pTSP->numVerts; j+=16 )
		{
			memcpy(pher, m_pher[i]+j, 16*sizeof(float) );
			memcpy(weights, m_weights[i]+j, 16*sizeof(float) );
			memcpy(iDist,  m_iDistSq[i]+j, 16*sizeof(float) );
			memcpy(nnFac,  m_fNN[i]+j, 16*sizeof(float) );

			for ( int k = 0; k < 16; k++ )
			{
				pher[k] = pher[k] * evapFac;
				pher[k] = fmax( pher[k], pherMin );
				pher[k] = fmin( pher[k], pherMax );
				
				weights[k] = pher[k] * iDist[k];
#ifndef VANILLA
				nnFac[k] = nnFac[k] * nnMult;
				nnFac[k] += 1.0f;
				weights[k] *= nnFac[k];
#endif
			}
			memcpy( m_pher[i]+j, pher, 16*sizeof(float) );
			memcpy( m_weights[i]+j, weights, 16*sizeof(float) );
		}
	}
#endif
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
  DoTours(); 
  timers->StopTimer(0);

  timers->StartTimer(1);
  Deposit();
  timers->StopTimer(1);
}

void AntSystem::Solve( int maxIterations, int maxStagnantIterations, bool continueStagnant )
{

	//performs an initial clear of the any system
	Clear();
	float shortestSoFar = 1e20f;
	int sinceChange = 0;
	
	int i;
	bool stagnated = false;
	timers->Clear();
	iStagnation = -1; // sentinel value indicates no stagnation

	//loop will continue until the max number of iterations are reached
	for ( i = 0; i < maxIterations && !(stagnated && !continueStagnant); i++ )
	{
		Iterate();
#ifdef EMULATE
		if ( i%20 == 0 )
		  printf("Iteration: %d, Shortest Distance: %f, Timers: %f %f \n",i,m_shortestDist,timers->GetTimer(0), timers->GetTimer(1));
#endif
		if ( m_shortestDist < shortestSoFar )
		{
			shortestSoFar = m_shortestDist;
			sinceChange = 0;
		}
		else
			sinceChange++;
		if ( !stagnated && (sinceChange >= maxStagnantIterations) )
		{
			iStagnation = i;
			CalcStagnationMetrics();
			stagnationEntropy = meanEntropy;
			stagnationLambda = meanLambda;
			stagnationTourLength = m_shortestDist;
			stagnated = true;
		}
	}
#ifdef EMULATE
	printf("Iteration: %d, Shortest Distance: %f, Shortest Tour: %d, Timers: %f %f \n", i, m_shortestDist, m_shortestTour, timers->GetTimer(0), timers->GetTimer(1));
	printf("Timers: tour %e Pheromone %e\n",timers->GetTimer(0),timers->GetTimer(1));
#endif
//	CalcStagnationMetrics();
}
