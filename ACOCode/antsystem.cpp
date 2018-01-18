#include "antsystem.h"
#include "ranluxgen.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <immintrin.h>

void AntSystem::Init( int nAnts, TSP *tsp, int seed )
{
	m_nAnts = nAnts;
	m_pTSP = tsp;
	m_shortestDist = 1e20f;

	timers = new Timers(5);

	int i;

	// intiailise pheromone matrix 

	m_pher = (float**)ALLOC( (m_pTSP->numVerts * sizeof( float* )) );
	m_weights = (float**)ALLOC( (m_pTSP->numVerts * sizeof( float* )) );
	m_iDistSq = (float**)ALLOC( (m_pTSP->numVerts * sizeof( float* )) );
	m_fNN = (float **)ALLOC( (m_pTSP->numVerts * sizeof( float* )) );
	// for Xeon Phi, each row of the array needs to be padded to a whole number * 16 float and 64-byte aligned.
	int nRowAlloc = m_pTSP->numVerts;

	if ( nRowAlloc%16 )
		nRowAlloc = 16 * (nRowAlloc/16 + 1 );


	for ( i = 0; i < m_pTSP->numVerts; i++ )
	{
	    m_pher[i] = (float*)ALLOC( (nRowAlloc * sizeof( float )) ); 
	    m_weights[i] = (float*)ALLOC(( nRowAlloc * sizeof( float ) )); 
	    m_iDistSq[i] = (float*)ALLOC(( nRowAlloc * sizeof( float ) )); 
	    m_fNN[i] = (float*)ALLOC(( nRowAlloc * sizeof( float ) )); 
	}
	// create ants
	m_pAnts = (Ant*)ALLOC( (m_nAnts * sizeof( Ant )) );
	// create a ranlux generator to generate seeds for the
	// ants (which use a numerical recipes quick and dirty generator).
	RanluxGenerator rlgen;
	rlgen.init( seed );
	for ( i = 0; i < m_nAnts; i++ )
	{
		m_pAnts[i].Init( this, rlgen.irand(0xFFFF) );
	}
	m_shortestTour = (int*)ALLOC( (m_pTSP->numVerts * sizeof(int)) );

	rho = 0.98f;

	Clear();
}

void AntSystem::Clear( void )
{
	m_shortestDist = 1e20f;

	float aDist = 0.0f;
	int i, j;
	float val;

	int *tour = (int*)malloc(m_pTSP->numVerts * sizeof(int));
	for ( i = 0; i < m_pTSP->numVerts; i++ )
	{
		tour[i] = i;
	}

	for ( i = 0; i < m_pTSP->numVerts-1; i++ )
	{
		float nearD = 1e20f;
		int nearI;
		for ( j = i+1; j < m_pTSP->numVerts; j++ )
		{
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
	aDist += m_pTSP->edgeDist[tour[m_pTSP->numVerts-2]][tour[m_pTSP->numVerts-1]];
	aDist += m_pTSP->edgeDist[tour[0]][tour[m_pTSP->numVerts-1]];
	float sanityCheck = 0.0f;
	for ( i = 0; i < m_pTSP->numVerts; i++ )
	{
		int i0 = tour[i];
		int i1 = tour[(i+1)%m_pTSP->numVerts];
		sanityCheck += m_pTSP->edgeDist[i0][i1];
	}
	//memcpy( m_shortestTour, tour, m_pTSP->numVerts * sizeof( int ) );
	//m_shortestDist = aDist;
	free(tour);

	val = 1.0f/((float)aDist*rho);

	for ( i = 0; i < m_pTSP->numVerts; i++ )
	{
		if ( m_pTSP->numNN != 0 )
		{
			memset( m_fNN[i], 0, m_pTSP->numVerts*sizeof(float));
			for ( j = 0; j < m_pTSP->numNN; j++ )
			{
				m_fNN[i][m_pTSP->nnList[i][j]] = 1.0f;
			}
		}
		else
		{
			for ( j = 0; j < m_pTSP->numVerts; j++ )
				m_fNN[i][j] = 1.0f;
		}

		for ( j = 0; j < m_pTSP->numVerts; j++ )
		{
			m_pher[i][j] = val;
			m_iDistSq[i][j] = 1.0f/(m_pTSP->edgeDist[i][j]*m_pTSP->edgeDist[i][j]);
			m_weights[i][j] = m_pher[i][j]/(m_pTSP->edgeDist[i][j]*m_pTSP->edgeDist[i][j]);
		}
	}

	float a = exp( -2.995732274 / (float)m_pTSP->numVerts );
	mmasConst = (1.0f - a ) / (( ((float)m_pTSP->numVerts  * 0.5f)-1.0f) * a );

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
	      m_pher[i][j] *= rho;
	    }
	}
}

void AntSystem::DoTours( void )
{
	int iNewShortest = -1;

#ifdef USE_OMP // run an ant per thread in the tour construction phase
#pragma omp parallel for 
#endif
	for ( int i = 0; i < m_nAnts; i++ )
	{
		m_pAnts[i].ConstructTour();
	}

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
		memcpy( m_shortestTour, m_pAnts[iNewShortest].tour, m_pTSP->numVerts * sizeof(int) );
	}
}

void AntSystem::DepositFromTour( int *tour, float tourLength )
{
	float deltaPher = 1.0f / tourLength;
	
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
//		DepositFromTour( m_pAnts[iBest].tour, dBest );
//	}
//	else
//	{
		DepositFromTour( m_shortestTour, m_shortestDist );
//	}
	// calculate the max and min pheromone values here

	__declspec(align(64)) float pherMin;
	__declspec(align(64)) float pherMax;
	__declspec(align(64)) float evapFac = rho;
	__declspec(align(64)) float one = 1.0f;
	__declspec(align(64)) float nnMult = 1000.0f; // nearest-neighbour boost factor

	pherMax = 1.0f/(rho * m_shortestDist );
	pherMin = pherMax / (2 * m_pTSP->numVerts);//* mmasConst;

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
	Clear();
	float shortestSoFar = 1e20f;
	int sinceChange = 0;

	int i;
	bool stagnated = false;
	iStagnation = -1; // sentinel value indicates no stagnation
	timers->Clear();
	for ( i = 0; i < maxIterations && !(stagnated && !continueStagnant); i++ )
	{
		Iterate();
//		if ( i%20 == 0 )
//		  printf("Iteration %d shortest %e %e %e\n",i,m_shortestDist,timers->GetTimer(0), timers->GetTimer(1));
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

	//printf("Iteration %d shortest %e %e %e\n",i,m_shortestDist,timers->GetTimer(0), timers->GetTimer(1));
	//printf("Timers: tour %e Pheromone %e\n",timers->GetTimer(0),timers->GetTimer(1));
	CalcStagnationMetrics();
}
