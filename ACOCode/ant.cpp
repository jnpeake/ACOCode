#include "ant.h"
#include "antsystem.h"
#include "tsp.h"


void Ant::Init( AntSystem *as, int seed )
{
	m_as = as;
	remaining = (int*)ALLOC( (m_as->GetTSP()->numVerts * sizeof(int)) );
	tour = (int*)ALLOC( (m_as->GetTSP()->numVerts * sizeof(int)) );
	rouletteVals = (float*)ALLOC( m_as->GetTSP()->numVerts * sizeof(float) );
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
	int *seeds = (int *)ALLOC( 16 * sizeof(int));
	srand(seed); 
	for ( int i = 0; i < 16; i++ )
	{
		seeds[i] = rand();
	}
	seedAvxRandom( seeds );
	FREE( seeds );

	// handy constant
	float f1 = 1.0f;
	ones = _mm512_extload_ps( &f1, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, 0 );

}

void Ant::PrintTour( void )
{
	int i;
	if ( m_as->GetTSP()->numVerts < 7 )
	{
		printf("[ ");
		for ( i = 0; i < m_as->GetTSP()->numVerts; i++ )
			printf("%d ",tour[i] );
		printf("]\n");
	}
	else
	{
		printf("[ %d %d %d ... %d %d %d ]\n",tour[0],tour[1],tour[2],
			   tour[m_as->GetTSP()->numVerts-3],tour[m_as->GetTSP()->numVerts-2],tour[m_as->GetTSP()->numVerts-1]);
	}
}

inline void maxLocStep( __m512 &oldWeights, __m512 &oldIndices, __m512 &newWeights, __m512 &newIndices )
{
	__mmask16 maxMask = _mm512_cmp_ps_mask( newWeights, oldWeights, _MM_CMPINT_GT );
	oldWeights = _mm512_mask_mov_ps( oldWeights, maxMask, newWeights );
	oldIndices = _mm512_mask_mov_ps( oldIndices, maxMask, newIndices );
}

int Ant::vRoulette( float *weights, int *tabu, int nWeights )
{
	__declspec(align(64)) float indexSeed[16] = {0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 
												 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f};
	__declspec(align(64)) float indexStep[16] = { 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f,
												  16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f };
	__declspec(align(64)) float zeroes[16] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
												  0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
 
	__m512 zero = _mm512_load_ps( zeroes );
	__m512 runningIndex = _mm512_load_ps( indexSeed );
	__m512 delta16 = _mm512_load_ps( indexStep );
	__m512 curIndices = zero;
	__m512 curWeights = zero;
 	__mmask16 tabuMask = _mm512_int2mask( tabu[0] );

	for ( int i = 0; i < nWeights/16; i++ )
	{
		__m512 nextWeights = _mm512_load_ps( weights + i*16 );
		__m512 nextIndices =  runningIndex;
		__m512 randoms = avxRandom();
		nextWeights = _mm512_mask_mov_ps( nextWeights, tabuMask, zero );
		curWeights = _mm512_add_ps( curWeights, nextWeights );
		__m512 roulette = _mm512_mul_ps( curWeights, randoms );
		__mmask16 ltMask = _mm512_cmp_ps_mask( roulette, nextWeights, _MM_CMPINT_LT );
		curIndices = _mm512_mask_mov_ps( curIndices, ltMask, runningIndex );
		tabuMask = _mm512_int2mask( tabu[i+1] );
		runningIndex = _mm512_add_ps( runningIndex, delta16 );
	}
	// now perform a reduction
	__declspec(align(64)) float finalWeights[16];
	__declspec(align(64)) float finalIndices[16];
	__declspec(align(64)) float rVals[16]; // random numbers
	_mm512_store_ps( finalWeights, curWeights );
	_mm512_store_ps( finalIndices, curIndices );
	_mm512_store_ps( rVals, avxRandom() );

	int nt = 8;
	int iRand = 0;
	for ( int i = 0; i < 4; i++ )
	{
		for ( int j = 0; j < nt; j++ )
		{
			finalWeights[j*2] += finalWeights[j*2+1];
			if ( finalWeights[j*2]*rVals[iRand++] < finalWeights[j*2+1] )
				finalIndices[j*2] = finalIndices[j*2+1];
			finalWeights[j] = finalWeights[j*2];
			finalIndices[j] = finalIndices[j*2];
		} 
		nt >>= 1;
	}
	return (int)finalIndices[0];
}

int Ant::iRoulette( float *weights, int *tabu, int nWeights )
{
	__declspec(align(64)) float indexSeed[16] = {0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 
												 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f};
	__declspec(align(64)) float indexStep[16] = { 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f,
												  16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f, 16.0f };
	__declspec(align(64)) float minusOnes[16] = { -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f,
												  -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f, -1.0f };

	__m512 minusOne = _mm512_load_ps( minusOnes );
	__m512 runningIndex = _mm512_load_ps( indexSeed );
	__m512 delta16 = _mm512_load_ps( indexStep );
	__m512 curIndices = minusOne;
	__m512 curWeights = minusOne;
 	__mmask16 tabuMask = _mm512_int2mask( tabu[0] );

	for ( int i = 0; i < nWeights/16; i++ )
	{
		__m512 nextWeights = _mm512_load_ps( weights + i*16);
		__m512 nextIndices =  runningIndex;
		__m512 randoms = avxRandom();
		nextWeights = _mm512_mul_ps( nextWeights, randoms );
		nextWeights = _mm512_mask_mov_ps( nextWeights, tabuMask, minusOne );
		tabuMask = _mm512_int2mask( tabu[i+1] );
		maxLocStep( curWeights, curIndices, nextWeights, nextIndices );
		runningIndex = _mm512_add_ps( runningIndex, delta16 );
	}

	// now reduce the elements of curWeights
	// maybe could do this with shuffle instructions to construct a parallel reduction ( needs 5 steps vs. 16 )
	__declspec(align(64)) float finalWeights[16];
	__declspec(align(64)) float finalIndices[16];
	_mm512_store_ps( finalWeights, curWeights );
	_mm512_store_ps( finalIndices, curIndices );
	float maxWeight = 0.0f;
	int indexMax = 0;
	for ( int i = 0; i < 16; i++ )
	{
		if ( finalWeights[i] > maxWeight )
		{
			maxWeight = finalWeights[i];
			indexMax = (int)finalIndices[i];
		}
	}
	return indexMax;
}

void Ant::seedAvxRandom( int *seeds )
{
    unsigned int c0;
	unsigned int c1;
	// constants for Numerical Recipes quick and dirty RNG
	c0 = 1664525L;
	c1 = 1013904223L;
	
	rSeed = _mm512_load_epi32( seeds );
	rC0 = _mm512_extload_epi32( &c0, _MM_UPCONV_EPI32_NONE, _MM_BROADCAST_1X16, 0 );
	rC1 = _mm512_extload_epi32( &c1, _MM_UPCONV_EPI32_NONE, _MM_BROADCAST_1X16, 0 );
}

inline __m512 Ant::avxRandom( void )
{
	// use AVX fused multiply-add to iterate RNG
	rSeed = _mm512_fmadd_epi32( rC0, rSeed, rC1 );
	// convert to float in range 0 to 1 and return
	return _mm512_cvtfxpnt_round_adjustepu32_ps( rSeed, _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_32 );
}


void Ant::ConstructTour( void )
{
	int i, j;
	tourDist = 0.0f;
	nTour = 0;
	TSP *tsp = m_as->GetTSP();
	// zero the tabu list
	memset( tabu, 0, nVert16*sizeof(int));

	// pick a random start city
	__m512 randoms = avxRandom();
	__declspec(align(64)) float r[16];
	_mm512_store_ps( r, randoms );
	int iLast = m_as->GetTSP()->numVerts * r[0]; 
	tour[0] = iLast;

	int iTabu = (iLast/16);
	int jTabu = iLast%16;
	tabu[iTabu] |= (1<<jTabu);

	for ( i = 1; i < tsp->numVerts; i++ )
	{
#pragma noinline
#if USE_VROULETTE
		tour[i] = vRoulette( m_as->m_weights[tour[i-1]], tabu, nVert16 );
#else
		tour[i] = iRoulette( m_as->m_weights[tour[i-1]], tabu, nVert16 );
#endif
		int iTabu = (tour[i]/16);
		int jTabu = tour[i]%16;
		tabu[iTabu] |= (1<<jTabu);
		tourDist += tsp->edgeDist[tour[i]][tour[i-1]];
	}	
	tourDist += tsp->edgeDist[tour[0]][tour[tsp->numVerts-1]];

}


