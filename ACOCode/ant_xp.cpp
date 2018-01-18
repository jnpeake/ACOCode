#include "ant_xp.h"
#include "antsystem_xp.h"
#include "tsp.h"


void Ant::Init( AntSystem *as, int *seeds )
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
	seedAvxRandom( seeds );

	// handy constant
#ifndef EMULATE
	float f1 = 1.0f;
	ones = _mm512_extload_ps( &f1, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, 0 );
#endif
}

void Ant::TestRandom( void )
{
#ifdef EMULATE
	float rVals[16];
	avxRandom( rVals );
#else
	__declspec(align(64)) float rVals[16]; // random numbers
	_mm512_store_ps( rVals, avxRandom() );
#endif
	for ( int i = 0; i < 16; i++ )
		printf("%f ",rVals[i] );
	printf("\n");

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

#ifndef EMULATE
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

inline __m512 ReduceMax( __m512 valvec, __m512 ivec )
{
	// return a vector with all elements equal to ivec[imax] where
	// valvec[imax] is largest element of valvec
	__m512 permVal;
	__m512 permInd;
	__mmask16 maxMask;
	// swap with neighbour 
	permVal = _mm512_swizzle_ps( valvec, _MM_SWIZ_REG_CDAB );
	permInd = _mm512_swizzle_ps( ivec, _MM_SWIZ_REG_CDAB );
	maxMask = _mm512_cmp_ps_mask( valvec, permVal, _MM_CMPINT_GT );
	valvec = _mm512_mask_mov_ps( permVal, maxMask, valvec );
	ivec = _mm512_mask_mov_ps( permInd, maxMask, ivec );
	// swap pairs
	permVal = _mm512_swizzle_ps( valvec, _MM_SWIZ_REG_BADC );
	permInd = _mm512_swizzle_ps( ivec, _MM_SWIZ_REG_BADC );
	maxMask = _mm512_cmp_ps_mask( valvec, permVal, _MM_CMPINT_GT );
	valvec = _mm512_mask_mov_ps( permVal, maxMask, valvec );
	ivec = _mm512_mask_mov_ps( permInd, maxMask, ivec );
	// swap lanes
	permVal = _mm512_permute4f128_ps( valvec, 0xB1 ); // 2, 3, 0, 1
	permInd = _mm512_permute4f128_ps( ivec, 0xB1 );
	maxMask = _mm512_cmp_ps_mask( valvec, permVal, _MM_CMPINT_GT );
	valvec = _mm512_mask_mov_ps( permVal, maxMask, valvec );
	ivec = _mm512_mask_mov_ps( permInd, maxMask, ivec );
	// swap pairs of lanes
	permVal = _mm512_permute4f128_ps( valvec, 0x4E ); // 1, 0, 3, 2
	permInd = _mm512_permute4f128_ps( ivec, 0x4E );
	maxMask = _mm512_cmp_ps_mask( valvec, permVal, _MM_CMPINT_GT );
	ivec = _mm512_mask_mov_ps( permInd, maxMask, ivec );
	// all elements of ivec now contain index of maximum
	return ivec;

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
#define VECTOR_REDUCTION
#ifdef VECTOR_REDUCTION
	__m512 reduced = ReduceMax( curWeights, curIndices );
	__declspec(align(64)) float finalIndices[16];
	_mm512_store_ps( finalIndices, reduced );
	return (int)finalIndices[0];
#else
	// serial reduction
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
#endif
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
	int iLast = tsp->numVerts * r[0]; 
	if ( iLast == tsp->numVerts )
		iLast = 0;
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

#else
#ifdef VANILLA
void Ant::ConstructTour( void )
{
	int i, j;
	tourDist = 0.0f;
	nTour = 0;
	TSP *tsp = m_as->GetTSP();
	// zero the tabu list
	memset( tabu, 0, nVert16*sizeof(int));

	// pick a random start city
	float r[16];
	avxRandom( r );
	int iLast = tsp->numVerts * r[0]; 
	if ( iLast == tsp->numVerts )
		iLast = 0;
	tour[0] = iLast;

	int iTabu = (iLast/16);
	int jTabu = iLast%16;
	tabu[iTabu] |= (1<<jTabu);
	for ( i = 1; i < tsp->numVerts; i++ )
	{
		// first look through the nearest neighbours
		float tot = 0.0f;
		for ( j = 0; j < tsp->numVerts; j++ )
		{
			if ( (j != tour[i-1]) && !(tabu[j/16]&(1<<(j%16))) && m_as->m_fNN[tour[i-1]][j] > 0.0f )
			{
				tot += m_as->m_weights[tour[i-1]][j];
		  	}
		}
		if ( tot > 0.0f )
		{
			avxRandom(r);
			float roulette = tot * r[0];
			int iSel = -1;
			int iLast = -1;
			tot = 0.0f;
			for ( j = 0; j < tsp->numVerts; j++ )
			{
				if ( (j != tour[i-1]) && !(tabu[j/16]&(1<<(j%16))) && m_as->m_fNN[tour[i-1]][j] > 0.0f )
				{
					tot += m_as->m_weights[tour[i-1]][j];
					iLast = j;
					if ( tot >= roulette )
					{
						iSel = j;
						break;
					}
				}
			}
			if ( iSel == -1 )
				iSel = iLast;
			tour[i] = iSel;
		}
		else
		{
			//pick nearest neighbour
			float near = 1e20f;
			int iSel = -1;
			for ( j = 0; j < tsp->numVerts; j++ )
			{
				if ( (j != tour[i-1]) && !(tabu[j/16]&(1<<(j%16))) )
				{
					float dist = tsp->edgeDist[tour[i-1]][j];
					if ( dist < near )
					{
						near = dist;
						iSel = j;
					}
				}
			}
			tour[i] = iSel;
		}
		  
		int iTabu = (tour[i]/16);
		int jTabu = tour[i]%16;
		tabu[iTabu] |= (1<<jTabu);
		tourDist += tsp->edgeDist[tour[i]][tour[i-1]];
	}	
	tourDist += tsp->edgeDist[tour[0]][tour[tsp->numVerts-1]];
}
#else
void Ant::ConstructTour( void )
{
	int i, j;
	tourDist = 0.0f;
	nTour = 0;
	TSP *tsp = m_as->GetTSP();
	// zero the tabu list
	memset( tabu, 0, nVert16*sizeof(int));

	// pick a random start city
	float r[16];
	avxRandom( r );
	int iLast = tsp->numVerts * r[0]; 
	if ( iLast == tsp->numVerts )
		iLast = 0;
	tour[0] = iLast;

	int iTabu = (iLast/16);
	int jTabu = iLast%16;
	tabu[iTabu] |= (1<<jTabu);
	for ( i = 1; i < tsp->numVerts; i++ )
	{
#ifdef USE_VROULETTE
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
	if (tourDist == 0)
	{
		printf("Distance 0");
	}
}
#endif
int Ant::iRoulette( float *weights, int *tabu, int nWeights )
{
	return 0;
}

int Ant::vRoulette( float *weights, int *tabu, int nWeights )
{
	int i, j;
	
	float runningIndex[16];
	float nextWeights[16];
	float nextIndices[16];
	float curIndices[16];
	float curWeights[16];
	float randoms[16];

	int tabuMask = tabu[0];

	for ( i = 0; i < 16; i++ )
	{
		curIndices[i] = curWeights[i] = 0.0f;
		runningIndex[i] = (float)i;
	}	

	for ( int i = 0; i < nWeights/16; i++ )
	{
		memcpy( nextWeights, weights+i*16, 16*sizeof(float) );
		memcpy( nextIndices, runningIndex, 16*sizeof(float) );
		avxRandom( randoms );
		for ( j = 0; j < 16; j++ )
		{
			if ( tabuMask&(1<<j) )
				nextWeights[j] = 0.0f;
			curWeights[j] += nextWeights[j];
			float roulette = curWeights[j] * randoms[j];
			bool ltMask = roulette < nextWeights[j];
			if ( ltMask )
			{
				curIndices[j] = runningIndex[j];
			}
			runningIndex[j] += 16.0f;
		} 
		tabuMask = tabu[i+1];
	}
	// now perform a reduction
	float finalWeights[16];
	float finalIndices[16];
	float rVals[16]; // random numbers
	memcpy( finalWeights, curWeights, 16*sizeof(float) );
	memcpy( finalIndices, curIndices, 16*sizeof(float) );
	avxRandom( rVals );

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
void Ant::seedAvxRandom( int *seeds )
{
	memcpy( rSeed, seeds, 16*sizeof(int) );
}

void Ant::avxRandom( float *r )
{
	for ( int i = 0; i < 16; i++ )
	{
		rSeed[i] = rSeed[i]*1664525L + 1013904223L;
		r[i] = (float)rSeed[i] * 2.328306437087974e-10;
	}
}
#endif
