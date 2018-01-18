#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>

#define ROULETTE_THREADS 1
__m512i rSeed[ROULETTE_THREADS];
	__m512i rC0;
	__m512i rC1;
	__m512 ones; // constant

void seedAvxRandom( int *seeds )
{
    unsigned int c0;
	unsigned int c1;
	// constants for Numerical Recipes quick and dirty RNG
	c0 = 1664525L;
	c1 = 1013904223L;
	
	for ( int i = 0; i < ROULETTE_THREADS; i++ )
	{
		rSeed[i] = _mm512_load_epi32( seeds + i*16 );
	}
	rC0 = _mm512_extload_epi32( &c0, _MM_UPCONV_EPI32_NONE, _MM_BROADCAST_1X16, 0 );
	rC1 = _mm512_extload_epi32( &c1, _MM_UPCONV_EPI32_NONE, _MM_BROADCAST_1X16, 0 );
}
__m512 avxRandom( int index )
{
	// use AVX fused multiply-add to iterate RNG
	rSeed[index] = _mm512_fmadd_epi32( rC0, rSeed[index], rC1 );
	// convert to float in range 0 to 1 and return
	return _mm512_cvtfxpnt_round_adjustepu32_ps( rSeed[index], _MM_FROUND_TO_NEAREST_INT, _MM_EXPADJ_32 );
}

int vRoulette( float *weights, int *tabu, int nWeights )
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
		__m512 randoms = avxRandom(0);
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
	_mm512_store_ps( rVals, avxRandom(0) );

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

#define ALLOC( _x ) _mm_malloc( (_x), 64 )
#define FREE( _x ) _mm_free( _x )

#define NWEIGHTS 512

float frand()
{
	return (float)rand()/(float)RAND_MAX;
}

int main( int argc, char *argv[] )
{
	float *weights = (float*)ALLOC(NWEIGHTS * sizeof(float));
	int *counts = (int*)ALLOC(NWEIGHTS * sizeof(int));
	int *tabu = (int*)ALLOC(NWEIGHTS * sizeof(int));
	int *seeds = (int*)ALLOC( 16 * sizeof(int));

	// fill with random values, and zero counts
	float tabuProb = 0.3f;
	float totWeight = 0.0f;
	for ( int i = 0; i < NWEIGHTS/16; i++ )
		tabu[i] = 0;
	for ( int i = 0; i < NWEIGHTS; i++ )
	{
		weights[i] = rand();
		if ( frand() < tabuProb )
		{
			int iTabu = i/16;
			int jTabu = i%16;
			tabu[iTabu] |= (1<<jTabu);
		}
		else
			totWeight += weights[i];
		counts[i] = 0;
	}

	// seed RNG
	for ( int i = 0; i < 16; i++ )
	{
		seeds[i] = rand();
	}
	seedAvxRandom( seeds );

	// do the monte carlo
	int nTrials = 1000000;
	for ( int i = 0; i < nTrials; i++ )
	{
		int selected = vRoulette( weights, tabu, NWEIGHTS );
		counts[selected]++;
	}

	// output results
	FILE *fp = fopen("roulette.dat", "w" );
	for ( int i = 0; i < NWEIGHTS; i++ )
	{
		int iTabu, jTabu;
		iTabu = i/16;
		jTabu = i%16;
		if ( tabu[iTabu]&(1<<jTabu) )
			weights[i] = 0.0f;
		fprintf(fp, "%f %f\n", weights[i]/totWeight, (float)counts[i]/(float)nTrials );
	}
	fclose(fp);

	FREE( tabu );
	FREE( counts );
	FREE( weights );

	return 0;
}
