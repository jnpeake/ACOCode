#include "vector.h"
#include "antsystemhelp.h"
#include <iostream>
#define SISD
//#define AVX512
Vector int2mask(int maskInt)
{
#ifdef SISD

	Vector mask;
	for (int i = 0; i < 16; ++i) {
		mask.values[i] = (maskInt >> i) & 1;
	}


	return mask;

#elif defined AVX512
	Vector resultVector;
	resultVector.maskVec = _mm512_int2mask(maskInt);
	return resultVector;
#endif // AVX512

}

Vector mask_mov(const Vector& v1, Vector bitMask, const Vector& v2)
{
	Vector resultVector;
#ifdef SISD
	for (int i = 0; i < 16; i++)
	{
		if (bitMask.values[i] == 0)
		{
			resultVector.values[i] = v1.values[i];
		}

		else
		{
			resultVector.values[i] = v2.values[i];
		}
	}

	return resultVector;

#elif defined AVX512
	resultVector.AVXVec = _mm512_mask_mov_ps(v1.AVXVec, bitMask.maskVec, v2.AVXVec);
	return resultVector;
#endif // AVX512
}

Vector gtMask(const Vector& v1, const Vector& v2)
{
Vector resultMask;
#ifdef SISD

	for (int i = 0; i < 16; i++)
	{
		if (v1.values[i] > v2.values[i])
		{
			resultMask.values[i] = 1;
		}

		else
		{
			resultMask.values[i] = 0;
		}
	}

	return resultMask;

#elif defined AVX512
	resultMask.maskVec = _mm512_cmp_ps_mask(v1.AVXVec, v2.AVXVec, _MM_CMPINT_GT);
	return resultMask;
#endif 
}

Vector ltMask(const Vector& v1, const Vector& v2)
{
Vector resultMask;
#ifdef SISD
	for (int i = 0; i < 16; i++)
	{
		if (v1.values[i] < v2.values[i])
		{
			resultMask.values[i] = 1;
		}

		else
		{
			resultMask.values[i] = 0;
		}
	}
	return resultMask;
#elif defined AVX512
	resultMask.maskVec = _mm512_cmp_ps_mask(v1.AVXVec,v2.AVXVec,_MM_CMPINT_LT);
	return resultMask;
#endif // AVX512

}

Vector vecRandom(Vector rC0, Vector rC1, Vector factor, Vector rSeed)
{
	Vector r;
#ifdef SISD
	for (int i = 0; i < 16; i++)
	{
		
		rSeed.iValues[i] = rSeed.iValues[i] * 1664525L + 1013904223L;
		r.values[i] = (float)rSeed.iValues[i] * 2.328306437087974e-10;
	}

	return r;

#elif defined AVX512
	rSeed.AVXIntVec = _mm512_mullo_epi32(rC0.AVXIntVec, rSeed.AVXIntVec);
	rSeed.AVXIntVec = _mm512_add_epi32(rC1.AVXIntVec, rSeed.AVXIntVec);
	// convert to float in range 0 to 1 and return
	__m512 returnValue = _mm512_cvt_roundepu32_ps(rSeed.AVXIntVec, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
	r.AVXVec = _mm512_mul_ps(returnValue, factor.AVXVec);
	return r;
#endif // AVX512

	
}

void seedVecRandom(Vector& rC0, Vector& rC1, Vector& factor, int *seeds, Vector& rSeed)
{
#ifdef SISD
	for(int i = 0; i < 16; i++)
	{
		rSeed.iValues[i] = seeds[i];
	}


#elif defined AVX512
		__declspec(align(64)) unsigned int c0[16] = { 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L,1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L };
	__declspec(align(64)) unsigned int c1[16] = { 1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L, 	1013904223L, 1013904223L, 1013904223L,1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L };
	__declspec(align(64)) float factors[16] = { 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f,2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f };
	
	rSeed.AVXIntVec = _mm512_load_epi32( seeds );
	rC0.AVXIntVec = _mm512_load_epi32(&c0); //m512i
	rC1.AVXIntVec = _mm512_load_epi32(&c1); //m512i
	factor.AVXVec = _mm512_load_ps(factors);


#endif // AVX512


	
}

void maxLocStep(Vector &oldWeights, Vector &oldIndices, Vector newWeights, Vector newIndices)
{
#ifdef SISD
	Vector maxMask = gtMask(newWeights, oldWeights);
	oldWeights = mask_mov(oldWeights, maxMask, newWeights);
	oldIndices = mask_mov(oldIndices, maxMask, newIndices);
#elif defined AVX512
	Vector maxMask;
	maxMask.maskVec = _mm512_cmp_ps_mask(newWeights.AVXVec, oldWeights.AVXVec, _MM_CMPINT_GT);
	oldWeights.AVXVec = _mm512_mask_mov_ps(oldWeights.AVXVec, maxMask.maskVec, newWeights.AVXVec);
	oldIndices.AVXVec = _mm512_mask_mov_ps(oldIndices.AVXVec, maxMask.maskVec, newIndices.AVXVec);
#endif // AVX512



}

int reduceMax(Vector curWeights, Vector curIndices)
{
#ifdef SISD

	int highestIndex = -1;
	float highestValue = 0;
	for (int i = 0; i < 16; i++)
	{
		if (curWeights.values[i] > highestValue)
		{
			highestValue = curWeights.values[i];
			highestIndex = curIndices.values[i];
		}
	}
	return highestIndex;
#elif defined AVX512
	// return a vector with all elements equal to ivec[imax] where
	// valvec[imax] is largest element of valvec
	__m512 permVal;
	__m512 permInd;
	__mmask16 maxMask;
	// swap with neighbour 
	permVal = _mm512_swizzle_ps(curWeights.AVXVec, _MM_SWIZ_REG_CDAB);
	permInd = _mm512_swizzle_ps(curIndices.AVXVec, _MM_SWIZ_REG_CDAB);
	maxMask = _mm512_cmp_ps_mask(curWeights.AVXVec, permVal, _MM_CMPINT_GT);
	curWeights.AVXVec = _mm512_mask_mov_ps(permVal, maxMask, curWeights.AVXVec);
	curIndices.AVXVec = _mm512_mask_mov_ps(permInd, maxMask, curIndices.AVXVec);
	// swap pairs
	permVal = _mm512_swizzle_ps(curWeights.AVXVec, _MM_SWIZ_REG_BADC);
	permInd = _mm512_swizzle_ps(curIndices.AVXVec, _MM_SWIZ_REG_BADC);
	maxMask = _mm512_cmp_ps_mask(curWeights.AVXVec, permVal, _MM_CMPINT_GT);
	curWeights.AVXVec = _mm512_mask_mov_ps(permVal, maxMask, curWeights.AVXVec);
	curIndices.AVXVec = _mm512_mask_mov_ps(permInd, maxMask, curIndices.AVXVec);
	// swap lanes
	permVal = _mm512_permute4f128_ps(curWeights.AVXVec, 0xB1); // 2, 3, 0, 1
	permInd = _mm512_permute4f128_ps(curIndices.AVXVec, 0xB1);
	maxMask = _mm512_cmp_ps_mask(curWeights.AVXVec, permVal, _MM_CMPINT_GT);
	curWeights.AVXVec = _mm512_mask_mov_ps(permVal, maxMask, curWeights.AVXVec);
	curIndices.AVXVec = _mm512_mask_mov_ps(permInd, maxMask, curIndices.AVXVec);
	// swap pairs of lanes
	permVal = _mm512_permute4f128_ps(curWeights.AVXVec, 0x4E); // 1, 0, 3, 2
	permInd = _mm512_permute4f128_ps(curIndices.AVXVec, 0x4E);
	maxMask = _mm512_cmp_ps_mask(curWeights.AVXVec, permVal, _MM_CMPINT_GT);
	curIndices.AVXVec = _mm512_mask_mov_ps(permInd, maxMask, curIndices.AVXVec);
	// all elements of ivec now contain index of maximum
	return curIndices.AVXVec[0];
#endif

}

void store(float* loc, Vector v1)
{
#ifdef SISD

	for (int i = 0; i < 16; i++)
	{
		loc[i] = v1.values[i];
	}
	
#elif defined AVX512
	_mm512_store_ps(loc, v1.AVXVec);
#endif
}
