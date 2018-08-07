#include "vector.h"
#include "ant_xp.h"

Vector int2mask(int maskInt)
{
#ifdef _WIN32

	Vector mask;
	for (int i = 0; i < 16; ++i) {
		mask.values[i] = (maskInt >> i) & 1;
	}


	return mask;
#endif

#ifdef AVX512
	resultVector = _mm512_int2mask(maskInt);
	return resultVector;
#endif // AVX512

}

Vector mask_mov(const Vector& v1, Vector bitMask, const Vector& v2)
{
	Vector resultVector;
#ifdef _WIN32
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
#endif

#ifdef AVX512
	resultVector = _mm512_mask_mov_ps(v1, bitmask, v2);
	return resultVector;
#endif // AVX512
}

Vector gtMask(const Vector& v1, const Vector& v2)
{
#ifdef _WIN32

	Vector resultMask;

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
#endif

#ifdef AVX512
	return _mm_cmp_ps_mask(v1, v2, _MM_CMPINT_GT)
#endif // AVX512
}

Vector ltMask(const Vector& v1, const Vector& v2)
{
#ifdef _WIN32

	Vector resultMask;

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
#endif

#ifdef AVX512
	return _mm_cmp_ps_mask(v1,v2, _MM_CMPINT_LT)
#endif // AVX512

}

Vector vecRandom(unsigned int* rSeed)
{
#ifdef _WIN32
	Vector r;
	for (int i = 0; i < 16; i++)
	{
		rSeed[i] = rSeed[i] * 1664525L + 1013904223L;
		r.values[i] = (float)rSeed[i] * 2.328306437087974e-10;
	}

	return r;
#endif // _WIN32

#ifdef AVX512
	__declspec(align(64)) unsigned int c0[16] = { 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L,1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L };
	__declspec(align(64)) unsigned int c1[16] = { 1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L, 	1013904223L, 1013904223L, 1013904223L,1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L };
	__declspec(align(64)) float factors[16] = { 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f,2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f };

	rC0 = _mm512_load_epi32(&c0); //m512i
	rC1 = _mm512_load_epi32(&c1); //m512i
	factor = _mm512_load_ps(factors);
#endif // AVX512

	
}

void seedVecRandom(int *seeds, unsigned int *rSeed)
{
#ifdef _WIN32
	memcpy(rSeed, seeds, 16 * sizeof(int));
#endif // DEBUG

#ifdef AVX512
	__declspec(align(64)) unsigned int c0[16] = { 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L,1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L };
	__declspec(align(64)) unsigned int c1[16] = { 1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L, 	1013904223L, 1013904223L, 1013904223L,1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L };
	__declspec(align(64)) float factors[16] = { 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f,2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f };

	rC0 = _mm512_load_epi32(&c0); //m512i
	rC1 = _mm512_load_epi32(&c1); //m512i
	factor = _mm512_load_ps(factors);
#endif // AVX512


	
}

void maxLocStep(Vector &oldWeights, Vector &oldIndices, Vector &newWeights, Vector &newIndices)
{
#ifdef _WIN32
	Vector maxMask = gtMask(newWeights, oldWeights);
	oldWeights = mask_mov(oldWeights, maxMask, newWeights);
	oldIndices = mask_mov(oldIndices, maxMask, newIndices);
#endif // _WIN32

#ifdef AVX512
	__mmask16 maxMask = _mm512_cmp_ps_mask(newWeights, oldWeights, _MM_CMPINT_GT);
	oldWeights = _mm512_mask_mov_ps(oldWeights, maxMask, newWeights);
	oldIndices = _mm512_mask_mov_ps(oldIndices, maxMask, newIndices);
#endif // AVX512



}

int reduceMax(Vector curWeights, Vector curIndices)
{
#ifdef _WIN32

	int highestIndex = -1;
	float highestValue = -1;
	for (int i = 0; i < 16; i++)
	{
		if (curWeights.values[i] > highestValue)
		{
			highestValue = curWeights.values[i];
			highestIndex = i;
		}
	}
	return curIndices.values[highestIndex];
#endif 

#ifdef AVX512
	// return a vector with all elements equal to ivec[imax] where
	// valvec[imax] is largest element of valvec
	__m512 permVal;
	__m512 permInd;
	__mmask16 maxMask;
	// swap with neighbour 
	permVal = _mm512_swizzle_ps(curWeights, _MM_SWIZ_REG_CDAB);
	permInd = _mm512_swizzle_ps(curIndices, _MM_SWIZ_REG_CDAB);
	maxMask = _mm512_cmp_ps_mask(curWeights, permVal, _MM_CMPINT_GT);
	curWeights = _mm512_mask_mov_ps(permVal, maxMask, valvec);
	curIndices = _mm512_mask_mov_ps(permInd, maxMask, curIndices);
	// swap pairs
	permVal = _mm512_swizzle_ps(curWeights, _MM_SWIZ_REG_BADC);
	permInd = _mm512_swizzle_ps(curIndices, _MM_SWIZ_REG_BADC);
	maxMask = _mm512_cmp_ps_mask(curWeights, permVal, _MM_CMPINT_GT);
	curWeights = _mm512_mask_mov_ps(permVal, maxMask, valvec);
	curIndices = _mm512_mask_mov_ps(permInd, maxMask, curIndices);
	// swap lanes
	permVal = _mm512_permute4f128_ps(curWeights, 0xB1); // 2, 3, 0, 1
	permInd = _mm512_permute4f128_ps(curIndices, 0xB1);
	maxMask = _mm512_cmp_ps_mask(curWeights, permVal, _MM_CMPINT_GT);
	curWeights = _mm512_mask_mov_ps(permVal, maxMask, valvec);
	curIndices = _mm512_mask_mov_ps(permInd, maxMask, curIndices);
	// swap pairs of lanes
	permVal = _mm512_permute4f128_ps(curWeights, 0x4E); // 1, 0, 3, 2
	permInd = _mm512_permute4f128_ps(curIndices, 0x4E);
	maxMask = _mm512_cmp_ps_mask(curWeights, permVal, _MM_CMPINT_GT);
	curIndices = _mm512_mask_mov_ps(permInd, maxMask, curIndices);
	// all elements of ivec now contain index of maximum
	return curIndices[0];
#endif

}

void store(float* loc, Vector v1)
{
#ifdef _WIN32

	for (int i = 0; i < 16; i++)
	{
		loc[i] = v1.values[i];
	}
#endif

#ifdef AVX512
	_mm512_store_ps(v1, loc);
#endif
}