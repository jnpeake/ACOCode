#include "vector.h"
#include <iostream>
#include "platform.h"


Vector int2mask(int maskInt) // NO AVX --------------------------
{
#ifdef SISD

	Vector mask;
	for (int i = 0; i < _VECSIZE; ++i) {
		mask.values[i] = (maskInt >> i) & 1;
	}


	return mask;

#elif defined AVX512
	Vector resultVector;
	resultVector.maskVec = _mm512_int2mask(maskInt);
	return resultVector;

#elif (defined AVX || defined AVX2)
	Vector resultVector;
	float mask[8];
	for (int i = 0; i < 8; i++) {
		if ((maskInt >> i) & 1)
		{
			mask[i] = 1.0f;
		}

		else
		{
			mask[i] = -1.0f;
		}
	}

	resultVector.AVXVec = _mm256_setr_ps(mask[0], mask[1], mask[2], mask[3], mask[4], mask[5], mask[6], mask[7]);

	return resultVector;

#endif // AVX512

}


Vector mask_mov(const Vector& v1, const Vector& bitMask, const Vector& v2) // NO AVX -----------------
{
	Vector resultVector;
#ifdef SISD
	for (int i = 0; i < _VECSIZE; i++)
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
#elif (defined AVX || defined AVX2)
	resultVector.AVXVec = _mm256_blendv_ps(v2.AVXVec, v1.AVXVec, bitMask.AVXVec);
	return resultVector;
#endif // AVX512
}


Vector gtMask(const Vector& v1, const Vector& v2) // NO AVX -------------------
{
	Vector resultMask;
#ifdef SISD

	for (int i = 0; i < _VECSIZE; i++)
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
#elif (defined AVX || defined AVX2)
	resultMask.AVXVec = _mm256_cmp_ps(v1.AVXVec, v2.AVXVec, _CMP_GT_OS);
	return resultMask;
#endif
}

Vector ltMask(const Vector& v1, const Vector& v2) // NO AVX --------------------------
{
	Vector resultMask;
#ifdef SISD
	for (int i = 0; i < _VECSIZE; i++)
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
	resultMask.maskVec = _mm512_cmp_ps_mask(v1.AVXVec, v2.AVXVec, _MM_CMPINT_LT);
	return resultMask;
#elif (defined AVX || defined AVX2)
	resultMask.AVXVec = _mm256_cmp_ps(v1.AVXVec, v2.AVXVec, _CMP_LT_OS);
	return resultMask;
#endif

}

void seedVecRandom(Vector& rC0, Vector& rC1, Vector& factor, int *seeds, Vector& rSeed) // SIGNIFICANTLY DIFFERENT TO AVX512 -----------------------
{
#ifdef SISD
	for (int i = 0; i < 8; i++)
	{
		rSeed.iValues[i] = seeds[i];
	}


#elif defined AVX512
	__declspec(align(64)) unsigned int c0[16] = { 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L,1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L, 1664525L };
	__declspec(align(64)) unsigned int c1[16] = { 1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L, 	1013904223L, 1013904223L, 1013904223L,1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L, 1013904223L };
	__declspec(align(64)) float factors[16] = { 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f,2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f };

	rSeed.AVXIntVec = _mm512_load_epi32(seeds);
	rC0.AVXIntVec = _mm512_load_epi32(&c0); //m512i
	rC1.AVXIntVec = _mm512_load_epi32(&c1); //m512i
	factor.AVXVec = _mm512_load_ps(factors);


#elif (defined AVX || defined AVX2)
	int c0 = 1664525L;
	int c1 = 1013904223L;
	float factors[16] __attribute__ ((aligned (32))) = { 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f,2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f, 2.3283064e-10f };

	rSeed.AVXIntVec = _mm256_load_si256((__m256i*)seeds);
	rC0.AVXIntVec = _mm256_set1_epi32(c0); //m512i
	rC1.AVXIntVec = _mm256_set1_epi32(c1); //m512i
	factor.AVXVec = _mm256_load_ps(factors);

#endif // AVX512




}

Vector vecRandom(Vector& rC0, Vector& rC1, Vector& factor, Vector& rSeed) // SIGNIFICANTLY DIFFERENT TO AVX512 --------------------------
{
	Vector r;
#ifdef SISD
	for (int i = 0; i < 8; i++)
	{

		rSeed.iValues[i] = rSeed.iValues[i] * 1664525L + 1013904223L;
		r.values[i] = (float)rSeed.iValues[i] * 2.328306437087974e-10;
		r.values[i] += 0.5f;
	}

	return r;

#elif defined AVX512
	rSeed.AVXIntVec = _mm512_mullo_epi32(rC0.AVXIntVec, rSeed.AVXIntVec);
	rSeed.AVXIntVec = _mm512_add_epi32(rC1.AVXIntVec, rSeed.AVXIntVec);
	// convert to float in range 0 to 1 and return
	__m512 returnValue = _mm512_cvt_roundepu32_ps(rSeed.AVXIntVec, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
	r.AVXVec = _mm512_mul_ps(returnValue, factor.AVXVec);
	return r;

#elif defined AVX2
	__m256 addedVal;
	addedVal =_mm256_set1_ps(0.5f);
	rSeed.AVXIntVec = _mm256_mullo_epi32(rC0.AVXIntVec, rSeed.AVXIntVec);
	rSeed.AVXIntVec = _mm256_add_epi32(rC1.AVXIntVec, rSeed.AVXIntVec);
	// convert to float in range 0 to 1 and return
	__m256 returnValue = _mm256_cvtepi32_ps(rSeed.AVXIntVec);
	
	//returnValue = _mm256_round_ps(returnValue, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
	returnValue = _mm256_mul_ps(returnValue, factor.AVXVec);
	r.AVXVec = _mm256_add_ps(returnValue,addedVal);
	return r;
	

#endif // AVX512


}

Vector EucDistanceVec(Vector x0, Vector y0, Vector x1, Vector y1)
	{
		float zeroPointFive = 0.5f;
		Vector zeroPointFiveVec;
		zeroPointFiveVec.AVXVec = _mm256_set1_ps(zeroPointFive);
		Vector xMinus = x0 - y0;
		Vector yMinus = y0 - y1;
		Vector xMinusSq = xMinus * xMinus;
		Vector yMinusSq = yMinus * yMinus;

		Vector d = xMinusSq + yMinusSq;
		d.AVXVec = _mm256_sqrt_ps(d.AVXVec);
		d = d + zeroPointFiveVec;
		Vector convertedVec;
		convertedVec.AVXIntVec = _mm256_cvtps_epi32(d.AVXVec);
		return convertedVec;
	}



void maxLocStep(Vector &oldWeights, Vector &oldIndices, Vector &newWeights, Vector &newIndices) // SIGNIFICANTLY DIFFERENT TO AVX512 --------------------------
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
#elif (defined AVX || defined AVX2)
	Vector maxMask;
	maxMask = gtMask(newWeights, oldWeights);
	oldWeights = mask_mov(newWeights, maxMask, oldWeights);
	oldIndices = mask_mov(newIndices, maxMask, oldIndices);
#endif


}


int reduceMax(Vector &curWeights, Vector &curIndices) // SIGNIFICANTLY DIFFERENT TO AVX512 ------------------------------
{
#ifdef SISD

	int highestIndex = -1;
	float highestValue = -1.0f;
	for (int i = 0; i < _VECSIZE; i++)
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
#elif (defined AVX || defined AVX2)
	
	__m256 permVal;
	__m256 permInd;
	__m256 maxMask;
	float result[8] __attribute__ ((aligned (32)));
	permVal = _mm256_permute_ps(curWeights.AVXVec, _MM_SHUFFLE(2,3,0,1)); //01001110
	permInd = _mm256_permute_ps(curIndices.AVXVec, _MM_SHUFFLE(2, 3, 0, 1)); //01001110
	maxMask = _mm256_cmp_ps(curWeights.AVXVec, permVal, _CMP_GT_OS);
	curWeights.AVXVec = _mm256_blendv_ps(permVal, curWeights.AVXVec, maxMask);
	curIndices.AVXVec = _mm256_blendv_ps(permInd, curIndices.AVXVec, maxMask);

	permVal = _mm256_permute_ps(curWeights.AVXVec, _MM_SHUFFLE(1,0,3,2)); //01001110
	permInd = _mm256_permute_ps(curIndices.AVXVec, _MM_SHUFFLE(1, 0, 3, 2)); //01001110
	maxMask = _mm256_cmp_ps(curWeights.AVXVec, permVal, _CMP_GT_OS);
	curWeights.AVXVec = _mm256_blendv_ps(permVal, curWeights.AVXVec, maxMask);
	curIndices.AVXVec = _mm256_blendv_ps(permInd, curIndices.AVXVec, maxMask);

	permVal = _mm256_permute2f128_ps(curWeights.AVXVec, curWeights.AVXVec, 0b0001);
	permInd = _mm256_permute2f128_ps(curIndices.AVXVec, curIndices.AVXVec, 0b0001);
	maxMask = _mm256_cmp_ps(curWeights.AVXVec, permVal, _CMP_GT_OS);
	curWeights.AVXVec = _mm256_blendv_ps(permVal, curWeights.AVXVec, maxMask);
	curIndices.AVXVec = _mm256_blendv_ps(permInd, curIndices.AVXVec, maxMask);
	store(result, curIndices);
	return result[0];
#endif

}

void store(float* loc, const Vector& v1)
{
#ifdef SISD

	for (int i = 0; i < _VECSIZE; i++)
	{
		loc[i] = v1.values[i];
	}

#elif defined AVX512
	_mm512_store_ps(loc, v1.AVXVec);
#elif (defined AVX || defined AVX2)
	_mm256_store_ps(loc, v1.AVXVec);
#endif
}

void store(int* loc, const Vector& v1)
{
#ifdef SISD

	for (int i = 0; i < _VECSIZE; i++)
	{
		loc[i] = v1.values[i];
	}

#elif defined AVX512
	_mm512_store_ps(loc, v1.AVXVec);
#elif (defined AVX || defined AVX2)
	_mm256_store_si256((__m256i *)loc, v1.AVXIntVec);
#endif
}

#ifdef SISD
void printVec(Vector v1)
{
	std::cout << " \n\n Printing Test Vector \n";
	for (int i = 0; i < v1.vectorSize; i++)
	{
		std::cout << v1.values[i] << " ";
	}
}
#endif
