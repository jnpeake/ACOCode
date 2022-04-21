#pragma once

#include <immintrin.h>
#include <cstdio> 
#include "platform.h"


class Vector
{
public:
#ifdef SISD
	float values[8];
	int iValues[8];
	int vectorSize = 8;
#elif defined AVX512
	__m512 AVXVec;
	__m512i AVXIntVec;
	__mmask16 maskVec;
#elif (defined AVX || defined AVX2)
	__m256 AVXVec;
	__m256i AVXIntVec;
	int mask;
#endif
	Vector operator+(const Vector& v) const
	{
		Vector vector;
#ifdef SISD

		for (int i = 0; i < vectorSize; i++)
		{
			vector.values[i] = this->values[i] + v.values[i];
		}

		return vector;
#elif defined AVX512
		vector.AVXVec = _mm512_add_ps(this->AVXVec, v.AVXVec);
		return vector;
#elif (defined AVX || defined AVX2)
		vector.AVXVec = _mm256_add_ps(this->AVXVec, v.AVXVec);
		return vector;

#endif
	}


	Vector operator-(const Vector& v) const
	{
		Vector vector;
#ifdef SISD
		for (int i = 0; i < vectorSize; i++)
		{
			vector.values[i] = (this->values[i] - v.values[i]);
		}

		return vector;
#elif defined AVX512
		vector.AVXVec = _mm512_sub_ps(this->AVXVec, v.AVXVec);
#elif (defined AVX || defined AVX2)
		vector.AVXVec = _mm256_sub_ps(this->AVXVec, v.AVXVec);
		return vector;

#endif
	}

	Vector operator*(const Vector& v) const
	{
		Vector vector;
#ifdef SISD	
		for (int i = 0; i < vectorSize; i++)
		{
			vector.values[i] = (this->values[i] * v.values[i]);
		}

		return vector;
#elif defined AVX512
		vector.AVXVec = _mm512_mul_ps(this->AVXVec, v.AVXVec);
		return vector;
#elif (defined AVX || defined AVX2)
		vector.AVXVec = _mm256_mul_ps(this->AVXVec, v.AVXVec);
		return vector;

#endif
	}

	Vector operator/(const Vector& v) const
	{
		Vector vector;
#ifdef SISD	
		for (int i = 0; i < vectorSize; i++)
		{
			vector.values[i] = (this->values[i] / v.values[i]);
		}

		return vector;
#elif defined AVX512
		vector.AVXVec = _mm512_div_ps(this->AVXVec, v.AVXVec);
		return vector;
#elif (defined AVX || defined AVX2)
		vector.AVXVec = _mm256_div_ps(this->AVXVec, v.AVXVec);
		return vector;

#endif
	}

	void vecMax(float maxVal)
	{
#ifdef SISD

		for (int i = 0; i < vectorSize; i++)
		{
			if (this->values[i] > maxVal)
			{
				this->values[i] = maxVal;
			}
		}
#elif defined AVX512
		__declspec(align(64)) float maxValAr[16] = { maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal };
		__m512 maxValVec = _mm512_load_ps(maxValAr);
		this->AVXVec = _mm512_min_ps(this->AVXVec, maxValVec);
#elif (defined AVX || defined AVX2)
		float maxValAr[8] __attribute__ ((aligned (32))) = { maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal};
		__m256 maxValVec = _mm256_load_ps(maxValAr);
		this->AVXVec = _mm256_min_ps(this->AVXVec, maxValVec);
#endif

	}

	void vecMin(float minVal)
	{
#ifdef SISD

		for (int i = 0; i < vectorSize; i++)
		{
			if (this->values[i] < minVal)
			{
				this->values[i] = minVal;
			}
		}
#elif defined AVX512
		__declspec(align(64)) float minValAr[16] = { minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal };
		__m512 minValVec = _mm512_load_ps(minValAr);
		this->AVXVec = _mm512_max_ps(this->AVXVec, minValVec);
#elif (defined AVX || defined AVX2)
		float minValAr[8] __attribute__ ((aligned (32))) = { minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal};
		__m256 minValVec = _mm256_load_ps(minValAr);
		this->AVXVec = _mm256_max_ps(this->AVXVec, minValVec);
#endif
	}


	void set1(float setValue)
	{
#ifdef SISD
		for (int i = 0; i < vectorSize; i++)
		{
			this->values[i] = setValue;
		}
#elif defined AVX512
		this->AVXVec = _mm512_set1_ps(setValue);
#elif (defined AVX || defined AVX2)
		this->AVXVec = _mm256_set1_ps(setValue);

#endif
	}


	void load(float* source)
	{
#ifdef SISD

		for (int i = 0; i < vectorSize; i++)
		{

			this->values[i] = source[i];


		}
#elif defined AVX512
		this->AVXVec = _mm512_load_ps(source);
#elif (defined AVX || defined AVX2)
		this->AVXVec = _mm256_load_ps(source);

#endif
	}

	void loadMask(float* source)
	{
#ifdef SISD

		for (int i = 0; i < vectorSize; i++)
		{

			this->values[i] = source[i];


		}
#elif defined AVX512
		this->AVXVec = _mm512_load_ps(source);
#elif (defined AVX || defined AVX2)
		this->AVXVec = _mm256_load_ps(source);

#endif
	}
	
	void load(int* source) //NO AVX  ------------------
	{
#ifdef SISD

		for (int i = 0; i < vectorSize; i++)
		{
			this->values[i] = source[i];
		}
#elif defined AVX512
		this->AVXIntVec = _mm512_load_epi32(source);
#elif (defined AVX || defined AVX2)
		this->AVXIntVec = _mm256_load_si256(((const __m256i *)source));
#endif
	}

	void load(unsigned int * source) //NO AVX --------------
	{
#ifdef SISD

		for (int i = 0; i < vectorSize; i++)
		{
			this->values[i] = source[i];
		}
#elif defined AVX512
		this->AVXIntVec = _mm512_load_epi32(source);
#elif defined AVX2
#endif
	}
	
	void maskLoad(Vector &src, Vector &mask, float* mem_addr) //AVX SIGNIFICANTLY DIFFERENT ------------------
	{
#ifdef SISD

		for (int i = 0; i < vectorSize; i++)
		{
			if (mask.values[i] == 1)
				this->values[i] = mem_addr[i];

			else
				this->values[i] = src.values[i];
		}
#elif defined AVX512
		this->AVXVec = _mm512_mask_load_ps(src.AVXVec, mask.maskVec, mem_addr);
#elif (defined AVX || defined AVX2)

		Vector resultVector;
		resultVector.AVXVec = _mm256_load_ps(mem_addr);
		this->AVXVec = _mm256_blendv_ps(resultVector.AVXVec, src.AVXVec, mask.AVXVec);

#endif
	}

	


private:


};

Vector int2mask(int maskInt);
Vector mask_mov(const Vector& v1, const Vector& bitMask, const Vector& v2);
Vector gtMask(const Vector& v1, const Vector& v2);
Vector ltMask(const Vector& v1, const Vector& v2);
Vector vecRandom(Vector& rC0, Vector& rC1, Vector& factors, Vector& rSeed);
void seedVecRandom(Vector& rC0, Vector& rC1, Vector& factors, int *seeds, Vector& rSeed);
void maxLocStep(Vector &oldWeights, Vector &oldIndices, Vector &newWeights, Vector &newIndices);
int reduceMax(Vector &curWeights, Vector &curIndices);
void store(float* loc, const Vector& v1);
void store(int* loc, const Vector& v1);
void printVec(Vector v1);