#pragma once

#include <immintrin.h>
#include <cstdio> 


#define SISD
//#define AVX512

class Vector
{
public:
#ifdef SISD
	float values[16];
	unsigned int iValues[16];
	int vectorSize = 16;
#elif defined AVX512
	__m512 AVXVec;
	__m512i AVXIntVec;
	__mmask16 maskVec;
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
#endif

#ifdef AVX2

#endif
	}

	Vector operator-(const Vector& v) const
	{
#ifdef SISD

		Vector vector;
		for (int i = 0; i < vectorSize; i++)
		{
			vector.values[i] = (this->values[i] - v.values[i]);
		}

		return vector;
#elif defined AVX512
		_mm512_sub_ps(this->AVXVec, v.AVXVec);
#endif

#ifdef AVX2

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
#endif

#ifdef AVX2

#endif
	}



	void extload(float* value)
	{
#ifdef SISD

		for (int i = 0; i < vectorSize; i++)
		{
			this->values[i] = value[0];
		}
#elif defined AVX512
		this->AVXVec = _mm512_extload_ps(value, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, 0);
#endif

#ifdef AVX2

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
		__declspec(align(64)) float maxValAr[16] ={maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal,maxVal};
		__m512 maxValVec = _mm512_load_ps(maxValAr);
		this->AVXVec = _mm512_min_ps(this->AVXVec, maxValVec);
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
		__declspec(align(64)) float minValAr[16] ={minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal,minVal};
		__m512 minValVec = _mm512_load_ps(minValAr);
		this->AVXVec = _mm512_max_ps(this->AVXVec, minValVec);
#endif

#ifdef AVX2

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
#endif

#ifdef AVX2

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
#endif

#ifdef AVX2

#endif
	}

	void load(int* source)
	{
#ifdef SISD

		for (int i = 0; i < vectorSize; i++)
		{
			this->values[i] = source[i];
		}
#elif defined AVX512
		this->AVXIntVec = _mm512_load_epi32(source);
#endif

#ifdef AVX2

#endif
	}

	void load(unsigned int * source)
	{
#ifdef SISD

		for (int i = 0; i < vectorSize; i++)
		{
			this->values[i] = source[i];
		}
#elif defined AVX512
		this->AVXIntVec = _mm512_load_epi32(source);
#endif

#ifdef AVX2

#endif
	}

	void maskLoad(Vector &src, Vector mask, float* mem_addr)
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
#endif

#ifdef AVX2

#endif
	}

	


private:


};

Vector int2mask(int maskInt);

Vector mask_mov(const Vector& v1, Vector bitMask, const Vector& v2);

Vector gtMask(const Vector& v1, const Vector& v2);

Vector ltMask(const Vector& v1, const Vector& v2);

Vector vecRandom(Vector rC0, Vector rC1, Vector factors, Vector rSeed);

void seedVecRandom(Vector& rC0, Vector& rC1, Vector& factors, int *seeds, Vector& rSeed);

void maxLocStep(Vector &oldWeights, Vector &oldIndices, Vector newWeights, Vector newIndices);
int reduceMax(Vector curWeights, Vector curIndices);

void store(float* loc, Vector v1);
