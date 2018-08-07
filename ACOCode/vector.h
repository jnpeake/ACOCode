#pragma once

#include <vector>


class Vector
{
public:
	float values[16];
	int vectorSize = 16;
	Vector operator+(const Vector& v) const
	{
#ifdef _WIN32
		Vector vector;
		for (int i = 0; i < vectorSize; i++)
		{
			vector.values[i] = this->values[i] + v.values[i];
		}

		return vector;
#endif

#ifdef AVX512
		_mm512_add_ps(this->values[i] + v.values[i]);
#endif

#ifdef AVX2

#endif
	}

	Vector operator-(const Vector& v) const
	{
#ifdef _WIN32

		Vector vector;
		for (int i = 0; i < vectorSize; i++)
		{
			vector.values[i] = (this->values[i] - v.values[i]);
		}

		return vector;
#endif
#ifdef AVX512
		_mm512_sub_ps(this->values[i] + v.values[i]);
#endif

#ifdef AVX2

#endif
	}

	Vector operator*(const Vector& v) const
	{
#ifdef _WIN32

		Vector vector;
		for (int i = 0; i < vectorSize; i++)
		{
			vector.values[i] = (this->values[i] * v.values[i]);
		}

		return vector;
#endif
#ifdef AVX512
		_mm512_mul_ps(this->values[i] + v.values[i]);
#endif

#ifdef AVX2

#endif
	}



	void extload(float value)
	{
#ifdef _WIN32

		for (int i = 0; i < vectorSize; i++)
		{
			this->values[i] = (value);
		}
#endif
#ifdef AVX512
		this->values[i] = _mm512_extload_ps(value, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, 0);
#endif

#ifdef AVX2

#endif
	}

	void vecMax(float maxVal)
	{
#ifdef _WIN32

		for (int i = 0; i < vectorSize; i++)
		{
			if (this->values[i] < maxVal)
			{
				this->values[i] = maxVal;
			}
		}
#endif
#ifdef AVX512
		this->values[i] = _mm512_max_ps(this->values[i], maxVal)
#endif

#ifdef AVX2

#endif
	}

	void vecMin(float minVal)
	{
#ifdef _WIN32

		for (int i = 0; i < vectorSize; i++)
		{
			if (this->values[i] > minVal)
			{
				this->values[i] = minVal;
			}
		}
#endif
#ifdef AVX512
		this->values[i] = _mm512_max_ps(this->values[i], minval)
#endif

#ifdef AVX2

#endif
	}


	void set1(float setValue)
	{
#ifdef _WIN32

		for (int i = 0; i < vectorSize; i++)
		{
			this->values[i] = setValue;
		}
#endif
#ifdef AVX512
		this->values[i] = _mm512_set1_ps(setValue);
#endif

#ifdef AVX2

#endif
	}


	void load(float* source)
	{
#ifdef _WIN32

		for (int i = 0; i < vectorSize; i++)
		{

				this->values[i] = source[i];


		}
#endif
#ifdef AVX512
		this->values[i] = _mm512_load_ps(source);
#endif

#ifdef AVX2

#endif
	}

	void load(int* source)
	{
#ifdef _WIN32

		for (int i = 0; i < vectorSize; i++)
		{
			this->values[i] = source[i];
		}
#endif
#ifdef AVX512
		this->values[i] = _mm512_load_epi32(source);
#endif

#ifdef AVX2

#endif
	}

	void load(unsigned int * source)
	{
#ifdef _WIN32

		for (int i = 0; i < vectorSize; i++)
		{
			this->values[i] = source[i];
		}
#endif
#ifdef AVX512
		this->values[i] = _mm512_load_epi32(source);
#endif

#ifdef AVX2

#endif
	}

	void maskLoad(Vector src, Vector mask, float* mem_addr)
	{
#ifdef _WIN32

		for (int i = 0; i < vectorSize; i++)
		{
			if (mask.values[i] == 1)
				this->values[i] = mem_addr[i];

			else
				this->values[i] = src.values[i];
		}
#endif
#ifdef AVX512
		_mm512_mask_load_ps(src, mask, mem_addr);
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

Vector vecRandom(unsigned int* rSeed);

void seedVecRandom(int* seeds, unsigned int* rSeed);

void maxLocStep(Vector &oldWeights, Vector &oldIndices, Vector &newWeights, Vector &newIndices);
int reduceMax(Vector curWeights, Vector curIndices);

void store(float* loc, Vector v1);
