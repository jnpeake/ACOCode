#pragma once

#include <vector>


class Vector
{
public:
	float values[16];
	Vector operator+(const Vector& v) const
	{
#ifdef _WIN32
		Vector vector;
		for (int i = 0; i < 16; i++)
		{
			vector.values[i] = this->values[i] + v.values[i];
		}

		return vector;
#endif
	}

	Vector operator-(const Vector& v) const
	{
#ifdef _WIN32

		Vector vector;
		for (int i = 0; i < 16; i++)
		{
			vector.values[i] = (this->values[i] - v.values[i]);
		}

		return vector;
#endif
	}

	Vector operator*(const Vector& v) const
	{
#ifdef _WIN32

		Vector vector;
		for (int i = 0; i < 16; i++)
		{
			vector.values[i] = (this->values[i] * v.values[i]);
		}

		return vector;
#endif
	}



	void extload(float value)
	{
#ifdef _WIN32

		for (int i = 0; i < 16; i++)
		{
			this->values[i] = (value);
		}
#endif
	}

	void vecMax(float maxVal)
	{
#ifdef _WIN32

		for (int i = 0; i < 16; i++)
		{
			if (this->values[i] < maxVal)
			{
				this->values[i] = maxVal;
			}
		}
#endif
	}

	void vecMin(float minVal)
	{
#ifdef _WIN32

		for (int i = 0; i < 16; i++)
		{
			if (this->values[i] > minVal)
			{
				this->values[i] = minVal;
			}
		}
#endif
	}


	void set1(float setValue)
	{
#ifdef _WIN32

		for (int i = 0; i < 16; i++)
		{
			this->values[i] = setValue;
		}
#endif
	}


	void load(float* source)
	{
#ifdef _WIN32

		for (int i = 0; i < 16; i++)
		{

				this->values[i] = source[i];


		}
#endif
	}

	void load(int* source)
	{
#ifdef _WIN32

		for (int i = 0; i < 16; i++)
		{
			this->values[i] = source[i];
		}
#endif
	}

	void load(unsigned int * source)
	{
#ifdef _WIN32

		for (int i = 0; i < 16; i++)
		{
			this->values[i] = source[i];
		}
#endif
	}

	void maskLoad(Vector src, Vector mask, float* mem_addr)
	{
#ifdef _WIN32

		for (int i = 0; i < 16; i++)
		{
			if (mask.values[i] == 1)
				this->values[i] = mem_addr[i];

			else
				this->values[i] = src.values[i];
		}
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
