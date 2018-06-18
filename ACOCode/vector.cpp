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
}

Vector mask_mov(const Vector& v1, Vector bitMask, const Vector& v2)
{
#ifdef _WIN32

	Vector resultVector;

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
}

Vector vecRandom(unsigned int* rSeed)
{
	Vector r;
	for (int i = 0; i < 16; i++)
	{
		rSeed[i] = rSeed[i] * 1664525L + 1013904223L;
		r.values[i] = (float)rSeed[i] * 2.328306437087974e-10;
	}

	return r;
}

void seedVecRandom(int *seeds, unsigned int *rSeed)
{
	memcpy(rSeed, seeds, 16 * sizeof(int));
}

void maxLocStep(Vector &oldWeights, Vector &oldIndices, Vector &newWeights, Vector &newIndices)
{
	Vector maxMask = gtMask(newWeights, oldWeights);
	oldWeights = mask_mov(oldWeights, maxMask, newWeights);
	oldIndices = mask_mov(oldIndices, maxMask, newIndices);
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
}

void store(float* loc, Vector v1)
{
#ifdef _WIN32

	for (int i = 0; i < 16; i++)
	{
		loc[i] = v1.values[i];
	}
#endif
}