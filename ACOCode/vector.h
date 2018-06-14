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

	void max(float maxVal)
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

	void min(float minVal)
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

	/*
	void swizzle(std::string swizzleParam)
	{
	#ifdef _WIN32
	float A[] = { this->values[0], this->values[1], this->values[2], this->values[3] };
	float B[] = { this->values[4], this->values[5], this->values[6], this->values[7] };
	float C[] = { this->values[8], this->values[9], this->values[10], this->values[11] };
	float D[] = { this->values[12], this->values[13], this->values[14], this->values[15] };
	float result[16];
	if (swizzleParam == "CDAB")
	{
	std::copy(C, C + 4, result);
	std::copy(D, D + 4, result + 4);
	std::copy(A, A + 4, result + 8);
	std::copy(B, B + 4, result + 12);
	}
	else if (swizzleParam == "BADC")
	{
	std::copy(B, B + 4, result);
	std::copy(A, A + 4, result + 4);
	std::copy(D, D + 4, result + 8);
	std::copy(C, C + 4, result + 12);
	}
	for (int i = 0; i < 16; i++)
	{
	this->values[i] = result[i];
	}
	#endif

	}
	void permute(int imm8)
	{
	#ifdef _WIN32
	float result[16];
	for (int i = 0; i < 4; i++)
	{
	for (int j = 0; j < 4; j++)
	{
	//imm8 = this->values[i];
	}
	}
	#endif
	}
	*/

	int reduceMax()
	{
#ifdef _WIN32

		int highestIndex = -1;
		int highestValue = -1;
		for (int i = 0; i < 16; i++)
		{
			if (this->values[i] > highestValue)
			{
				highestValue = this->values[i];
				highestIndex = i;
			}
		}

		return highestIndex;

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

	float* store()
	{
#ifdef _WIN32

		float* storeValues = (float*)malloc((16 * sizeof(float)));

		for (int i = 0; i < 16; i++)
		{
			storeValues[i] = this->values[i];
		}

		return storeValues;
#endif
	}


private:


};

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
		if (bitMask.values[i] == 1)
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

Vector random(int* seeds)
{
	Vector randomVector;

	for (int i = 0; i < 16; i++)
	{
		srand(seeds[i]);
		randomVector.values[i] = ((double)rand() / (RAND_MAX));
	}

	return randomVector;

}
