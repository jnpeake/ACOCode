#include <iostream>
#include "Vector.h"


using namespace std;

void plusTest()
{
	Vector testVector1;
	Vector testVector2;
	Vector testVector3;
	std::vector<int> testVec1;
	std::vector<int> testVec2;

	for (int i = 0; i < 16; i++)
	{
		testVec1.push_back(i);
		testVec2.push_back(i * 3);
	}

	//testVector1.values = testVec1;
	//testVector2.values = testVec2;

	for (int i = 0; i < 16; i++)
	{
		cout << " " << testVec1[i];
	}
	cout << "\n";
	for (int i = 0; i < 16; i++)
	{
		cout << " " << testVec2[i];
	}

	testVector3 = testVector1 + testVector2;

	cout << "\n";
	for (int i = 0; i < 16; i++)
	{
		cout << " " << testVector3.values[i];
	}
}

void minusTest()
{
	Vector testVector1;
	Vector testVector2;
	Vector testVector3;
	std::vector<int> testVec1;
	std::vector<int> testVec2;

	for (int i = 0; i < 16; i++)
	{
		testVec1.push_back(i);
		testVec2.push_back(i * 3);
	}

	//testVector1.values = testVec1;
	//testVector2.values = testVec2;

	for (int i = 0; i < 16; i++)
	{
		cout << " " << testVec1[i];
	}
	cout << "\n";
	for (int i = 0; i < 16; i++)
	{
		cout << " " << testVec2[i];
	}

	testVector3 = testVector1 - testVector2;

	cout << "\n";
	for (int i = 0; i < 16; i++)
	{
		cout << " " << testVector3.values[i];
	}
}

void multTest()
{
	Vector testVector1;
	Vector testVector2;
	Vector testVector3;
	int testVec1[16];
	int testVec2[16];

	for (int i = 0; i < 16; i++)
	{
		testVector1.values[i] = i;
		testVector2.values[i] = (i * 3);
	}

	testVector3 = testVector1 * testVector2;

	cout << "\n";
	for (int i = 0; i < 16; i++)
	{
		cout << " " << testVector3.values[i];
	}
}

void extLoadTest()
{
	Vector testVector1;

	testVector1.extload(5);


	for (int i = 0; i < 16; i++)
	{
		cout << " " << testVector1.values[i];
	}
}

void gtMaskTest()
{
	Vector testVector1;
	for (int i = 0; i < 16; i++)
	{
		testVector1.values[i] = i;
	}
	Vector testVector2;

	for (int i = 0; i < 16; i++)
	{
		if (i == 0 || i == 5 || i == 9 | i == 15)
		{
			testVector2.values[i] = 100;
		}

		else
		{
			testVector2.values[i] = 0;
		}
	}


	Vector resultMask = gtMask(testVector1, testVector2);
	for (int i = 0; i < 16; i++)
	{
		cout << " " << resultMask.values[i];
	}

}

void ltMaskTest()
{
	Vector testVector1;
	for (int i = 0; i < 16; i++)
	{
		testVector1.values[i] = i;
	}
	Vector testVector2;

	for (int i = 0; i < 16; i++)
	{
		if (i % 2 == 0)
			testVector2.values[i] = i * 2;
		else
			testVector2.values[i] = 0;
	}

	Vector resultMask = ltMask(testVector1, testVector2);
	for (int i = 0; i < 16; i++)
	{
		cout << " " << resultMask.values[i];
	}

}

void maxTest()
{
	Vector testVector1;
	for (int i = 0; i < 16; i++)
	{
		testVector1.values[i] = i;
	}

	testVector1.max(5);
	for (int i = 0; i < 16; i++)
	{
		cout << " " << testVector1.values[i];
	}

}

void minTest()
{
	Vector testVector1;
	for (int i = 0; i < 16; i++)
	{
		testVector1.values[i] = i;
	}

	testVector1.min(5);
	for (int i = 0; i < 16; i++)
	{
		cout << " " << testVector1.values[i];
	}

}

void maskTest()
{
	Vector maskVector = int2mask(1873);
	for (int i = 15; i >= 0; i--)
	{
		cout << " " << maskVector.values[i];
	}

}


void setTest()
{
	Vector testVector;
	testVector.set1(50);
	for (int i = 0; i < 16; i++)
	{
		cout << " " << testVector.values[i];
	}

}

void maskMovTest()
{
	int mask[] = { 0,0,0,0,0,1,1,0,0,0,0,1,1,1,1,0 };

	Vector testVector1;
	Vector testVector2;
	for (int i = 0; i < 16; i++)
	{
		testVector1.values[i] = i;
	}

	testVector2.set1(99);

	Vector resultVector = mask_mov(testVector1, mask, testVector2);
	for (int i = 0; i < 16; i++)
	{
		cout << " " << resultVector.values[i];
	}

}

void storeTest()
{
	Vector testVector;
	testVector.set1(12);

	float* testArray = (float*)malloc((16 * sizeof(float)));

	testArray = testVector.store();

	for (int i = 0; i < 16; i++)
	{
		cout << " " << testArray[i];
	}
}



void loadTest()
{
	Vector testVector1;
	float* testArray;
	testArray = (float*)malloc((16 * sizeof(float)));

	for (int i = 0; i < 16; i++)
	{
		testArray[i] = i;
	}

	testVector1.load(testArray);
	for (int i = 0; i < 16; i++)
	{
		cout << " " << testVector1.values[i];
	}
}

void reduceTest()
{
	Vector testVector1;
	for (int i = 0; i < 16; i++)
	{
		testVector1.values[i] = i*4;
	}

	testVector1.values[5] = 100000;

	int highestIndex = testVector1.reduceMax();

	cout << highestIndex;
}

void randomTest()
{
	Vector testVector1 = random();

	for (int i = 0; i < 16; i++)
	{
		cout << " " << testVector1.values[i];
	}

}


int main()
{
	randomTest();
	return 0;
}
