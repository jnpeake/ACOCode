#include <cstdio> 
#include <cstdlib>
#include <iostream>
#include "vector.h"
#include "antsystemhelp.h"

class VectorTesting
{
	public void compareVec512(Vector eVec, Vector AVXVec)
	{
		float[16] storedVec;
		store(storedVec, AVXVec);
	
		for(i = 0; i < 16; i++)
		{
			if(storedVec[i] == eVec[i])
			{
				printf("ERROR");
			}
		}

	}
		
	int main( int argc, char *argv[] )
	{
		__declspec(align(64)) float[16] testArr = { 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f,
		8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f};
		
		Vector eVec;
		Vector AVXVec;
		
		eVec.Values = { 0.0f, 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f,
		8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f};
		
		load(AVXVec.AVXVec, testArr);
		
	}
	
}
