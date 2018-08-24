
#include <immintrin.h>

class AntSystemHelp;


inline void printVectorMM512(__m512 vector)
{
	printf("Printing MM512 Vector ***************\n");
	__declspec(align(64)) float rVals[16]; // random numbers
	_mm512_store_ps( rVals, vector );

	for ( int i = 0; i < 16; i++ )
		printf("%e ",rVals[i] );
	printf("\n");
}

inline void printVectorMM512i(__m512i vector)
{
	printf("Printing MM512i Vector ***************\n");
	__declspec(align(64)) int rVals[16]; // random numbers
	_mm512_store_epi32( rVals, vector );

	for ( int i = 0; i < 16; i++ )
		printf("%u ",rVals[i] );
	printf("\n");
}

inline void printVectorMM512Decimal(__m512 vector)
{
	printf("Printing MM512 Vector ***************\n");
	__declspec(align(64)) float rVals[16]; // random numbers
	_mm512_store_ps( rVals, vector );

	for ( int i = 0; i < 16; i++ )
		printf("%f ",rVals[i] );
	printf("\n");
}

inline void printVectorMask(__mmask16 vector)
{
	printf("Printing Mask Vector ***************\n");
	int maskInt; // random numbers
	maskInt = _mm512_mask2int( vector );
	printf("%d ",maskInt );
	printf("\n");
}
