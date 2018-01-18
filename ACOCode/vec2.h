#ifndef _VEC2_INC_
#define _VEC2_INC_
// a 2-d vector class
#include <math.h>

class V2
{
public:
	float x;
	float y;

	V2()
	{
		x = y = 0.0f;
	}
	~V2(){}

	V2( float nx, float ny )
	{
		x = nx;
		y = ny;
	}

	V2 operator+( const V2 &v )
	{
		return V2( x+v.x, y+v.y );
	}

    V2 operator-( const V2 &v )
	{
		return V2( x-v.x, y-v.y );
	}


	V2 operator*( float f )
	{
		return V2( x*f, y*f );
	}

	float operator*( const V2 &v )
	{
		// dot product
		return x*v.x + y*v.y;
	}

	void operator~()
	{
		// normalize
		float mag = magnitude();
		if ( mag != 0.0f )
		{
			float iMag = 1.0f/magnitude();
			x *= iMag;
			y *= iMag;
		}
	}

	float magnitude( void )
	{
		return sqrt( x*x + y*y );
	}

	float magnsqrd( void )
	{
		return x*x + y*y;
	}
};
#endif
