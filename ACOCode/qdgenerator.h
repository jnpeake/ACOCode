#ifndef _QDGENERATOR_INC_
#define _QDGENERATOR_INC_

class QDGenerator
{
private:
	unsigned int m_seed;
	unsigned int gen()
	{
		m_seed = 1664525L*m_seed + 1013904223L;
		return m_seed;
	}
public:
	void init( int seed )
	{
		m_seed = seed;
	}
	float frand( void )
	{
		return (float)gen() * 2.318306447e-10f;
	}

	int irand( int iMax )
	{
		return ((int)( (iMax+1)*frand() )%iMax);
	}
};
#endif
