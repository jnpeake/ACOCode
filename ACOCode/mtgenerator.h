#ifndef _MTGENERATOR_INC_
#define _MTGENERATOR_INC_
#include <random>

class MTGenerator
{
private:
	std::mt19937 m_gen;
public:
	void init( int seed )
	{
		m_gen.seed( seed );
	}
	float frand( void )
	{
		return (float)m_gen() / (float)m_gen.max();
	}

	int irand( int iMax )
	{
		return ((int)( (iMax+1)*frand() )%iMax);
	}
};
#endif