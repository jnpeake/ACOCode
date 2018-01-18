#ifndef _TIMERS_INC_
#define _TIMERS_INC_

#ifdef _WIN32
// don't need this on windows for now - it just needs to build
class Timers 
{
	double microsecs( void )
	{
		return 0.0;
	}
public:
	Timers( int nTimers )
	{
	}
	void Clear( void )
	{
	}
	void StartTimer( int timer )
	{
	}
	void StopTimer( int timer )
	{
	}
	double GetTimer( int timer )
	{
		return 0.0;
	}
};

#else
#include <sys/time.h>

class Timers
{
	// cumulative timer class for profiling
	int numTimers;
	double *timers;

	double microsecs( void )
	{
		timeval t;
		gettimeofday( &t, 0 );
		return (double)(t.tv_sec*1e6 + t.tv_usec);
	}

public:
	Timers( int nTimers )
	{
		numTimers = nTimers;
		timers = new double[numTimers];
	}
	void Clear( void )
	{
		for ( int i = 0; i < numTimers; i++ )
			timers[i] = 0.0;
	}

	void StartTimer( int timer )
	{
		timers[timer] -= microsecs();
	}
	void StopTimer( int timer )
	{
		timers[timer] += microsecs();
	}
	double GetTimer( int timer )
	{
		return timers[timer] * 1e-6;
	}
};
#endif
#endif
