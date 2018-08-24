#ifndef _ANTSYSTEM_INC_
#define _ANTSYSTEM_INC_
#include "ant_xp.h"
#include "antsystemhelp.h"
#include "timers.h"
#include "tsp.h"
#include "platform.h"

class AntSystem
{
	friend class Ant;
public:	
	Ant *m_pAnts;
	int m_nAnts;

	TSP *m_pTSP;

	float m_shortestDist;
	int *m_shortestTour;

	float meanEntropy;
	float meanLambda;

	float alpha;
	float beta;
	float rho;
	float mmasConst; // pherMin/pherMax for MMAS
	Timers *timers;

	// iteration count at which stagnation occured
	int iStagnation;
	// stagnation metrics at point of stagnation
	float stagnationEntropy;
	float stagnationLambda;
	float stagnationTourLength;
	void DoTours( void );
	void DepositFromTour( int *tour, float tourLength );
	void Deposit( void );
	void Evaporate( void );
	void Iterate( void );
	void CalcStagnationMetrics( void );

//protected:
	// pheromone information
	float (**m_pher);     // pheromone value
	float (**m_weights);  // edge weights
	float (**m_iDistSq);  // precalculated inverse square of edge distances
	float (**m_fNN);		// nearest neighbour flags

//public:
	void Init( int nAnts, TSP *tsp, int seed );
	void Clear( void );
	void Solve( int maxIterations, int maxStagnantIterations, bool continueStagnant = false );
	TSP *GetTSP( void ) { return m_pTSP; }
	int fallbackCount;
	int usingNNCount;

	// return results
	float GetShortestTourLength() { return m_shortestDist; }
	float GetFinalEntropy( void ) { return meanEntropy; }
	float GetFinalLambdaBranching( void ) { return meanLambda; }
	int GetStagnationIteration( void ) { return iStagnation; }
	float GetStagnationEntropy( void ) { return stagnationEntropy; }
	float GetStagnationLambdaBranching( void ) { return stagnationLambda; }
	float GetStagnationTourLength( void ) { return stagnationTourLength; }
	float GetPheromoneTime() { return (float)timers->GetTimer(1); }
	float GetTourTime() { return (float)timers->GetTimer(0); }

};

#endif
