#include "tsp.h"
#include "antsystem_xp.h"
#include "timers.h"
#include <stdio.h>

int main( int argc, char *argv[] )
{
#ifdef VTUNECAPTURE
	__itt_pause(); // pause capture so we aren't picking up all the setup
#endif
 if ( argc != 4 && argc != 5 && argc != 6 )
    {
      fprintf(stderr,"parameters: problemfile, number of iterations, nNeighbours, number of ants, (seed)\n");
    }

  TSP *tsp = new TSP();
  int nNeighbours;
  sscanf( argv[3], "%d", &nNeighbours );


  //creates edge matrix and NN list
  tsp->Init( argv[1], nNeighbours );
  
  int nIter;
  int nAnts;
  //determines the number of iterations
  sscanf( argv[2],"%d", &nIter );

  sscanf( argv[4],"%d", &nAnts );

  Timers *timers = new Timers(10);

  AntSystem *as = new AntSystem();

  int seed = 12345;

  //determines seed, if specified
  if ( argc == 6 )
	  sscanf( argv[5], "%d", &seed );

  as->Init( nAnts, tsp, seed );
  timers->StartTimer(0);
#ifdef VTUNECAPTURE
  __itt_resume();
#endif

  as->Solve( nIter, nIter );
  timers->StopTimer(0);
  printf("shortest tour: %f time %lf\n",as->GetShortestTourLength(),timers->GetTimer(0) );
}
