#include "tsp.h"
#include "antsystem.h"
#include "timers.h"
#include <stdio.h>

int main( int argc, char *argv[] )
{
	// parameters:
	// instance file, num trials, iterations per trial, number of near neighbours, seed, output filename

	if ( argc != 7 )
	{
		fprintf(stderr, "need 6  arguments: instance file, num trials, num iterations, num near neighbours, seed, out file");
	}

	TSP *tsp = new TSP();
	int nNeighbours;
	sscanf( argv[4], "%d", &nNeighbours );
	tsp->Init( argv[1], nNeighbours );
	
	int nIter;
	sscanf( argv[3],"%d", &nIter );
	
	int nTrials;
	sscanf( argv[2], "%d", &nTrials );

	AntSystem *as = new AntSystem();
	int seed;
	sscanf( argv[5], "%d", &seed );
	
	as->Init( tsp->numVerts, tsp, seed );


	char *outputFileName = (char*)malloc( strlen( argv[6] ) + 1 );
	strcpy( outputFileName, argv[6] );

	FILE *fp;
	fp = fopen(outputFileName,"w");
	fclose(fp);

	for ( int i = 0; i < nTrials; i++ )
	{
		as->Solve( nIter, nIter );
		float t0, t1;
		float length;
		t0 = as->GetTourTime();
		t1 = as->GetPheromoneTime();
		length = as->GetShortestTourLength();

		printf("iteration %d: %e %e %e\n",i,length, t0, t1 ); fflush(stdout);
		fp = fopen( outputFileName, "a" );
		fprintf(fp, "%e %e %e\n",length, t0, t1 );
		fclose(fp);
	}

}
