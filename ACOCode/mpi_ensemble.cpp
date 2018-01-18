#include "tsp.h"
#include "antsystem.h"
#include "timers.h"
#include <stdio.h>
#include <mpi.h>

typedef struct
{
	float tourLength;
	float finalEntropy;
	float finalLambda;
	int nTrialsToStagnation;
	float stagnationLength;
	float stagnationEntropy;
	float stagnationLambda;
	float _pad;
} TRIAL_RESULT;

int main( int argc, char *argv[] )
{
	int mpiWorldSize;
	int mpiWorldRank;

	// Initialize the MPI environment
	MPI_Init( &argc, &argv );
	// Get the number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &mpiWorldSize);
	// Get the rank of the process
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiWorldRank);

    char outputFile[64];
	char problemFile[64];
	int numTrials, maxIter, stagIter;

	if ( argc == 6 )
	{
		strcpy( problemFile, argv[1] );
		strcpy( outputFile, argv[2] );
		sscanf(argv[3],"%d",&numTrials);
		sscanf(argv[4],"%d",&maxIter);
		sscanf(argv[5],"%d",&stagIter);
	}
	else
	{
		fprintf(stderr, "need 6 arguments\n");
	}

	// load the problem
	TSP *tsp = new TSP();
	tsp->Init( problemFile );

	// intialise the ant system
	AntSystem *as = new AntSystem();
	// number of ants = number of vertices, seed based on MPI rank
	as->Init( tsp->numVerts, tsp, 1000 + mpiWorldRank );

	int nTrialsPerRank = numTrials / mpiWorldSize;
	TRIAL_RESULT *pRes = (TRIAL_RESULT*)malloc( nTrialsPerRank * 2 * sizeof( TRIAL_RESULT ) );

	// do the roulette trials
	as->SetRouletteScheme( RS_ROULETTE );
	for ( int iTrial = 0; iTrial < nTrialsPerRank; iTrial++ )
	{
		if ( mpiWorldRank == 0 )
			printf("problem: %s roulette iter %d\n",problemFile, iTrial);
		as->Solve( maxIter, stagIter, true );
		// get the results of the trial
		pRes[iTrial].tourLength = as->GetShortestTourLength();
		pRes[iTrial].finalEntropy = as->GetFinalEntropy();
		pRes[iTrial].finalLambda = as->GetFinalLambdaBranching();
		int iStag = as->GetStagnationIteration();
		if ( iStag == -1 )
		{
			// didn't reach stagnation
			pRes[iTrial].nTrialsToStagnation = maxIter;
			pRes[iTrial].stagnationLength = as->GetShortestTourLength();
			pRes[iTrial].stagnationEntropy = as->GetFinalEntropy();
			pRes[iTrial].stagnationLambda = as->GetFinalLambdaBranching();
		}
		else
		{
			pRes[iTrial].nTrialsToStagnation = iStag;
			pRes[iTrial].stagnationLength = as->GetStagnationTourLength();
			pRes[iTrial].stagnationEntropy = as->GetStagnationEntropy();
			pRes[iTrial].stagnationLambda = as->GetStagnationLambdaBranching();
		}
	}
	// do the iRoulette trials
	as->SetRouletteScheme( RS_IROULETTE );
	for ( int iTrial = nTrialsPerRank; iTrial < 2*nTrialsPerRank; iTrial++ )
	{
		if ( mpiWorldRank == 0 )
			printf("problem: %s iroulette iter %d\n",problemFile,iTrial);
		as->Solve( maxIter, stagIter, true );
		// get the results of the trial
		pRes[iTrial].tourLength = as->GetShortestTourLength();
		pRes[iTrial].finalEntropy = as->GetFinalEntropy();
		pRes[iTrial].finalLambda = as->GetFinalLambdaBranching();
		int iStag = as->GetStagnationIteration();
		if ( iStag == -1 )
		{
			// didn't reach stagnation
			pRes[iTrial].nTrialsToStagnation = maxIter;
			pRes[iTrial].stagnationLength = as->GetShortestTourLength();
			pRes[iTrial].stagnationEntropy = as->GetFinalEntropy();
			pRes[iTrial].stagnationLambda = as->GetFinalLambdaBranching();
		}
		else
		{
			pRes[iTrial].nTrialsToStagnation = iStag;
			pRes[iTrial].stagnationLength = as->GetStagnationTourLength();
			pRes[iTrial].stagnationEntropy = as->GetStagnationEntropy();
			pRes[iTrial].stagnationLambda = as->GetStagnationLambdaBranching();
		}
	}

	// dump the results...
	if ( mpiWorldRank == 0 )
	{
		// open a file, dump the results of this rank
		FILE *fp = fopen( outputFile, "w" );

		for ( int i = 0; i < 2*nTrialsPerRank; i++ )
		{
			fprintf(fp, "%d %e %e %e %e %e %e\n",
					pRes[i].nTrialsToStagnation,
					pRes[i].tourLength,
					pRes[i].finalEntropy, 
					pRes[i].finalLambda,
					pRes[i].stagnationLength,
					pRes[i].stagnationEntropy,
					pRes[i].stagnationLambda ); 
		}

		// now receive results from the other ranks in turn
		for ( int iRank = 1; iRank < mpiWorldSize; iRank++ )
		{
			MPI_Recv( pRes, 2*nTrialsPerRank*sizeof(TRIAL_RESULT)/4, MPI_FLOAT, iRank, 0, MPI_COMM_WORLD, NULL );

			for ( int i = 0; i < 2*nTrialsPerRank; i++ )
			{
				fprintf(fp, "%d %e %e %e %e %e %e\n",
						pRes[i].nTrialsToStagnation,
						pRes[i].tourLength,
						pRes[i].finalEntropy, 
						pRes[i].finalLambda,
						pRes[i].stagnationLength,
						pRes[i].stagnationEntropy,
						pRes[i].stagnationLambda ); 
			}
		}
		fclose(fp);
	}
	else
	{
		// send the results to rank 0
		MPI_Send( pRes, 2*nTrialsPerRank*sizeof(TRIAL_RESULT)/4, MPI_FLOAT, 0, 0, MPI_COMM_WORLD );
	}

	MPI_Finalize();
	
	return 0;
}
