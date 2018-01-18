#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <omp.h>
#include "ranluxgen.h"
#include "mtgenerator.h"
#include "qdgenerator.h"
//#include "seeds.h"
#include "vec2.h"


int ALPHA = 1;
int BETA = 2;
float EVAP = 0.5f;

#define RS_CROULETTE 0 // classic
#define RS_IROULETTE 1 // independent
#define RS_PROULETTE 2 // pairs

#define SEED 1234

double timers[5] = {0.0};
#define START_TIME( _i ) {timers[_i] -= microsecs();}
#define STOP_TIME( _i ) {timers[_i] += microsecs();}

int g_numThreads;

#ifndef _WIN32
#include <sys/time.h>
double microsecs( void )
{
	timeval t;
	gettimeofday( &t, NULL );
	return (double)(t.tv_sec*1e6 + t.tv_usec);
}
#else
#include <Windows.h>
double pcFreq = 0.0;
double microsecs( void )
{
	LARGE_INEGER count;
	QueryPerformanceCounter( &count );
	if ( pcFreq == 0.0 )
	{
		LARGE_INTEGER freq;
		QueryPerformanceFrequency( &freq );
		pcFreq = 1e6 / (double)freq.QuadPart;
	}
	return (double)count.QuadPart * pcFreq;
}
#endif

#define QD_RAND 0
inline unsigned long qdrand( void )
{
	// Numerical Recipes quick and dirty RNG
	static int seed = 0x12344321;
	seed = seed*1664525 + 1013904223;
	return seed;
}

float frand( void )
{
#if QD_RAND
	return (float)qdrand() * ( 1.0f / (float)0xFFFFFFFF ); // compiler should do the divide
#else
	return (float)rand()/(float)RAND_MAX;
#endif
}

int iRand( int iMax )
{
	return (int) (iMax*frand());
}

int g_iSortTemp;
void *g_tsp;
int distComp( const void *v1, const void *v2 );

class TSP
{
public:
	typedef struct
	{
		float dist;
		float pher;
		float cval;
		float idist;
	} EDGE;

	EDGE **edges;
	int numVerts;
	V2 *verts;

	int **m_nearList;
	int m_nearListSize;

	void Init( const char *fileName )
	{
		int i, j, k;

		numVerts= 0;

		FILE *fp = fopen( fileName, "r");

		if ( fp == NULL )
		{
			fprintf(stderr, "could not open %s for read\n", fileName );
			return;
		}

		char *line = (char*)malloc( 1024 );

		// read the header 
		int stillHeader = 1;
		while ( stillHeader )
		{
			fscanf( fp, "%s", line );
			if ( strcmp(line,"NODE_COORD_SECTION") == 0 )
				stillHeader = 0;
			if ( strcmp(line,"DIMENSION") == 0 )
			{
				fscanf( fp, "%s", line );
				fscanf( fp, "%s", line );
				sscanf(line,"%d", &numVerts );
			}
			else if ( strcmp(line,"DIMENSION:") == 0 )
			{
				fscanf( fp, "%s", line );
				sscanf(line,"%d", &numVerts );
			}

		}
		//printf("read the header, number of vertices is %d\n",numVerts );
		// read in the vertices
		verts = (V2*)malloc( numVerts * sizeof(V2) );
		int iVert = 0;
		for ( i = 0; i < numVerts; i++ )
		{
			float x, y;
			fscanf( fp, "%d %e %e", &j, &x, &y );
			V2 v = V2( x, y );
			// check for duplicates
			int found = 0;
			for ( int j = 0; j < iVert; j++ )
			{
				int dist = (int)((v - verts[j]).magnitude()+0.5f);
				if ( dist == 0 )
				{
					found = 1;
					break;
				}
			}
			if ( !found )
				verts[iVert++] = V2( x, y );
		}
		fclose(fp);
		numVerts = iVert;
#if 0
		fp = fopen("vertsread.dat", "w");
		for ( i = 0; i < numVerts; i++ )
		{
			fprintf(fp,"%f %f\n",verts[i].x,verts[i].y);
		}
		fclose(fp);
#endif
		// construct the edges
		edges = (EDGE**)malloc( numVerts * sizeof( EDGE* ) );
		for ( i = 0; i < numVerts; i++ )
			edges[i] = (EDGE*)malloc( numVerts * sizeof( EDGE ) );

		for ( i = 0; i < numVerts; i++ )
		{
			for ( j = 0; j < numVerts; j++ )
			{
				edges[i][j].dist = (int)((verts[i] - verts[j]).magnitude()+0.5f);
			}
		}
		// don't need the vertices any more
//		free( verts );
		// nearest neighbour lists
#if 0
		m_nearListSize = numVerts > 20 ? 20 : numVerts;
		m_nearList = (int**)malloc( numVerts * sizeof(int*));
		int *distBuffer =  (int*)malloc( numVerts - 1 );

		for ( i = 0; i < numVerts; i++ )
		{
			m_nearList[i] = (int*)malloc( m_nearListSize * sizeof(int) );
			k = 0;
			for ( j = 0; j < numVerts; j++ )
			{
				if ( j != i )
					distBuffer[k++] = j;
			}
			g_iSortTemp = i;
			g_tsp = this;
			qsort( distBuffer, numVerts-1, sizeof(int), &distComp ); // ugly
			for ( j = 0; j < m_nearListSize; j++ )
			{
				m_nearList[i][j] = distBuffer[j];
			}
		}
#endif
	}
	
	float GetPher( int i, int j )
	{
		return edges[i][j].pher;
	}

	float GetDist( int i, int j )
	{
		return edges[i][j].dist;
	}

	void Evaporate( void )
	{
		int i, j;
		for ( i = 0; i < numVerts; i++ )
		{
			for ( j = 0; j < numVerts; j++ )
			{
				edges[i][j].pher *= EVAP;
			}
		}
	}
};

int distComp( const void *v1, const void *v2 )
{
	int i1, i2;
	i1 = *(int*)v1;
	i2 = *(int*)v2;
	float d0, d1;
	d0 = ((TSP*)g_tsp)->edges[g_iSortTemp][i1].dist;
	d1 = ((TSP*)g_tsp)->edges[g_iSortTemp][i2].dist;
	if ( d0 < d1 ) return -1;
	if ( d0 > d1 ) return 1;
	return 0;
}

int g_debugTour = 0;


class ANT
{
public:
	float tourDist;
	float *rouletteVals;
	int *tour;
	int *remaining;
	int nTour;
	TSP *m_tsp;

	int *rIndices;
	//RanluxGenerator m_gen;
	//MTGenerator m_gen;
	QDGenerator m_gen;

	void Init( TSP *tsp )
	{
		m_tsp = tsp;
		remaining = (int*)malloc( tsp->numVerts * sizeof(int) );
		tour = (int*)malloc( tsp->numVerts * sizeof(int) );
		rouletteVals = (float*)malloc( tsp->numVerts * sizeof(float) );
		nTour = 0;

		rIndices = (int*)malloc( tsp->numVerts * sizeof(int) );
	}

	void PrintTour( void )
	{
		int i;
		if ( m_tsp->numVerts < 7 )
		{
			printf("[ ");
			for ( i = 0; i < m_tsp->numVerts; i++ )
				printf("%d ",tour[i] );
			printf("]\n");
		}
		else
		{
			printf("[ %d %d %d ... %d %d %d ]\n",tour[0],tour[1],tour[2],
				   tour[m_tsp->numVerts-3],tour[m_tsp->numVerts-2],tour[m_tsp->numVerts-1]);
		}
	}

	void ConstructTour( TSP *tsp, int rouletteScheme = RS_CROULETTE )
	{
		int i, j;
		tourDist = 0.0f;
		nTour = 0;
		memset( tour, 0, tsp->numVerts * sizeof( int ) );

		// pick a random starting point
		tour[0] = m_gen.irand( tsp->numVerts ) % tsp->numVerts;

		nTour = 1;

		j = 1;
		for ( i = 0; i < tsp->numVerts; i++ )
		{
			if ( i != tour[0] )
				tour[j++] = i;
		}
		while ( nTour < tsp->numVerts-1 )
		{
			// set up the roulette selection of the next city
			int iNext;
			if ( rouletteScheme == RS_IROULETTE )
			{
				float fMax = -1e20f;
				int iMax = -1;
				for ( iNext = 0; iNext < tsp->numVerts - nTour; iNext++ )
				{
					float thisVal = tsp->edges[tour[nTour-1]][tour[iNext+nTour]].cval * m_gen.frand();
					if ( thisVal > fMax )
					{
						fMax = thisVal;
						iMax = iNext;
					}
				}
				int swap = tour[nTour];
				tour[nTour] = tour[iMax+nTour];
				tour[iMax+nTour] = swap;
				// add to the distance
				tourDist += tsp->GetDist( tour[nTour], tour[nTour-1] );
				nTour++;
			}
			else if ( rouletteScheme == RS_PROULETTE )
			{
				int nLeft = tsp->numVerts - nTour;
				for ( i = 0; i < nLeft; i++ )
				{
					rIndices[i] = i;
					rouletteVals[i] = tsp->edges[tour[nTour-1]][tour[i+nTour]].cval;
				}

				int istep = 1;
				while ( istep < nLeft )
				{
					float randomVal = m_gen.frand();
					for ( int i = 0; i < nLeft; i += 2*istep )
					{
						float w0, w1;
						w0 = rouletteVals[i];
						if ( i+istep < nLeft )
							w1 = rouletteVals[i+istep];
						else 
							w1 = 0.0f;
						if ( randomVal*(w0+w1) > w0 )
							rIndices[i] = rIndices[i+istep];
						rouletteVals[i] = w0+w1;
					}
					istep *= 2;
				}
				int iMax = rIndices[0];
				int swap = tour[nTour];
				tour[nTour] = tour[iMax+nTour];
				tour[iMax+nTour] = swap;
				// add to the distance
				tourDist += tsp->GetDist( tour[nTour], tour[nTour-1] );
				nTour++;
			}
			else
			{
				float rouletteTot = 0.0f;
				for ( iNext = 0; iNext < tsp->numVerts - nTour; iNext++ )
				{
					float thisVal = tsp->edges[tour[nTour-1]][tour[iNext+nTour]].cval;
					rouletteVals[iNext] = thisVal;
					rouletteTot += thisVal;
				}
				// pick a city
				float rouletteVal = m_gen.frand() * rouletteTot;
				rouletteTot = 0.0f;
				for ( iNext = 0; iNext < tsp->numVerts - nTour; iNext++ )
				{
					rouletteTot += rouletteVals[iNext];
					if ( rouletteTot >= rouletteVal )
					{
						// choose this one - swap the chosen one into the right slot in the tour
						int swap = tour[nTour];
						tour[nTour] = tour[iNext+nTour];
						tour[iNext+nTour] = swap;
						// add to the distance
						tourDist += tsp->GetDist( tour[nTour], tour[nTour-1] );
						nTour++;
						break;
					}
				}
			}
		}
		// last one should be already in place
		tourDist += tsp->GetDist( tour[nTour], tour[nTour-1] );
		nTour++; 
		tourDist += tsp->GetDist( tour[0], tour[nTour-1] );

	 }
} ;

class AS
{
public:

	ANT *m_pAnts;
	int m_nAnts;
	TSP *m_pTSP;
	int m_converged;
	float m_shortestDist;
	int *m_shortestTour;
	int m_rouletteScheme;
	float meanEntropy;
	float meanLambda;

	void Init( int nAnts, TSP *tsp )
	{
		m_nAnts = nAnts;
		m_pTSP = tsp;
		m_converged = 0;
		m_shortestDist = 1e20f;
		m_rouletteScheme = RS_CROULETTE;

		m_pAnts = (ANT*)malloc( m_nAnts * sizeof( ANT ) );
		int i;

		for ( i = 0; i < m_nAnts; i++ )
		{
			m_pAnts[i].Init( tsp );
		}

		m_shortestTour = (int*)malloc( m_pTSP->numVerts * sizeof(int) );

		Clear();

	}

	void Clear( void )
	{
		m_shortestDist = 1e20f;
		m_converged = 0;

		float aDist = 0.0f;

		int i, j;

		for ( i = 0; i < m_pTSP->numVerts; i++ )
		{
			aDist += m_pTSP->GetDist(i, (i+1)%m_pTSP->numVerts );;
		}
		float val = (float)m_nAnts / aDist;
		for ( i = 0; i < m_pTSP->numVerts; i++ )
		{
			for ( j = 0; j < m_pTSP->numVerts; j++ )
			{
				m_pTSP->edges[i][j].pher = val;
				m_pTSP->edges[i][j].idist = 1.0f/m_pTSP->edges[i][j].dist;
				m_pTSP->edges[i][j].cval = pow(m_pTSP->edges[i][j].pher, ALPHA)/pow(m_pTSP->edges[i][j].dist,BETA);
			}
		}

	}
	void DoTours( void )
	{
		int i;
		int iNewShortest = -1;

//#pragma omp parallel for num_threads(4)
		for ( i = 0; i < m_nAnts; i++ )
		{
			//g_debugTour = (i==0);
			m_pAnts[i].ConstructTour( m_pTSP, m_rouletteScheme );
		}

		for ( i = 0; i < m_nAnts; i++ )
		{
			if ( m_pAnts[i].tourDist < m_shortestDist )
			{
				m_shortestDist = m_pAnts[i].tourDist;
				iNewShortest = i;
			}
		}
		// copy new shortest tour, if any
		if ( iNewShortest != -1 )
		{
			memcpy( m_shortestTour, m_pAnts[iNewShortest].tour, m_pTSP->numVerts * sizeof(int) );
		}
	}

	float CalcConvergence()
	{
		float totPher = 0.0f;
		float tourPher = 0.0f;
		int i, j;
		for ( i = 0; i < m_pTSP->numVerts-1; i++ )
		{
			for ( j = i+1; j < m_pTSP->numVerts; j++ )
			{
				totPher += m_pTSP->edges[i][j].pher;
			}
		}
		for ( j = 0; j < m_pTSP->numVerts; j++ )
		{
			tourPher += m_pTSP->edges[m_shortestTour[j]][m_shortestTour[(j+1)%m_pTSP->numVerts]].pher;
		}
		return 1.0f - tourPher/totPher;
	}
	
	int Converged( void )
	{
		return m_converged;
	}

	void Deposit( void )
	{
		int i, j;

		for ( i = 0; i < m_nAnts; i++ )
		{
			float deltaPher = 1.0f / m_pAnts[i].tourDist;

			for ( j = 0; j < m_pTSP->numVerts; j++ )
			{
				m_pTSP->edges[m_pAnts[i].tour[j]][m_pAnts[i].tour[(j+1)%m_pTSP->numVerts]].pher += deltaPher;
				m_pTSP->edges[m_pAnts[i].tour[(j+1)%m_pTSP->numVerts]][m_pAnts[i].tour[j]].pher += deltaPher;
			}
		}
		// calculate edge probabilities, entropy and lambda branching
		meanEntropy = 0.0f;
		meanLambda = 0.0f;
		for ( i = 0; i < m_pTSP->numVerts; i++ )
		{
			float entropy = 0.0f;
			float lambda = 0.0f;
			float pMin, pMax;
			pMin = 1e20f;
			pMax = 0.0f;
			for ( j = 0; j < m_pTSP->numVerts; j++ )
			{
				float p;
				if ( m_pTSP->edges[i][j].dist > 0.0f )
					p = m_pTSP->edges[i][j].pher / (m_pTSP->edges[i][j].dist * m_pTSP->edges[i][j].dist );
				else
					p = 0.0f;

				m_pTSP->edges[i][j].cval = m_pTSP->edges[i][j].pher / (m_pTSP->edges[i][j].dist * m_pTSP->edges[i][j].dist );;

				if ( m_pTSP->edges[i][j].pher < pMin )
					pMin = m_pTSP->edges[i][j].pher;
				if ( m_pTSP->edges[i][j].pher > pMax )
					pMax = m_pTSP->edges[i][j].pher;
				if ( p > 0.0f )
					entropy -= p*log(p);
			}
			float pThresh = pMin + 0.05f*(pMax - pMin);
			for ( j = 0; j < m_pTSP->numVerts; j++ )
			{
				float p = m_pTSP->edges[i][j].pher;
				if ( p > pThresh &&  i != j)
					lambda += 1.0f;
				if ( p > pThresh )
				{
					lambda += 1.0f;			
				}
			}
			meanEntropy += entropy;
			meanLambda += lambda;
		}
		meanEntropy /= (float)m_pTSP->numVerts;
		meanLambda /= (float)m_pTSP->numVerts;
	}

	void SetRouletteScheme( int nr )
	{
		m_rouletteScheme = nr;
	}

	void Iterate( void )
	{
		DoTours();
		m_pTSP->Evaporate();
		Deposit();
	}
};

#define NITERS 1000
#define NTRIALS 128

#define USE_MPI
#ifdef USE_MPI

#include <mpi.h>
#endif

typedef struct
{
	float shortestTour;
	float meanEntropy;
	float meanLambda;
	int nIterations;
} TRIAL_RESULT;

int main( int argc, char *argv[] )
{
	int mpiWorldSize;
	int mpiWorldRank;
	// Initialize the MPI environment
#ifdef USE_MPI
	MPI_Init(&argc, &argv);
#else
	mpiWorldSize = 1;
	mpiWorldRank = 0;
#endif
    char fileName[64];
	char problemFile[64];

	if ( argc == 6 )
	{
		sscanf(argv[1],"%d",&ALPHA);
		sscanf(argv[2],"%d",&BETA);
		sscanf(argv[3],"%f",&EVAP);
		strcpy(fileName,argv[4]);
		strcpy(problemFile,argv[5]);
	}
#ifdef USE_MPI	
	// Get the number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &mpiWorldSize);
	
	// Get the rank of the process
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiWorldRank);
#endif
	TSP tsp;
	AS as;

	tsp.Init(problemFile);
	as.Init( tsp.numVerts, &tsp );

	// seed the ants' RNGs
	int iRandInt = (64 + mpiWorldRank) * as.m_nAnts;
	int i;

	unsigned int *randomInts;
	FILE *fp = fopen("randoms.dat", "rb" );

	fseek( fp, 0, SEEK_END );
	int randSize = ftell( fp );
	fseek( fp, 0, SEEK_SET );

	randomInts = (unsigned int*)malloc( randSize );
	fread( randomInts, randSize, 1, fp );
	fclose( fp );


	for ( i = 0; i < as.m_nAnts; i++ )
	{
		as.m_pAnts[i].m_gen.init( randomInts[iRandInt++] );
	}

	free( randomInts );
	
	int numTrials = NTRIALS/mpiWorldSize;

	TRIAL_RESULT *trialResults0 = (TRIAL_RESULT*)malloc( numTrials * sizeof(TRIAL_RESULT) );
	TRIAL_RESULT *trialResults1 =  (TRIAL_RESULT*)malloc( numTrials * sizeof(TRIAL_RESULT) );

	int dumpAtEnd = 0;

// allocate space for the answers
	TRIAL_RESULT *allRes0 = (TRIAL_RESULT*)malloc( NTRIALS * sizeof (TRIAL_RESULT) );
	TRIAL_RESULT *allRes1 = (TRIAL_RESULT*)malloc( NTRIALS * sizeof (TRIAL_RESULT) );

	for ( int iTrial = 0; iTrial < numTrials; iTrial++ )
	{
		float lengths[2];
		int step;
		double t0 = microsecs();
		as.Clear();
		as.SetRouletteScheme( RS_IROULETTE );
		int sinceChange = 0;
		float shortest = 1e20f;

		for ( step = 0; step < NITERS && sinceChange < 200; step++ )
		{
			as.Iterate();
			if ( as.m_shortestDist < shortest )
			{
				shortest = as.m_shortestDist;
				sinceChange = 0;
			}
			else
			{
				sinceChange++;
			}
		}

		trialResults0[iTrial].shortestTour = as.m_shortestDist;
		trialResults0[iTrial].meanLambda = as.meanLambda;
		trialResults0[iTrial].meanEntropy = as.meanEntropy;
		trialResults0[iTrial].nIterations = step;

		as.Clear();
		as.SetRouletteScheme( RS_CROULETTE );
			 
		sinceChange = 0;
		shortest = 1e20f;

		for ( step = 0; step < NITERS && sinceChange < 200; step++ )
		{
			as.Iterate();
			if ( as.m_shortestDist < shortest )
			{
				shortest = as.m_shortestDist;
				sinceChange = 0;
			}
			else
			{
				sinceChange++;
			}
		}

		trialResults1[iTrial].shortestTour = as.m_shortestDist;
		trialResults1[iTrial].meanLambda = as.meanLambda;
		trialResults1[iTrial].meanEntropy = as.meanEntropy;
		trialResults1[iTrial].nIterations = step;

		double t1 = microsecs();

		if (mpiWorldRank == 0 )
			printf("Trial %d: %e seconds %d %d\n",iTrial,(t1-t0)*1e-6,trialResults0[iTrial].nIterations,trialResults1[iTrial].nIterations );

		if ( !dumpAtEnd || iTrial == numTrials-1 )
		{
			// gather the data so far and dump to a file
			if ( mpiWorldRank == 0 )
			{
				// copy this rank's results in
				memcpy( allRes0, trialResults0, (iTrial+1) * sizeof(TRIAL_RESULT) );
				memcpy( allRes1, trialResults1, (iTrial+1) * sizeof(TRIAL_RESULT) );
				// receive the rest
#ifdef USE_MPI
				for ( int iRank = 1; iRank < mpiWorldSize; iRank++ )
				{
					MPI_Recv( allRes0 + (iTrial+1)*iRank, (iTrial+1)* sizeof(TRIAL_RESULT) , MPI_FLOAT, iRank, 0, MPI_COMM_WORLD, NULL );
					MPI_Recv( allRes1 + (iTrial+1)*iRank, (iTrial+1)* sizeof(TRIAL_RESULT) , MPI_FLOAT, iRank, 1, MPI_COMM_WORLD, NULL );
				}
#endif
				// write the file
				char thisFileName[64];
				sprintf(thisFileName,"%s%3.3d.dat",fileName,iTrial );
				FILE *fp = fopen(thisFileName, "w");
				for ( int i = 0; i < (iTrial+1)*mpiWorldSize; i++ )
					fprintf(fp, "%e %e %e %d %e %e %e %d\n",allRes0[i].shortestTour,allRes0[i].meanEntropy,allRes0[i].meanLambda,
							allRes0[i].nIterations, allRes1[i].shortestTour,allRes1[i].meanEntropy,allRes1[i].meanLambda, allRes1[i].nIterations );
				fclose(fp);
			}
			else
			{
#ifdef USE_MPI
				// send to rank 0
				MPI_Send( trialResults0, (iTrial+1)* sizeof(TRIAL_RESULT), MPI_FLOAT, 0, 0, MPI_COMM_WORLD );
				MPI_Send( trialResults1, (iTrial+1)* sizeof(TRIAL_RESULT), MPI_FLOAT, 0, 1, MPI_COMM_WORLD );
#endif
			}
		}
		MPI_Barrier( MPI_COMM_WORLD );
	}
#ifdef USE_MPI
	MPI_Finalize();
#endif

	free( allRes0 );
	free( allRes1 );
return 0;


}
