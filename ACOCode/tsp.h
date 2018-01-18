#ifndef _TSP_INC_
#define _TSP_INC_
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

typedef struct
{
	float dist;
	int index;
} distSort;

typedef enum 
{
	EDGE_WEIGHT_EUC2D,
	EDGE_WEIGHT_CEIL2D,
	EDGE_WEIGHT_GEO,
	EDGE_WEIGHT_ATT,
	EDGE_WEIGHT_UNSUPPORTED
} EdgeType;

static int nnComp( const void *p0, const void *p1 )
{
	if ( ((distSort*)p0)->dist < ((distSort*)p1)->dist )
		return -1;
	else if ( ((distSort*)p0)->dist > ((distSort*)p1)->dist )
		return 1;
	else
		return 0;
}

class TSP
{
public:
	float **edgeDist;
	int numVerts;
	int **nnList;
	int numNN;

	int (*distanceFunc)(float, float, float, float);

	void FillNNList( int iList )
	{
		distSort *tempList = (distSort*)malloc( (numVerts-1) * sizeof( distSort ) );
		int count = 0;
		for ( int i = 0; i < numVerts; i++ )
		{
			if ( i != iList )
			{
				tempList[count].index = i;
				tempList[count].dist = edgeDist[i][iList];
				count++;
			}
		}
		qsort( tempList, numVerts-1, sizeof(distSort), nnComp );
		for ( int i = 0; i < numNN; i++ )
		{
			nnList[iList][i] = tempList[i].index;
		}
		free( tempList );
	}

	static int EucDistance( float x0, float y0, float x1, float y1 )
	{
		float d = (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1);
		d = sqrt(d);
		return (int)(d+0.5f);
	}

	static int CeilDistance(float x0, float y0, float x1, float y1)
	{
		float d = (x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1);
		d = sqrt(d);
		return (int)ceil(d);
	}

	static int AttDistance(float x0, float y0, float x1, float y1)
	{
		float d = (x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1);
		d = sqrt(d * 0.1f);
		int td = ((int)(d+0.5f));

		if (td < d)
			return td + 1;
		else
			return td;
	}

	static int GeoDistance(float x0, float y0, float x1, float y1)
	{
		float lat0, long0, lat1, long1;
		int deg;
		float min;

		float PI = 3.141592f;
		
		deg = (int)x0;
		min = x0 - deg;
		lat0 = PI*(deg + 5.0f*min / 3.0f) / 180.0f;
		deg = (int)y0;
		min = y0 - deg;
		long0 = PI*(deg + 5.0f*min / 3.0f) / 180.0f;

		deg = (int)x1;
		min = x1 - deg;
		lat1 = PI*(deg + 5.0f*min / 3.0f) / 180.0f;
		deg = (int)y1;
		min = y1 - deg;
		long1 = PI*(deg + 5.0f*min / 3.0f) / 180.0f;

		float rad = 6378.388f;
		float q1 = cos(long0 - long1);
		float q2 = cos(lat0 - lat1);
		float q3 = cos(lat0 + lat1);
		int dist = (int)(rad * acos(0.5f*((1.0f + q1)*q2 - (1.0f - q1)*q3)) + 1.0f);
		return dist;
	}


	void Init( const char *fileName, int nNearNeighbours )
	{
		int i, j, k;

		float *vertX;
		float *vertY;

		numVerts= 0;

		FILE *fp = fopen( fileName, "r");

		if ( fp == NULL )
		{
			fprintf(stderr, "could not open %s for read\n", fileName );
			perror("Error");
			return;
		}

		char *line = (char*)malloc( 1024 );
		// read the header 
		int stillHeader = 1;
		EdgeType edgeType = EDGE_WEIGHT_UNSUPPORTED;

		while ( stillHeader )
		{
			fscanf( fp, "%s", line );
			if ( strcmp(line,"NODE_COORD_SECTION") == 0 )
				stillHeader = 0;
			else if ( strcmp(line,"DIMENSION") == 0 || strcmp(line,"DIMENSION:") == 0 ) 
			{
				if (strcmp(line, "DIMENSION:") != 0) // may or may not be a space between word and colon, read the colon if so
					fscanf( fp, "%s", line );
				fscanf( fp, "%s", line );
				sscanf(line,"%d", &numVerts );
			}
			else if (strcmp(line, "EDGE_WEIGHT_TYPE") == 0 || strcmp(line,"EDGE_WEIGHT_TYPE:")==0)
			{
				if (strcmp(line, "EDGE_WEIGHT_TYPE:") != 0) // may or may not be a space between word and colon, read the colon if so
					fscanf(fp, "%s", line);
				fscanf(fp, "%s", line);
				if (strcmp(line, "EUC_2D") == 0)
					edgeType = EDGE_WEIGHT_EUC2D;
				else if (strcmp(line, "CEIL_2D") == 0)
					edgeType = EDGE_WEIGHT_CEIL2D;
				else if (strcmp(line, "GEO") == 0)
					edgeType = EDGE_WEIGHT_GEO;
				else if (strcmp(line, "ATT") == 0)
					edgeType = EDGE_WEIGHT_ATT;
				else
					stillHeader = 0; // unsupported edge type, quit reading
			}
		}
		switch (edgeType)
		{
		case EDGE_WEIGHT_UNSUPPORTED:
			fprintf(stderr, "unsupported edge weight type in tsplib file\n");
			return;
		case EDGE_WEIGHT_EUC2D:
			distanceFunc = &EucDistance;
			break;
		case EDGE_WEIGHT_CEIL2D:
			distanceFunc = &CeilDistance;
			break;
		case EDGE_WEIGHT_ATT:
			distanceFunc = &AttDistance;
			break;
		case EDGE_WEIGHT_GEO:
			distanceFunc = &GeoDistance;
			break;
		}

		// read in the vertices
		vertX = (float*)malloc( numVerts * sizeof(float) );
		vertY = (float*)malloc( numVerts * sizeof(float) );
		int iVert = 0;
		for ( i = 0; i < numVerts; i++ )
		{
			float x, y;
			fscanf( fp, "%d %e %e", &j, &x, &y );

			int found = 0;
			for ( int j = 0; j < iVert; j++ )
			{
				int dist = distanceFunc( x, y, vertX[i], vertY[j] );
				if ( dist == 0 )
				{
					found = 1;
					break;
				}
			}
			if ( !found )
			{
				vertX[iVert] =  x;
				vertY[iVert] =  y;
				iVert++;
			}
		}
		fclose(fp);
		numVerts = iVert;

		// construct the edges
		edgeDist = (float**)malloc( numVerts * sizeof( float* ) );
		for ( i = 0; i < numVerts; i++ )
			edgeDist[i] = (float*)malloc( numVerts * sizeof( float ) );

		for ( i = 0; i < numVerts; i++ )
		{
			for ( j = 0; j < numVerts; j++ )
			{
				if (i == j)
					edgeDist[i][j] = 0;
				else
					edgeDist[i][j] = distanceFunc( vertX[i], vertY[i], vertX[j], vertY[j] );
			}
		}
		// don't need the vertex positions any more
		free( vertX );
		free( vertY );
		// initialize the nearest neighbour lists
		numNN = nNearNeighbours;

		nnList = (int **)malloc( numVerts * sizeof(int*) );

		for ( i = 0; i < numVerts; i++ )
		{
			nnList[i] = (int*)malloc( numNN * sizeof( int ) );
			FillNNList( i );
		}
	}

};
		
#endif
