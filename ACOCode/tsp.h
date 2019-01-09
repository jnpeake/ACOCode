#ifndef _TSP_INC_
#define _TSP_INC_
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <omp.h>
#include "platform.h"

typedef struct
{
	float dist;
	int index;
} distSort;

typedef struct
{
	int nnMask; //mask with bits set for vector lanes which contains NN
	int vectIndex = -1; //index of pheromone matrix vector which contains NNs
} nearestNeighbour;

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
	float *vertX;
	float *vertY;
	float nnDist;
	int numVerts;
	int **nnList;
	int **nnIndexes;
	int numNN;
	int nnHist[32] = {};
	int matrixSize;
	nearestNeighbour **neighbourVectors;
	nearestNeighbour *newNN;
	int (*distanceFunc)(float, float, float, float);
	float (*distanceSquaredFunc)(float,float,float,float);
	

	void CalcNNTour(float *&origVertX, float *&origVertY)
	{
		float aDist = 0.0f;
		int i, j;
		int *tour = (int*)malloc(numVerts * sizeof(int));
		float *newVertX = (float*)malloc( numVerts * sizeof(float) );
		float *newVertY = (float*)malloc( numVerts * sizeof(float) );


		for ( i = 0; i < numVerts; i++ )
		{
			tour[i] = i;
		}
		for ( i = 0; i < numVerts-1; i++ )
			{
				float nearD = 1e20f;
				int nearI;
				for ( j = i+1; j < numVerts; j++ )
				{
					if ( CalcEdgeDist(tour[i],tour[j]) < nearD )
					{
						nearD = CalcEdgeDist(tour[i],tour[j]);
						nearI = j;
					}
				}
				aDist += nearD;
				int swap = tour[i+1];
				tour[i+1] = tour[nearI];
				tour[nearI] = swap;
			}

			//the distance between the last two vertices is added to aDist, as well as the distance between the last and first vertex
			aDist += CalcEdgeDist(tour[numVerts-2],tour[numVerts-1]);
			aDist += CalcEdgeDist(tour[0],tour[numVerts-1]);

			nnDist = aDist;

			float sanityCheck = 0.0f;
			for ( i = 0; i < numVerts; i++ )
			{
				int i0 = tour[i];

				int i1 = tour[(i+1)%numVerts];
				//printf("\n SANITY CHECK i1: %d", i1);
				sanityCheck += CalcEdgeDist(i0,i1);
			}		
			for (i = 0; i < numVerts; i++)
			{
				newVertX[i] = origVertX[tour[i]];
				newVertY[i] = origVertY[tour[i]];
				/*
				for(j = 0; j < numVerts; j++)
				{ 
					if(i == j)
					{
						edgeDistNew[i][j] = 0;
					}
					else
					{
						//edgeDistNew[i][j] = CalcEdgeDist(tour[i],tour[j]);
					}	
				}*/
			}
			free(tour);
			//printf("X: %f, Y: %f\n", newVertX[10],newVertY[10]);
			//printf("X: %f, Y: %f\n", origVertX[10],origVertY[10]);
			//origVertX = newVertX;
			//origVertY = newVertY;			
			//printf("X: %f, Y: %f", origVertX[10],origVertY[10]);
	}
	
	void FillNNList(int iList)
	{	
		nearestNeighbour *newNN = (nearestNeighbour*)malloc(numNN * sizeof(nearestNeighbour));
		distSort *tempList = (distSort*)malloc((numVerts - 1) * sizeof(distSort));
		

			
		int count = 0;
		for (int i = 0; i < numVerts; i++)
		{
			if (i != iList)
			{				
				tempList[count].index = i;
				tempList[count].dist = CalcEdgeDist(i,iList);
				count++;					
			}
		}
		qsort(tempList, numVerts - 1, sizeof(distSort), nnComp);		
		for(int i = 0; i < numNN; i++)
		{
			nnIndexes[iList][i] = tempList[i].index;
		}
		int count2 = 0;
		int count3 = 0;
		for (int i = 0; i < numNN; i++)
		{		
			

			if (tempList[i].index != -1)
			{	
				nearestNeighbour _nn;
				_nn.nnMask = 0;
				_nn.vectIndex = tempList[i].index / _VECSIZE;;
				int remainder = tempList[i].index % _VECSIZE;
			
				int edgeDistCount = 0;
				//fills edgedistance - adds every value contained in a vectorIndex, not just nearest neighbours
				for(int j = _nn.vectIndex* _VECSIZE; j < (_nn.vectIndex * _VECSIZE) + _VECSIZE; j++ )
					{
					
						edgeDist[iList][(count3*_VECSIZE)+edgeDistCount] = CalcEdgeDist(iList,j);
						nnList[iList][(count3*_VECSIZE)+edgeDistCount] = j;

						//printf("%f \n", CalcEdgeDist(iList,j));
						//printf("%d, %d: %f \n",iList,(count3*_VECSIZE)+edgeDistCount),edgeDist[iList][(count3*_VECSIZE)+edgeDistCount];
						edgeDistCount++;
					}
				count3++;
					
				
				_nn.nnMask |=  1UL << remainder;
				tempList[i].index = -1;
				for (int j = 0; j < numNN; j++)
				{
					if (tempList[j].index / _VECSIZE == _nn.vectIndex && tempList[j].index != -1)
					{
						int remainder = tempList[j].index % _VECSIZE;
						_nn.nnMask |= 1UL << remainder;
						tempList[j].index = -1;
					}
				}		
				newNN[count2] = _nn;
				count2++;
			}
		}

		for (int i = count2; i < numNN; i++)
		{
			nearestNeighbour _nn;
			_nn.nnMask = -1;
			_nn.vectIndex = -1;
			newNN[count2] = _nn;
			count2++;
		}

		free(tempList);
		neighbourVectors[iList] = newNN;
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

	static float EucSqDistance( float x0, float y0, float x1, float y1 )
	{
		return (x0-x1)*(x0-x1) + (y0-y1)*(y0-y1);
	}

	static float CeilSqDistance(float x0, float y0, float x1, float y1)
	{
		return (x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1);
	}

	static int AttSqDistance(float x0, float y0, float x1, float y1)
	{
		float d = (x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1);
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

	float CalcEdgeDist(int pointA, int pointB)
	{
		return (float)distanceFunc(vertX[pointA], vertY[pointA], vertX[pointB], vertY[pointB]);
	}

	float CalcEdgeDistSquared(int x0, int y0, int x1, int y1)
	{
		return (float)distanceSquaredFunc(x0,y0,x1,y1);
	}

	int getXFromPoint(int point)
	{
		return vertX[point];
	}

	int getYFromPoint(int point)
	{
		return vertY[point];
	}


	void Init( const char *fileName, int nNearNeighbours )
	{
		int i, j, k;
		numVerts = 0;
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

		//this while loop reads through the TSP's header details to determine parameters
		while ( stillHeader )
		{
			
			//fscanf reads data from a stream (fp) and stores it in line. %s indicates that only strings are read from the file
			fscanf( fp, "%s", line );
			
			//ends the loop if NODE_COORD_SECTION is reached as this is the end of the header
			//strcmp compares the two strings given as parameters
			if ( strcmp(line,"NODE_COORD_SECTION") == 0 )
				stillHeader = 0;

			//determines the dimension of the TSP, dimension in this case being the number of nodes
			else if ( strcmp(line,"DIMENSION") == 0 || strcmp(line,"DIMENSION:") == 0 ) 
			{
				if (strcmp(line, "DIMENSION:") != 0) // may or may not be a space between word and colon, read the colon if so
					fscanf( fp, "%s", line );
				fscanf( fp, "%s", line );
				sscanf(line,"%d", &numVerts ); //saves the value of line, in this case the dimension, to &numverts
			}

			//determines the edge weight type, which specifies how the edge weights are determined
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

		free(line);

		//end of while loop - header has ended

		switch (edgeType)
		{
		case EDGE_WEIGHT_UNSUPPORTED:
			fprintf(stderr, "unsupported edge weight type in tsplib file\n");
			return;
		case EDGE_WEIGHT_EUC2D:
			distanceFunc = &EucDistance;
			distanceSquaredFunc = &EucSqDistance;
			break;
		case EDGE_WEIGHT_CEIL2D:
			distanceFunc = &CeilDistance;
			distanceSquaredFunc = &CeilSqDistance;
			break;
		case EDGE_WEIGHT_ATT:
			distanceFunc = &AttDistance;
			break;
		case EDGE_WEIGHT_GEO:
			distanceFunc = &GeoDistance;
			break;
		}

		// read in the vertices
		//allocates memory for the vertices - the space allocated is a matrix (number of vertices * sizeof(float), which is 4.
		vertX = (float*)malloc( numVerts * sizeof(float) );
		vertY = (float*)malloc( numVerts * sizeof(float) );
		int iVert = 0;
		//loop is performed for the previously specified number of vertices
		for ( i = 0; i < numVerts; i++ )
		{
			float x, y;
			//x and y are set as the x and y value of the next unread line in the TSP file
			fscanf( fp, "%d %e %e", &j, &x, &y );
			


			int found = 0;
			for ( int j = 0; j < iVert; j++ )
			{
				//printf("\n J:%d  X:%F  Y:%F", j, x, y);

				//performs distance function - this varies depending on edge_weight_type

				//printf("\n X:%f  Y:%f  vertX:%f  vertY:%f", x, y, vertX[i], vertY[j]);
				int dist = distanceFunc( x, y, vertX[i], vertY[j] );
				if ( dist == 0 )
				{
					//printf("\n FOUND X:%f  Y:%f  vertX:%f  vertY:%f", x, y, vertX[i], vertY[j]);
					found = 1;
					break;
				}
			}

			//adds vertices to vertX and vertY if they don't already exist in vertX and vertY
			if ( !found )
			{
				vertX[iVert] =  x;
				vertY[iVert] =  y;
				iVert++;
			}
		}
		fclose(fp);
		numVerts = iVert;
		

		
/*
		for ( i = 0; i < numVerts; i++ )
		{

			for ( j = 0; j < numVerts; j++ )
			{
				if (i == j)
					edgeDist[i][j] = 0;
				else
					edgeDist[i][j] = distanceFunc( vertX[i], vertY[i], vertX[j], vertY[j] );
			}
		}*/

		
		//printf("X: %f, Y: %f\n", vertX[0],vertY[0]);
		CalcNNTour(vertX,vertY);
		//printf("X: %f, Y: %f", vertX[0],vertY[0]);

		numNN = nNearNeighbours;
		if(numNN * _VECSIZE < numVerts)
		{
			matrixSize = numNN*_VECSIZE;
		}
		else
		{
			matrixSize = numVerts;
		}
		neighbourVectors = (nearestNeighbour**)malloc(numVerts * sizeof(nearestNeighbour**));
		nnList = (int**)malloc(numVerts * sizeof(int*));
		edgeDist = (float**)malloc( numVerts * sizeof( float* ) );
		nnIndexes = (int**)malloc(numVerts * sizeof(int*));

		


		#ifdef USE_OMP
		#pragma omp parallel for
		#endif
		for ( i = 0; i < numVerts; i++ )
		{
			edgeDist[i] = (float*)malloc( (matrixSize) * sizeof( float ) );
			//each entry in nnList is a pointer-to-int array of numNN * 4
			//nnList is filled with 20 nearest neighbours of each city
			nnList[i] = (int*)malloc( (matrixSize) * sizeof( int ) );
			nnIndexes[i] = (int*)malloc( numNN * sizeof( int ) );

			
  			FillNNList( i );
			
			  /*for(int j = 0; j < numNN * _VECSIZE; j++)
			  {
				  printf("\n%d, %d: %f",i,j,edgeDist[i][j]);
			  }*/
  			
			
		}


		/*for(int j = 0; j < 32; j++)
		{
			printf("%d ",nnHist[j]);
		}
		printf("\n");*/


		
	}

};
		
#endif
