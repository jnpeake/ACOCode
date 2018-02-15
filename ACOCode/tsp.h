#ifndef _TSP_INC_
#define _TSP_INC_
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <vector>

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
	int numVerts;
	int **nnList;
	int numNN;
	nearestNeighbour **neighbourVectors;
	nearestNeighbour *newNN;
	
	// nn list
	/*

	;
	nearestNeighbour *newNN;

	std::list<nearestNeighbour> nList;
	nList.push_back()
		for (nearestNeighbour item : nList)
		{

		}
	auto iter = nList.begin();
	while (iter != nList.end())
	{
		// item is *iter
		iter++;
	}
*/
	int (*distanceFunc)(float, float, float, float);
	
	/*
	void FillNNList( int iList )
	{
		//*tempList is a distSort array of size numverts-1 * 8
		//distSort is a struct declared earlier, with two parameters:
		//float dist;
		//int index;
		printf("\n Size of distSort %d ", sizeof(distSort));
		distSort *tempList = (distSort*)malloc( (numVerts-1) * sizeof( distSort ) );
		int count = 0;
		for ( int i = 0; i < numVerts; i++ )
		{
			//checks if i matches the iList parameter, which itself is the iterating value of another for loop
			//if false, adds a distSort struct to the array with the index i and a dist value of the value at i, iList of the edgeDist
			if ( i != iList )
			{
				tempList[count].index = i;
				tempList[count].dist = edgeDist[i][iList];
				//printf("\n Temp List %d data: Index: %d Dist: %f", count, i, edgeDist[i][iList]);
				count++;
			}
		}
		qsort( tempList, numVerts-1, sizeof(distSort), nnComp );
		for ( int i = 0; i < numNN; i++ )
		{
			//printf("\n Adding value %d to nnList at position %d,%d", tempList[i].index, iList, i);
			nnList[iList][i] = tempList[i].index;
		}
		free( tempList );
	}
	*/
	
	void FillNNList(int iList)
	{

			
			//emulated code

			//*tempList is a distSort array of size numverts-1 * 8
			//distSort is a struct declared earlier, with two parameters:
			//float dist;
			//int index;
			//printf("\n Size of distSort %d ", sizeof(distSort));
			nearestNeighbour *newNN = (nearestNeighbour*)malloc(numNN * sizeof(nearestNeighbour));
			distSort *tempList = (distSort*)malloc((numVerts - 1) * sizeof(distSort));
			int count = 0;
			for (int i = 0; i < numVerts; i++)
			{
				//checks if i matches the iList parameter, which itself is the iterating value of another for loop
				//if false, adds a distSort struct to the array with the index i and a dist value of the value at i, iList of the edgeDist

				//the purpose of this for loop is to get the list of weights for the vertex at iList - stay the same as default
				if (i != iList)
				{
						
					tempList[count].index = i;
					tempList[count].dist = edgeDist[i][iList];
					//printf("\n Temp List %d data: Index: %d Dist: %f", count, i, edgeDist[i][iList]);
					count++;
						
				}
			}
			qsort(tempList, numVerts - 1, sizeof(distSort), nnComp);

		
			
			int count2 = 0;
			/////now we have the sorted list - time to find the index of the matrices in the pheromone matrix
			for (int i = 0; i < numNN; i++)
			{

				if (tempList[i].index != -1)
				{
					nearestNeighbour _nn;
					_nn.nnMask = 0;
					_nn.vectIndex = tempList[i].index / 16;
					int remainder = tempList[i].index % 16;
					_nn.nnMask |=  1UL << remainder;
					tempList[i].index = -1;
					for (int j = 0; j < numNN; j++)
					{
						if (tempList[j].index / 16 == _nn.vectIndex && tempList[j].index != -1)
						{
							
							int remainder = tempList[j].index % 16;
							_nn.nnMask |= 1UL << remainder;
							tempList[j].index = -1;
						}
					}		
					//_nn.nnMask = (1 << _nn.nnMask) - 1;
					newNN[count2] = _nn;
					count2++;

					//printf("\n INDEX: %d \n", _nn.vectIndex);
					for (int j = 0; j < 16; j++)
					{
						//printf(" %d",mask[j]);	
					}
				}

				//printf("\n%d",i);
				
				//printf("\n Adding value %d to nnList at position %d,%d", tempList[i].index, iList, i);
				//nnList[iList][i] = tempList[i].index;
			}

			for (int i = count2; i < numNN; i++)
			{
				nearestNeighbour _nn;
				_nn.nnMask = NULL;
				_nn.vectIndex = -1;
				newNN[count2] = _nn;
				count2++;
			}
			free(tempList);
			neighbourVectors[iList] = newNN;
			//free(newNN);
		
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

		//fopen opens the specified file and associates it with a stream that can be identified by the specified pointer (*fp)
		FILE *fp = fopen( fileName, "r");

		if ( fp == NULL )
		{
			fprintf(stderr, "could not open %s for read\n", fileName );
			perror("Error");
			return;
		}

		//malloc = allocate memory block. In this case the block is 1KB (1024 Bytes). *line is a pointer to the beginning of this block.
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

		//end of while loop - header has ended

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
		//printf("\n Pre numVerts %d", numVerts);
		//number of vertices is updated in case duplicates were found in the previous step
		numVerts = iVert;
		//printf("\n Post numVerts %d", numVerts);
		

		// construct the edges
		// pointer-to-pointer-to-float matrix of numverts * 4
		edgeDist = (float**)malloc( numVerts * sizeof( float* ) );

		//each value in edgeDist array is a pointer to a float which is numverts * 4
		for ( i = 0; i < numVerts; i++ )
			edgeDist[i] = (float*)malloc( numVerts * sizeof( float ) );

		//this creates a matrix of numVerts*numVerts
		for ( i = 0; i < numVerts; i++ )
		{

			//distance between each city is calculated and added to edge matrix
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
		//numNN is set to the number of nearest neighbours passed into the program as a parameter
		numNN = nNearNeighbours;

		//nnList is a pointer-to-pointer-to-int array of numVerts * 4
		//nnList = (int **)malloc( numVerts * sizeof(int*) );
		neighbourVectors = (nearestNeighbour**)malloc(numVerts * sizeof(nearestNeighbour**));
		
		
		for ( i = 0; i < numVerts; i++ )
		{
			//each entry in nnList is a pointer-to-int array of numNN * 4
			//nnList is filled with 20 nearest neighbours of each city
			//nnList[i] = (int*)malloc( numNN * sizeof( int ) );
			FillNNList( i );
		}

		
	}

};
		
#endif
