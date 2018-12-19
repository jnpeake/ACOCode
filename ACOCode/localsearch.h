#pragma once
#include "tsp.h"
#include "ranluxgen.h"
#include <random>
#include <climits>
//
// DMACS local search. This code is adapted from the ACOTSP v1.03 code
// by Thomas Stuetzle - copyright 1999 Thomas Stuetzle, issued under
// Version 2 of the GNU General Public Licence.
//
// Adapted from ls.c by Huw Lloyd, June 2018.
//

class LocalSearch
{
	TSP *tsp;
	int n; // number of verts
	RanluxGenerator &rng; 
	
	int nn_ls;            /* maximal depth of nearest neighbour lists used in the 
			      local search */ 
	int dlb_flag = true;  /* flag indicating whether don't look bits are used. I recommend 
			      to always use it if local search is applied */

	int *random_vector;
    int *pos;               /* positions of cities in tour */ 
    int *dlb;               /* vector containing don't look bits */ 

void generate_random_permutation( int *r)
{
/*    
      FUNCTION:       generate a random permutation of the integers 0 .. n-1
      INPUT:          length of the array
   
   OUTPUT:         pointer to the random permutation
      (SIDE)EFFECTS:  the array holding the random permutation is allocated in this 
                      function. Don't forget to free again the memory!
      COMMENTS:       only needed by the local search procedures
*/
	int  i, help, node, tot_assigned = 0;
	double    rnd;

	std::uniform_real_distribution<float> dist(0.0f,1.0f);

	for ( i = 0 ; i < n; i++) 
		r[i] = i;
	for ( i = 0 ; i < n ; i++ ) {
		/* find (randomly) an index for a free unit */ 
		rnd  = rng.frand();
		if ( rnd >= 1.0f )
			rnd = 0.5f;
		node = (long int) (rnd  * (n - tot_assigned)); 
		help = r[i];
		r[i] = r[i+node];
		r[i+node] = help;
		tot_assigned++;
	}

	return;
}


public:
LocalSearch( TSP *tsp, RanluxGenerator &rng, int nn_ls ) : tsp(tsp), n(tsp->numVerts), rng(rng), nn_ls(nn_ls) 
{
    pos = (int*)malloc(n * sizeof(int));
    dlb = (int*)malloc(n * sizeof(int));
	random_vector = (int*)malloc(n * sizeof(int));  
}
~LocalSearch()
{
	free(pos); free(dlb); free(random_vector);
}

	void TwoOpt( int *tour );
	void TwoHalfOpt( int *tour );
	void ThreeOpt( int *tour );
};