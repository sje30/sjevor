#ifdef __cplusplus
extern "C" {
#endif



#include <stdio.h>
#include <math.h>



/* This is the expected maximum number of neighbours. */
/* 18 should be a conservative number, 13 normally works okay. */
#define MAX_NUM_NEIGHS 18

/* Definitions for the routines that will be called from Octave. */

void sjevor(Sfloat *xpts, Sfloat *ypts, Sfloat *dims, char **popts, 
	    Sfloat *info, int *sneighs, Sfloat *ias,
	    int *del_ids2, Sfloat *del_lens2, Sfloat *del_angs2,
	    int *npts);
//%input xpts(npts), ypts(npts), dims(4)
//%output info(npts,4), sneighs(npts, MAX_NUM_NEIGHS)

void sjevoradd(Sfloat *xpts, Sfloat *ypts, Sfloat *temp, int npts);
//%input xpts(npts), ypts(npts)
//%output temp(npts)  






/* General function definitions. */
void find_rejects(int npts);
int out_of_bounds(int v);

void init_neighs(int npts);
void find_neighs();
void add_neigh(int i, int j);
void write_neighs(int npts);
void find_nnd(Sfloat *xpts, Sfloat *ypts, int npts,
	      Sfloat *temp, int *sneighs);
void find_vertices(int npts);
void find_areas(int npts, Sfloat *temp);
void myfree(void *ptr);
void sje_readsites(Sfloat *xpts, Sfloat *ypts, int npts);

void find_internal_angles(Sfloat *xpts, Sfloat *ypts, int npts, Sfloat *ias);

/* To sort the neighbours by distance, we use qsort().  This uses
 * the following simple structure. */
typedef struct keydist {
  int   key;			/* number of neighbour cell  */
  Sfloat dist;			/* distance from reference cell. */
} Keydist;

int keydist_cmp (Keydist *c1, Keydist *c2);

/*need to specify this as a hard limit somehow MAX_NUM_NEIGHS];*/
#define MAX_DISTS MAX_NUM_NEIGHS
Keydist dists[MAX_DISTS]; 

int sje_debug;			/* non-zero if we want debug output. */


#ifdef __cplusplus
}
#endif

