#ifdef __cplusplus
extern "C" {
#endif



#include <stdio.h>
#include <math.h>
#include <R.h>

void sjevor(Sfloat *xpts, Sfloat *ypts, Sfloat *dims, char **popts,
	    Sfloat *info, int *sneighs,
	    Sfloat *ias,
	    int *del_ids2, Sfloat *del_lens2, Sfloat *del_angs2,
	    Sfloat *poly_pts2,
	    Sfloat *vx1, Sfloat *vy1, int *vertices,
	    int *pnpts,
	    int *limits, int *debug);


void sjevoradd(Sfloat *xpts, Sfloat *ypts, Sfloat *temp, int npts);






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
void find_internal_angles(Sfloat *xpts, Sfloat *ypts, int npts,
			  Sfloat *ias, int ias_max);


/* To sort the neighbours by distance, we use qsort().  This uses
 * the following simple structure. */
typedef struct keydist {
  int   key;			/* number of neighbour cell  */
  Sfloat dist;			/* distance from reference cell. */
} Keydist;

int keydist_cmp (Keydist *c1, Keydist *c2);
Keydist *dists;			/* dynamically allocated memory */

int sje_debug;			/* non-zero if we want debug output. */
int max_del_tris;		/* max# of del triangles per point (5). */

#ifdef __cplusplus
}
#endif

