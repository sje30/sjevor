#ifdef __cplusplus
extern "C" {
#endif

#ifdef silence_unused_code /* silence emacs for indenting purposes. */
}
#endif

#include <stdio.h>
#include <math.h>


void sjevor(float *xpts, float *ypts, 
	    float *temp, int *sneighs,
	    int npts);
//%input xpts(npts), ypts(npts)
//%output temp(npts,4), sneighs(npts, 13)

void sjevoradd(float *xpts, float *ypts, float *temp, int npts);
//%input xpts(npts), ypts(npts)
//%output temp(npts)  




void find_rejects(int npts);
int out_of_bounds(int v);

void init_neighs(int npts);
void find_neighs();
void add_neigh(int i, int j);
void write_neighs(int npts);
void find_nnd(float *xpts, float *ypts, int npts,
	      float *temp, int *sneighs);

typedef struct keydist {
  int   key;
  float dist;
} Keydist;


/*need to specify this as a hard limit somehow MAX_NUM_NEIGHS];*/
Keydist dists[20]; 


int keydist_cmp (Keydist *c1, Keydist *c2);



void find_vertices(int npts);
void find_areas(int npts, float *temp);
#ifdef __cplusplus
}
#endif

