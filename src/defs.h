#ifndef NULL
#define NULL 0
#endif
#define DELETED -2

typedef double Sfloat;		/* from S.h? */

int triangulate, sorted, plot, debug;

struct	Freenode	{
  struct Freenode *nextfree;
};
struct	Freelist {
  struct Freenode	*head;
  int	        	nodesize;
};

char *getfree();
/* char *malloc(); */ /* sje: comment out */
char *myalloc();

Sfloat xmin, xmax, ymin, ymax, deltax, deltay;


struct Point	{
Sfloat x,y;
};

/* structure used both for sites and for vertices */
struct Site	{
struct	Point	coord;
int		sitenbr;
int		refcnt;
};


struct	Site	*sites;
int		nsites;
int		siteidx;
int		sqrt_nsites;
int		nvertices;
struct 	Freelist sfl;
struct	Site	*bottomsite;


struct Edge	{
Sfloat		a,b,c;
struct	Site 	*ep[2];
struct	Site	*reg[2];
int		edgenbr;
};
#define le 0
#define re 1
int nedges;
struct	Freelist efl;

int has_endpoint(),right_of();
struct Site *intersect();
Sfloat dist();
struct Point PQ_min();
struct Halfedge *PQextractmin();
struct Edge *bisect();

struct Halfedge {
struct Halfedge	*ELleft, *ELright;
struct Edge	*ELedge;
int		ELrefcnt;
char		ELpm;
struct	Site	*vertex;
Sfloat		ystar;
struct	Halfedge *PQnext;
};

struct   Freelist	hfl;
struct	Halfedge *ELleftend, *ELrightend;
int 	ELhashsize;
struct	Halfedge **ELhash;
struct	Halfedge *HEcreate(), *ELleft(), *ELright(), *ELleftbnd();
struct	Site *leftreg(), *rightreg();


int PQhashsize;
struct	Halfedge *PQhash;
struct	Halfedge *PQfind();
int PQcount;
int PQmin;
int PQempty();



/******************************************************************/
/* sje defines. */

#define RAD_TO_DEG  57.29577951308232

Sfloat	*vx, *vy; int     vnum, vnum_max;

int lnum, lnum_max;
Sfloat *la, *lb, *lc;
int   *lb1, *lb2;

int ednum, ednum_max;
int *el, *ev1, *ev2;

int del_idn, del_idmax;
int *del_ids;			/* pointer where Delaunay info can be stored.*/
Sfloat *del_lens, *del_angs;

int     *reject;		/* reject[s] is 1 iff site S is a reject. */


int first_index;		/* for 0/1 offset problem when printing out. */


int *numneighs;			/* numneighs[s] = number of neighbours of S. */
int *neighs;			/* 2.d row-major array (normal C).
				 *  neighs(NIND(S,N)) stores the Nth
				 *  neighbour of site S.*/

/* S is the site number and N is the nth neighbour so far of that site. */
#define NIND(S,N) ( (S*MAX_NUM_NEIGHS) + N)


/* This is the index into the output `info'.  info is a 2-d array
 * such that info[RIND(S,N,NPTS)] stores the Nth piece of info for site S.
 * This is ordered column-wise, since we return this to Octave.
 */
/* S is the site number and N is the nth output value for that site. */
#define RIND(S,N,NPTS) ( (N*NPTS) + S)


int ignore_rejects;		/* non-zero if we want to calculate nnd
				 * and area for cells at the border. */

/* Indexing into the sorted neighs array. sneighs[SNIND(S,N,NPTS)]
 * S is the site number and N is the nth nearest neighbour of that site.
 * The number of sites, NPTS, is needed since this is a column-major array that * gets returned to Octave.
 */
#define SNIND(S,N,NPTS) ( (N*NPTS) + S)

/* vertices1[VIND(S,N)] stores the index number of the Nth vertice found
 * for site S.  Both vertices1, vertices2 are 2-d row-major arrays.
 */
int max_numvertices;
#define VIND(S,V) ( (S*max_numvertices) + V)



int *verticeso;
int max_numvertices_o;
#define VOIND(S,V) ( (S*max_numvertices) + V)
/* verticeso[VOIND(S,N)] stores the index number of the Nth (ordered)
 * vertice found for site S.
 */


int	*numvertices;
/* 1-d array: numvertices[s] stores the number of vertices found so far
 * for site S. */


Sfloat    sje_minx, sje_maxx, sje_miny, sje_maxy;
/* min and max values of the field being processed. */


/* Keep a record of pointers allocated by Fortune. */
int num_fortune_pointers;

/* For a dmin mosaic with 5000 pts, it used about 300 pointers, so this
 * should be more than enough.
 */
#define MAX_FORTUNE_POINTERS 1000

void *fortune_pointers[MAX_FORTUNE_POINTERS];


/* Switches for controlling what we calculate in voronoi code. */
int need_areas;
int sort_neighs;
