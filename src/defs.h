#ifndef NULL
#define NULL 0
#endif
#define DELETED -2

int triangulate, sorted, plot, debug;

struct	Freenode	{
struct	Freenode	*nextfree;
};
struct	Freelist	{
struct	Freenode	*head;
int			nodesize;
};
char *getfree();
/* char *malloc(); */ /* sje: comment out */
char *myalloc();

float xmin, xmax, ymin, ymax, deltax, deltay;


struct Point	{
float x,y;
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
float		a,b,c;
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
float dist();
struct Point PQ_min();
struct Halfedge *PQextractmin();
struct Edge *bisect();

struct Halfedge {
struct Halfedge	*ELleft, *ELright;
struct Edge	*ELedge;
int		ELrefcnt;
char		ELpm;
struct	Site	*vertex;
float		ystar;
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




/* sje defines. */
float	*vx, *vy; int     vnum, vnum_max;

int lnum, lnum_max;
float *la, *lb, *lc;
int   *lb1, *lb2;

int ednum, ednum_max;
int *el, *ev1, *ev2;

int	*numpoints;
int     *reject;		/* reject[i] is 1 iff site i is a reject. */


int first_index;		/* for 0/1 offset problem when printing out. */


int *numneighs;
#define MAX_NUM_NEIGHS 13
int *neighs;

/* S is the site number and N is the nth neighbour so far of that site. */
#define NIND(S,N) ( (S*MAX_NUM_NEIGHS) + N)


/* S is the site number and N is the nth output value for that site. */
#define RIND(S,N,NPTS) ( (N*NPTS) + S)


int ignore_rejects;

/* Indexing into the sorted neighs array.
 * S is the site number and N is the nth nearest neighbour of that site. */
#define SNIND(S,N,NPTS) ( (N*NPTS) + S)

/* S is the site number and V is the vth vertice for that site. */
int max_numvertices = 20;	/*  TODO just a guess! */
#define VIND(S,V) ( (S*max_numvertices) + V)




int *verticeso;
int max_numvertices_o = 20;
/* S is the site number and V is the vth vertice for that site. */
#define VOIND(S,V) ( (S*max_numvertices) + V)

int	*numvertices;
