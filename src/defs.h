#ifndef NULL
#define NULL 0
#endif
#define DELETED -2

#include <R.h>



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




struct Point	{
Sfloat x,y;
};

/* structure used both for sites and for vertices */
struct Site	{
struct	Point	coord;
int		sitenbr;
int		refcnt;
};

struct Edge	{
Sfloat		a,b,c;
struct	Site 	*ep[2];
struct	Site	*reg[2];
int		edgenbr;
};
#define le 0
#define re 1

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


struct	Halfedge *HEcreate(), *ELleft(), *ELright(), *ELleftbnd();
struct	Site *leftreg(), *rightreg();





/******************************************************************/
/* sje defines. */

#define RAD_TO_DEG  57.29577951308232
