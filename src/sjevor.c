/****************************************************************************
***
***
*** $RCSfile$
*** 
*** Stephen Eglen
*** 
***
*** Created 17 Apr 2000
***
*** $Revision$
*** $Date$
****************************************************************************/



/*********************** + Purpose of File + ************************/
/* take a set of x,y coordinates and then call Voronoi code. */
/*********************** * Purpose of File * ************************/


/* -  Include Files - */


#include "defs.h"
#include "sjevor.h"
/* - Defines - */

/* - Function Declarations - */

void freeinit(struct Freelist *fl, int size);
void sje_readsites(float *xpts, float *ypts, int npts);
struct Site *nextone();
void geominit();
void voronoi(int triangulate, struct Site *(*nextsite)());

/* - Global Variables - */ 



/* - Start of Code  - */


void sjevor(float *xpts, float *ypts, float *temp, int npts)
{


  struct Site *(*next)();
  printf("%f %f\n", xpts[2], ypts[2]);

  /* Now that we have the data, we need to send it to the Voronoi code.
   */


  /* Initialise the data. */
  freeinit(&sfl, sizeof *sites);

  sje_readsites(xpts, ypts, npts);
  next = nextone;


  siteidx = 0;			/* defined globally in defs.h */

  geominit();
/*    if(plot) plotinit(); */

  voronoi(triangulate, next); 
}



void sjevoradd(float *xpts, float *ypts, float *temp, int npts)
{
  /* test function to see that matwrap is working okay. */
  int i;
  /* printf("%f %f\n", xpts[2], ypts[2]);*/
  for (i=0; i< npts; i++)
    temp[i] = xpts[i] + ypts[i];

}
