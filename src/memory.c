#
#include "defs.h"
#include <stdio.h>
#include <stdlib.h>		/* stop malloc library compile error. */
void freeinit(fl, size)
     struct	Freelist *fl;
     int	size;
{
  fl -> head = (struct Freenode *) NULL;
  fl -> nodesize = size;
}

char *getfree(fl)
     struct	Freelist *fl;
{
  int i; struct Freenode *t;
  if(fl->head == (struct Freenode *) NULL)
    {	t =  (struct Freenode *) myalloc(sqrt_nsites * fl->nodesize);
    for(i=0; i<sqrt_nsites; i+=1) 	
      makefree((struct Freenode *)((char *)t+i*fl->nodesize), fl);
    };
  t = fl -> head;
  fl -> head = (fl -> head) -> nextfree;
  return((char *)t);
}



makefree(curr,fl)
     struct Freenode *curr;
     struct Freelist *fl;
{
  curr -> nextfree = fl -> head;
  fl -> head = curr;
}

int total_alloc;
char *myalloc(n)
     unsigned n;
{
  char *t;
  if ((t=R_alloc(n,1)) == (char *) 0)
    {    fprintf(stderr,"Insufficient memory processing site %d (%d bytes in use)\n",
		 siteidx, total_alloc);
    error("Error in myalloc from memory.c");
    };
  total_alloc += n;

  /* We used to keep a memory of the fortune pointers allocated, but
   * now instead we use R_alloc().  R takes care of freeing the memory
   * at the end of the .C call. This seems more sensible. */
  
  return(t);
}
