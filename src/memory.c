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
  if ((t=malloc(n)) == (char *) 0)
    {    fprintf(stderr,"Insufficient memory processing site %d (%d bytes in use)\n",
		 siteidx, total_alloc);
    exit(-1);			/* sje: include -1. */
    };
  total_alloc += n;


  /* Keep hold of memory allocated here.  Since the Fortune code
   * doesn't seem to free any memory itself, we keep a list of all the
   * pointers allocated and then we can free them after we have
   * finished with the Fortune code.  We store them in a list.
   * Luckily, this function seems to be the only place where memory is
   * allocated.
   */
  
  /*printf("sje: allocated %d to total %d\n", n, total_alloc); */
  
  if (num_fortune_pointers >= MAX_FORTUNE_POINTERS) {
    Rprintf("bye bye... reached limit of fortune pointers %d\n",
	   num_fortune_pointers);
    exit(-1);
  } else {
    fortune_pointers[num_fortune_pointers++] = t;
  }
  return(t);
}
