#
#include "defs.h"
#include <stdio.h>
#include <math.h>
/* for those who don't have Cherry's plot */
/* #include <plot.h> */
openpl(){}
line(){}
circle(){}
range(){}

Sfloat pxmin, pxmax, pymin, pymax, cradius;


#undef sje_print_vorinfo	/* Define this if you want the output
				 * to stdout that is normally shown. */
out_bisector(e)
     struct Edge *e;
{
  if(triangulate & plot &!debug)
    line(e->reg[0]->coord.x, e->reg[0]->coord.y, 
	 e->reg[1]->coord.x, e->reg[1]->coord.y);
  if(!triangulate & !plot &!debug)
#ifdef sjetemp
    printf("l %f %f %f\n", e->a, e->b, e->c); /* sje: add \n */
#endif

  /* sje -- printout bisecting info as well. */
#ifdef sje_print_vorinfo  
  printf("l %f %f %f %d %d\n", 
	 e->a, e->b, e->c, e->reg[le]->sitenbr, e->reg[re]->sitenbr);
#endif
  
  if(debug)
    printf("line(%d) %gx+%gy=%g, bisecting %d %d\n", e->edgenbr,
	   e->a, e->b, e->c, e->reg[le]->sitenbr, e->reg[re]->sitenbr);


  la[lnum] = e->a;
  lb[lnum] = e->b;
  lc[lnum] = e->c;
  lb1[lnum] = e->reg[le]->sitenbr;
  lb2[lnum] = e->reg[re]->sitenbr;

  lnum++;
  if (lnum >= lnum_max) {
    printf("%s:%d lnum_max (%d) reached\n",
	   __FILE__, __LINE__, lnum_max);
    exit(-1);
  }
  
}


out_ep(e)
     struct Edge *e;
{

  /* sje: always clip_line to get the polygon indexes computed.
  if(!triangulate & plot)
  */
  clip_line(e);

#ifdef sje_print_vorinfo
  if(!triangulate & !plot)
    {	printf("e %d", e->edgenbr);
    printf(" %d ", e->ep[le] != (struct Site *)NULL ? e->ep[le]->sitenbr : -1);
    printf("%d\n", e->ep[re] != (struct Site *)NULL ? e->ep[re]->sitenbr : -1);
    };
#endif
  
  el[ednum] = e ->edgenbr;
  ev1[ednum] = e->ep[le] != (struct Site *)NULL ? e->ep[le]->sitenbr : -1;
  ev2[ednum] = e->ep[re] != (struct Site *)NULL ? e->ep[re]->sitenbr : -1;

  ednum++;

  if (ednum >= ednum_max) {
    printf("%s:%d ednum_max (%d) reached\n",
	   __FILE__, __LINE__, ednum_max);
    exit(-1);
  }

}

out_vertex(v)
     struct Site *v;
{
  vx[vnum] = v->coord.x;   vy[vnum] = v->coord.y;
  vnum++;

  if (vnum >= vnum_max) {
    printf("%s:%d vnum_max (%d) reached\n",
	   __FILE__, __LINE__, vnum_max);
    exit(-1);
  }

#ifdef sje_print_vorinfo
  if(!triangulate & !plot &!debug)
    printf ("v %f %f\n", v->coord.x, v->coord.y);
  if(debug)
    printf("vertex(%d) at %f %f\n", v->sitenbr, v->coord.x, v->coord.y);
#endif
}


out_site(s)
     struct Site *s;
{
  ;

#ifdef sje_print_vorinfo  
  if(!triangulate & plot & !debug)
    circle (s->coord.x, s->coord.y, cradius);
  if(!triangulate & !plot & !debug)
    printf("s %f %f\n", s->coord.x, s->coord.y);
  if(debug)
    printf("site (%d) at %f %f\n", s->sitenbr, s->coord.x, s->coord.y);
#endif
}


#define MYSQR(A,B) ((A*A) + (B*B))
out_triple(s1, s2, s3)
     struct Site *s1, *s2, *s3;
{
  Sfloat a, b, c, a2, b2, c2;	/* side lengths and their squares */
  Sfloat theta_a, theta_b, theta_c;

  int offset = 1;		/* for converting from 0-based to 1-based. */
  /* Save the triangulation information */
  if (del_idn+3 > del_idmax ) {
    printf("%s:%d del_idmax (%d) reached\n",
	   __FILE__, __LINE__, del_idmax);
    exit(-1);
  } else {
    /* Store the ids of the triangle. */
    del_ids[del_idn] = s1->sitenbr + offset;
    del_ids[del_idn+1] = s2->sitenbr + offset;
    del_ids[del_idn+2] = s3->sitenbr + offset;

    /* Calculate the square of the lengths, and the lengths. */
    a2 = MYSQR( (s1->coord.x - s2->coord.x), (s1->coord.y - s2->coord.y));
    b2 = MYSQR( (s2->coord.x - s3->coord.x), (s2->coord.y - s3->coord.y));
    c2 = MYSQR( (s3->coord.x - s1->coord.x), (s3->coord.y - s1->coord.y));
    a = sqrt(a2); b = sqrt(b2); c = sqrt(c2);

    /* Calculate the angles of Delaunay triangle using cosine rule. */
    theta_b = acos( (a2 + c2 - b2) / (2*a*c)) * RAD_TO_DEG;
    theta_a = acos( (b2 + c2 - a2) / (2*b*c)) * RAD_TO_DEG;
    theta_c = acos( (a2 + b2 - c2) / (2*a*b)) * RAD_TO_DEG;

    del_lens[del_idn  ] = a; del_angs[del_idn  ] = theta_a;
    del_lens[del_idn+1] = b; del_angs[del_idn+1] = theta_b;
    del_lens[del_idn+2] = c; del_angs[del_idn+2] = theta_c;

    
    del_idn +=3;
    
  }
    
  if(triangulate & !plot &!debug)
    printf("%d %d %d\n", s1->sitenbr, s2->sitenbr, s3->sitenbr);
  if(debug)
    printf("circle through left=%d right=%d bottom=%d\n", 
	   s1->sitenbr, s2->sitenbr, s3->sitenbr);
}



plotinit()
{
  Sfloat dx,dy,d;

  dy = ymax - ymin;
  dx = xmax - xmin;
  d = ( dx > dy ? dx : dy) * 1.1;
  pxmin = xmin - (d-dx)/2.0;
  pxmax = xmax + (d-dx)/2.0;
  pymin = ymin - (d-dy)/2.0;
  pymax = ymax + (d-dy)/2.0;
  cradius = (pxmax - pxmin)/350.0;
  openpl();
  range(pxmin, pymin, pxmax, pymax);
}


int clip_line(e)
     struct Edge *e;
{
  struct Site *s1, *s2;
  Sfloat x1,x2,y1,y2;

  if(e -> a == 1.0 && e ->b >= 0.0)
    {	s1 = e -> ep[1];
    s2 = e -> ep[0];
    }
  else 
    {	s1 = e -> ep[0];
    s2 = e -> ep[1];
    };

  if(e -> a == 1.0)
    {
      y1 = pymin;
      if (s1!=(struct Site *)NULL && s1->coord.y > pymin)
	y1 = s1->coord.y;
      if(y1>pymax) return;
      x1 = e -> c - e -> b * y1;
      y2 = pymax;
      if (s2!=(struct Site *)NULL && s2->coord.y < pymax) 
	y2 = s2->coord.y;
      if(y2<pymin) return(0);
      x2 = e -> c - e -> b * y2;
      if ((x1> pxmax & x2>pxmax) | (x1<pxmin&x2<pxmin)) return;
      if(x1> pxmax)
	{	x1 = pxmax; y1 = (e -> c - x1)/e -> b;};
      if(x1<pxmin)
	{	x1 = pxmin; y1 = (e -> c - x1)/e -> b;};
      if(x2>pxmax)
	{	x2 = pxmax; y2 = (e -> c - x2)/e -> b;};
      if(x2<pxmin)
	{	x2 = pxmin; y2 = (e -> c - x2)/e -> b;};
    }
  else
    {
      x1 = pxmin;
      if (s1!=(struct Site *)NULL && s1->coord.x > pxmin) 
	x1 = s1->coord.x;
      if(x1>pxmax) return(0);
      y1 = e -> c - e -> a * x1;
      x2 = pxmax;
      if (s2!=(struct Site *)NULL && s2->coord.x < pxmax) 
	x2 = s2->coord.x;
      if(x2<pxmin) return(0);
      y2 = e -> c - e -> a * x2;
      if ((y1> pymax & y2>pymax) | (y1<pymin&y2<pymin)) return(0);
      if(y1> pymax)
	{	y1 = pymax; x1 = (e -> c - y1)/e -> a;};
      if(y1<pymin)
	{	y1 = pymin; x1 = (e -> c - y1)/e -> a;};
      if(y2>pymax)
	{	y2 = pymax; x2 = (e -> c - y2)/e -> a;};
      if(y2<pymin)
	{	y2 = pymin; x2 = (e -> c - y2)/e -> a;};
    };
  /* Rather than output line, just store the points.
  line(x1,y1,x2,y2);
  */


  if (poly_idn+4 > poly_idmax ) {
    Rprintf("%s:%d poly_idmax (%d) reached\n",
	   __FILE__, __LINE__, poly_idmax);
    Rprintf("oh oh -- can we get back to R?  \n");
    /* SJE: have a problem -- error() causes R to crash, not good. */

  } else {
    /* Store the ends of the line bisector. */
    poly_pts[poly_idn++] = x1;
    poly_pts[poly_idn++] = y1;
    poly_pts[poly_idn++] = x2;
    poly_pts[poly_idn++] = y2;
  }
}
