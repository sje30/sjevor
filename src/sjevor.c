/****************************************************************************
***
***
*** $RCSfile: sjevor.c,v $
*** 
*** Stephen Eglen
*** 
***
*** Created 17 Apr 2000
***
*** $Revision: 1.1 $
*** $Date: 2000/04/17 12:09:20 $
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


void sjevor(float *xpts, float *ypts,
	    float *temp, int *sneighs,
	    int npts)
{


  struct Site *(*next)();
  printf("%f %f\n", xpts[2], ypts[2]);

  /* Now that we have the data, we need to send it to the Voronoi code.
   */

  /* Set up my data structures for remembering things. */


  first_index = 1;		/* 0/1 offset problem. */
  ignore_rejects = 1;

  vnum_max = 4 * npts;		/* conservative guess. */
  
  vx = (float*)calloc(vnum_max, sizeof(float));
  if (! vx) { 
    printf("could not allocate space for vx\n");
    exit(-1);
  }

  vy = (float*)calloc(vnum_max, sizeof(float));
  if (! vy) { 
    printf("could not allocate space for vy\n");
    exit(-1);
  }
  vnum = 0;

  lnum_max = 5 * npts;		/* conservative guess. */

  la  = (float*)calloc(lnum_max, sizeof(float));
  lb  = (float*)calloc(lnum_max, sizeof(float));
  lc  = (float*)calloc(lnum_max, sizeof(float));
  lb1 = (int*)calloc(lnum_max, sizeof(int));
  lb2 = (int*)calloc(lnum_max, sizeof(int));
  lnum = 0;

  if ((!la) || (!lb) || (!lc) || (!lb1) || (!lb2)) {
    printf("could not allocate space for line data structures\n");
    exit(-1);
  }

  ednum_max = 5 * npts;
  el  = (int*)calloc(ednum_max, sizeof(int));
  ev1  = (int*)calloc(ednum_max, sizeof(int));
  ev2  = (int*)calloc(ednum_max, sizeof(int));

  if ((!el) || (!ev1) || (!ev2) ) {
    printf("could not allocate space for line data structures\n");
    exit(-1);
  }
  
  
  numpoints = (int*)calloc(npts, sizeof(int));
  if (! numpoints) { 
    printf("could not allocate space for numpoints\n");
    exit(-1);
  }

  reject = (int*)calloc(npts, sizeof(int));
  if (! reject) { 
    printf("could not allocate space for reject\n");
    exit(-1);
  }

  
  
  /*************************************************************/
  /* initialise the data. */
  freeinit(&sfl, sizeof *sites);

  sje_readsites(xpts, ypts, npts);
  next = nextone;


  siteidx = 0;			/* defined globally in defs.h */

  geominit();
/*    if(plot) plotinit(); */

  voronoi(triangulate, next);
  /*************************************************************/

  printf("read %d vertices\n", vnum);
  printf("read %d lines\n", lnum);
  printf("read %d edges\n", ednum);

  /* Post-process the memory structures. */



  /* find the rejects. */
  find_rejects(npts);


  init_neighs(npts); find_neighs(); write_neighs(npts);

  find_nnd(xpts, ypts, npts, temp, sneighs);


  /* Find the areas of the polygons. */
  find_vertices(npts);
  find_areas(npts, temp);
  
  /* Clear-up memory. */
  free(vx); free(vy);
  free(la); free(lb); free(lc); free(lb1); free(lb2);
  free(el); free(ev1); free(ev2);
  free(numpoints);
  free(reject);
}



void sjevoradd(float *xpts, float *ypts, float *temp, int npts)
{
  /* test function to see that matwrap is working okay. */
  int i;
  /* printf("%f %f\n", xpts[2], ypts[2]);*/
  for (i=0; i< npts; i++)
    temp[i] = xpts[i] + ypts[i];

}


void find_rejects(int npts)
{
  int i, edge;
  int line, v1, v2, p1, p2;
  int num_rejects = 0;
  
  for(i=0; i < npts; i++) {
    reject[i] = 0;
    numpoints[i] = 0;
  }

  /* better for me to set the min, max otherwise can get small problems in
     what is defined as a reject. */
  
  xmin = 0.0; xmax = 1.0; ymin = 0.0; ymax = 1.0;
  printf("find_rejects: boundary %f %f %f %f\n", xmin, xmax, ymin, ymax);


  for (edge=0; edge < ednum; edge++) {
    line = el[edge];
    v1 = ev1[edge];
    v2 = ev2[edge];

    p1 = lb1[line];     p2 = lb2[line];

    if ( (v1 == -1) || (v2 == -1) ||
	 out_of_bounds(v1) || out_of_bounds(v2)) {
      /* reject those cells. */
      reject[p1] = 1; reject[p2] = 1;
    } else {
      /* This edge is valid, so add vertices to both points. */
    }

      
  }

  /* print out how many rejects found. */
  for (i=0; i< npts; i++) {
    if (reject[i] ) {
      printf("%d ", (i + first_index));
      num_rejects++;
    }
  }
    printf("; total rejects: %d\n", num_rejects);
}



int out_of_bounds(int v)
{
  /* Return 1 if vertice v is out of bounds. */
  float x, y;
  x = vx[v]; y = vy[v];


  if ( ( x < xmin) || (x > xmax) || (y < ymin) || (y > ymax)) {
    return 1;
  } else {
    return 0;
  }
}

void init_neighs(int npts)
{
  /* initialise the neighbours data structures. */
  int i;

  numneighs = (int*)calloc(npts, sizeof(int));
  neighs = (int*)calloc(npts*MAX_NUM_NEIGHS, sizeof(int));
  for (i=0; i< npts; i++) {
    numneighs[i] = 0;
  }

}

void find_neighs()
{
  /* Find all of the neighbouring data points.
   * This is done by checking all the lines. */

  int l, pta, ptb;

    for (l=0; l < lnum; l++) {
	pta = lb1[l]; ptb = lb2[l];
	
	add_neigh(pta, ptb);
	add_neigh(ptb, pta);
    }
}


void add_neigh(int i, int j)
{
  /* Make points i and j neighbours of each other, if they are
   * not already neighbours.
   */
  
  int num;
  int looking, k;
  
  num = numneighs[i];
    
  looking = 1;
  k = 0;
  while(looking && k < num) {
    if (neighs[NIND(i,k)] == j) {
      looking = 0;
    } else {
      k++;
    }
  }
  /* We couldn't find j in the list, so we add it.*/
  if (looking == 1) {
    neighs[NIND(i,num)] = j;
    numneighs[i]++;
  }
}



void write_neighs(int npts)
{
  
  /* Write out the list of neighbours. */

  int i, n, v;
  FILE	*fp;
  char *file  = "neighs";
  fp = fopen( file, "w");
  if (! fp ) {
    printf("write_neighs: %s could not be opened for writing", file);
    exit(-1);
  }

  for(i=0; i<npts; i++) {
    if (ignore_rejects && reject[i]) {
      /* exclude this cell */
      fprintf(fp, "    \n");
    } else {
      for(n=0; n < numneighs[i]; n++) {
	v = neighs[NIND(i,n)] + first_index;
	fprintf(fp, "%d ", v);
      }
      fprintf(fp,"\n");
    }
  }
  fclose(fp);
}




void find_nnd(float *xpts, float *ypts, int npts,
	      float *temp, int *sneighs)
{
  /* Find the NND for each datapoint.
   * We can also optionally sort the neighbours according to distance.
   */

  FILE	*nndfp, *snfp;
  char  *file = "nnds";
  int   p;
  float mindist, dist;
  int   minidx, v;
  float px, py,  dx, dy, dist2,  x, y;
  int n, i;
  int   sort_neighs = 1;
  int   num_sneighs;
  nndfp = fopen( file, "w");
  if (!nndfp ) {
    printf("%s: %s could not be opened for writing",
	   "find_nnd", file);
    exit(-1);
  }
  /*
  if (sort_neighs) {
    dists = (float*)calloc(MAX_, sizeof(float));
    if (! dists) { 
      printf("%s: could not allocate space for dists\n", find_nnd);
      exit(-1);
    }
  }
  */

  if (sort_neighs) {
    snfp = fopen("sneighs", "w");
    if (!snfp) {
      printf("Could not open sneighs for writing\n");
      exit(-1);
    }

    /* initialise the sneighs array to store -1. */
    /* TODO: This loop can be optimised! */
    num_sneighs = npts*MAX_NUM_NEIGHS;
    for(i=0; i< num_sneighs; i++) {
      sneighs[i] = -1;
    }

  }
  
  for (p=0; p<npts; p++) {

    /* printf("finding nearest to point %d\n", p); */
    if (reject[p] && ignore_rejects) { 
      fputs("-1 -1\n", nndfp);
      	
      temp[RIND(p,0,npts)] =  p;
      temp[RIND(p,1,npts)] = -1;
      temp[RIND(p,2,npts)] = -1;
      temp[RIND(p,3,npts)] = -1;

      if (snfp) fputs("\n", snfp);
      
      continue;		/* jump to next p in the for loop. */
    }

    mindist = 99999999; minidx = -1;
    px = xpts[p];	py = ypts[p];


	/* For each cell, look amongst its neighbours
	 * to find its closest neighbour. */
	
	for(n=0; n < numneighs[p]; n++) {
	    i = neighs[NIND(p,n)];
	    x = xpts[i]; y = ypts[i];
	    dx = (x - px); dy = (y - py);
	    dist2 = (dx*dx) + (dy*dy);
	    
	    if (dist2 < mindist) {
		mindist = dist2;
		minidx = i;
	    }

  	    if (sort_neighs) {
	      dists[n].key = i + first_index;
	      dists[n].dist = dist2;
	    }
	}

	dist = (float)sqrt((double)mindist);
	v = minidx + first_index;
	fprintf(nndfp, "%.4f %d\n", dist, v);

	
	temp[RIND(p,0,npts)] = p;
	temp[RIND(p,1,npts)] = v;
	temp[RIND(p,2,npts)] = dist;
	temp[RIND(p,3,npts)] = -1;
	

	if (sort_neighs) {
	  /* Sort the neighbours according to distance from data point. */

	  qsort (dists, numneighs[p], sizeof (Keydist), keydist_cmp);

	  for(i=0; i<numneighs[p]; i++) {
	    if (snfp) fprintf(snfp,"%d ", dists[i].key);
	    sneighs[SNIND(p,i,npts)] = dists[i].key;
	  }
	  
	  fputs("\n", snfp);
	}

  }

  
/*      ##print "nearest to px py is pts[minidx][0] pts[minidx][1]\n"; */

/*      if (sort_neighs) { */
/*  	## output the sorted neighbour list too. */
/*  	&write_neighs(sorted_neighs_file); */
/*      } */

	fclose(nndfp);
	
	if (sort_neighs) {
	  free(dists);
	  if (snfp) fclose(snfp);
	}
}

/**********************************************************************/
/* sorting routines. */

int 
keydist_cmp (Keydist *c1, Keydist *c2)
{

  double temp =  (c1->dist - c2->dist);
  if (temp > 0)
    return 1;
  else if (temp < 0)
    return -1;
  else
    return 0;
  
}



/*********************************************************************/
/* Calculating area.
 *
 */

void find_vertices(int npts)
{

  int   e,p, l, v1, v2, p1, p2, numv;
  int   first, next, norder;
  int   done = -2;
  int looking, nverts, i;  
  int *vertices1, *vertices2;
  
  numvertices = (int*)calloc(npts, sizeof(int));
  vertices1 = (int*)calloc(max_numvertices*npts, sizeof(int));
  vertices2 = (int*)calloc(max_numvertices*npts, sizeof(int));
  verticeso = (int*)calloc(max_numvertices_o*npts, sizeof(int));

  if (!vertices1 || !vertices2) {
    printf("couldn't allocate vertices1/2\n");
    exit(-1);
  }

  
  if (! numvertices) { 
    printf("%s: could not allocate space for numvertices\n", __FUNCTION__);
    exit(-1);
  }
  /* todo: free vertices1/2. */
  /* free(numvertices); todo*/
  
  /* Initialise array to store the number of vertices for each data point. */
  for (p=0; p<npts; p++) {
    numvertices[p] = 0;
  }

    for (e=0; e < ednum; e++) {
      /* For each edge, there are two vertices, v1 and v2. */
      
      l = el[e];
      v1 = ev1[e]; v2 = ev2[e];
      
      /* These two vertices belong to two points, p1 and p2. */

      p1 = lb1[l]; p2 = lb2[l];

      p = p1;
      numv = numvertices[p];
      if ((numvertices[p] ++) > max_numvertices) {
	printf("exceeded max_numvertices: %d\n", max_numvertices);
	exit(-1);
      }

      vertices1[VIND(p,numv)] = v1;
      vertices2[VIND(p,numv)] = v2;

      p = p2;
      numv = numvertices[p];
      if ((numvertices[p] ++) > max_numvertices) {
	printf("exceeded max_numvertices: %d\n", max_numvertices);
	exit(-1);
      }

      vertices1[VIND(p,numv)] = v1;
      vertices2[VIND(p,numv)] = v2;
    }    

    for(i=0; i< numvertices[0]; i++) {
      printf("%d %d  ", vertices1[VIND(0,i)], vertices2[VIND(0,i)]);
    }
    printf("\n");
    
    /*
      vertices[p] stores the vertices for data point p.
      Each vertex will be included twice in the list, so we
      now need to order them.
      
      
      Polygon vertices must be ordered so that they are either clockwise
      or anticlockwise for the area calculation to work.  So for a polygon:
      A----B
      |    |
      C    D
      \E-/
      If the points are stored A D C B E, they must be reordered into:
      A B D E C A (making sure that the first point is also the last point.)
      
     The ordered vertices are stored in verticeso[p]
    */

    for (p=0; p< npts; p++) {
      printf("finding vertices for site %d\n", p);
      /* This won't work if the data point has any vertices 
       * out at infinity, so we should ignore them for now. */
	if (reject[p]) {
	  continue;		/* go on to next p. */
	} else {
	  ;			/* try point p */
	}

	/* Order the vertices for data point p.
	 * Go through the list of vertices, finding the next
	 * one until we return to the first point. */
	

	first = vertices1[VIND(p,0)]; vertices1[VIND(p,0)] = done;
  	 next = vertices2[VIND(p,0)]; vertices2[VIND(p,0)] = done;
	
	norder = 0;
	verticeso[VOIND(p,norder)] = first; norder++;
	looking = 1;
	nverts = numvertices[p];
	while (looking) {
	  printf("looking for %d\n", next);
	  verticeso[VOIND(p,norder)] = next; norder++;
	  for (i=0; i<nverts; i++) {
	    if (next == vertices1[VIND(p,i)]) {
		    next =  vertices2[VIND(p,i)];
		    vertices2[VIND(p,i)] = done;
		    i= nverts;
		} else if (next == vertices2[VIND(p,i)]) {
		    next =  vertices1[VIND(p,i)];
		    vertices1[VIND(p,i)] = done;
		    i= nverts;
		}
	    }
	    if (next == first) {
	      /* we've finished, since we're back to the first vertex. */
	      looking = 0;
	    }
	}
	/*  Now do the next data point. */
    }
    

    /*Print out the ordered vertices for each data point. */
#ifdef unused
    if (need_pvertices) {
	open(VERT, ">vertices_point_file");
	for (p=0; p<npts; p++) {
	    
	    if (reject{p}) { 
		print STDERR "vertice printing: ignoring point p\n";
		print VERT "-1\n";
		next;
	    }
	    n = numvertices[p];
	    for (i=0; i<n; i++) {
		v = verticeso[p][i] + first_index;
		print VERT "v ";
	    }
	    print VERT "\n";
	}
	close(VERT);
    }
#endif
}



void find_areas(int npts, float *temp)
{
  /* Find the area of each valid polygon.
   * Print out the ordered vertices for each data point.
   *
   * To get the areas, the vertices must first have been ordered,
   * and stored in verticeso.
   *
   * Algorithm taken from: 
   * http://www.mhri.edu.au/~pdb/geometry/polyarea/
   *
   * To check that the area calculation is correct, grab the
   * coordinates  in the following way and then check with
   * polyarea.m in Octave.
   * pfe 'd (`lines 150 150 vor_pvertices`); lines d d vertices'
   */

  int p, n, i, v, firstv, v1;
  float xi, yi, xi1, yi1;
  float sum, dsum;

  firstv = 0;			/* keep compiler quiet. */
  
  for (p=0; p < npts; p++) {
    if (reject[p]) {
      temp[RIND(p,3,npts)] = -1;
      continue;
    }
    sum = 0.0;
    n = numvertices[p];

    for (i=0; i < n; i++) {
      v = verticeso[VOIND(p,i)];
      
      if (i == 0) {
	firstv = v;
      }
	    
      /* When calculating area, the first and last point must be
       * the same.  Verticeso doesn't store the same point twice,
       * so we need to keep hold of the first vertice in firstv.
       */

      if (i == n-1 ) {
	v1 = firstv;
      } else {
	v1 = verticeso[VOIND(p,i+1)];
      }
      
      xi  = vx[v];  yi = vy[v];
      xi1 = vx[v1]; yi1 = vy[v1];
      
      dsum = (xi * yi1) - (xi1 * yi);
      sum += dsum;
      /*printf("%d %d %f\n", p, i, dsum);*/
    }
    sum /= 2.0;
    sum = fabs(sum);	 /* Ensure area is always positive.  Depends
			  * on ordering of coordinates. */

    /*v = p+first_index;*/
    temp[RIND(p,3,npts)] = sum;
  }
  
}
