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
*** $Revision: 1.8 $
*** $Date: 2002/01/04 01:57:24 $
****************************************************************************/



/*********************** + Purpose of File + ************************/
/* Take a set of x,y coordinates and then call Voronoi code. */
/*********************** * Purpose of File * ************************/


/* -  Include Files - */

#include <strings.h>
#include "defs.h"
#include "sjevor.h"
/* - Defines - */

/* - Function Declarations - */

/* These are the Fortune functions that I use. */
void freeinit(struct Freelist *fl, int size);
struct Site *nextone();
void geominit();
void voronoi(int triangulate, struct Site *(*nextsite)());



/* - Global Variables - */ 




/* - Start of Code  - */


void sjevor(Sfloat *xpts, Sfloat *ypts, Sfloat *dims, char **popts,
	    Sfloat *info, int *sneighs, Sfloat *ias,
	    int *del_ids2, Sfloat *del_lens2, Sfloat *del_angs2,
	    Sfloat *poly_pts2, int *poly_npts,
	    Sfloat *vx1, Sfloat *vy1, int *vertices,
	    int *pnpts)
{

  int i, npts, j;
  struct Site *(*next)();

  char *opts;

  opts = *popts;

  /* Set all pointers initially to be NULL; seems like sometimes they
   * could be initialised to rogue values.  Other pointers are okay,
   * since they will always get allocated.
   */

  npts = *pnpts;
  numvertices = NULL; verticeso = NULL;
  
  /* Parse options so that we set appropriate flags. */
  need_areas     = (strchr(opts, 'a')) ? 1:0;
  sort_neighs    = (strchr(opts, 's')) ? 1:0;
  ignore_rejects = (strchr(opts, 'i')) ? 0:1;
  plot           = (strchr(opts, 'p')) ? 0:1;

  /* printf("need areas: %d sort_neighs %d ignore_rejects %d\n",
     need_areas, sort_neighs, ignore_rejects); */


  
  /* Set up my data structures for remembering things. */


  first_index = 1;		/* 0/1 offset problem in C vs R/Octave arrays*/

  num_fortune_pointers = 0;	/* Set this to zero before we
				 * call and Fortune routines.
				 */

  /* Conservative guesses for the amount of space to allocate for
   * each structure.
   */

  /* Number of vertices should equal the number of neighbours. */
  max_numvertices = MAX_NUM_NEIGHS; 	/*  TODO just a guess! */
  max_numvertices_o = max_numvertices;
  vnum_max = MAX_DEL_TRIS * npts;
  lnum_max = MAX_DEL_TRIS * npts;		/* conservative guess. */
  ednum_max = MAX_DEL_TRIS * npts;

  /* e.g. for `t' data set of 100 points, we had 187 vertices, 286
   * lines and 286 edges. */

  /* Define "alloc_vx" if we wish the C program to allocate (& free) the
   * memory for the vertice positions (vx, vy) and the vertice numbers for
   * each site (verticeso)
   */
#ifdef alloc_vx
  vx = (Sfloat*)calloc(vnum_max, sizeof(Sfloat));
  if (! vx) { 
    printf("could not allocate space for vx %d\n", vnum_max);
    exit(-1);
  }

  vy = (Sfloat*)calloc(vnum_max, sizeof(Sfloat));
  if (! vy) { 
    printf("could not allocate space for vy\n");
    exit(-1);
  }
#else
  /* R already allocates this memory to store the vertices. */
  vx = vx1; vy = vy1; verticeso = vertices;
#endif
  
  vnum = 0;

  la  = (Sfloat*)calloc(lnum_max, sizeof(Sfloat));
  lb  = (Sfloat*)calloc(lnum_max, sizeof(Sfloat));
  lc  = (Sfloat*)calloc(lnum_max, sizeof(Sfloat));
  lb1 = (int*)calloc(lnum_max, sizeof(int));
  lb2 = (int*)calloc(lnum_max, sizeof(int));
  lnum = 0;

  if ((!la) || (!lb) || (!lc) || (!lb1) || (!lb2)) {
    printf("could not allocate space for line data structures\n");
    exit(-1);
  }

  el  = (int*)calloc(ednum_max, sizeof(int));
  ev1  = (int*)calloc(ednum_max, sizeof(int));
  ev2  = (int*)calloc(ednum_max, sizeof(int));
  ednum = 0;
  
  if ((!el) || (!ev1) || (!ev2) ) {
    printf("could not allocate space for line data structures\n");
    exit(-1);
  }
  
  

  reject = (int*)calloc(npts, sizeof(int));
  if (! reject) { 
    printf("could not allocate space for reject\n");
    exit(-1);
  }

  sje_minx = dims[0]; sje_maxx = dims[1];
  sje_miny = dims[2]; sje_maxy = dims[3];

  del_ids = del_ids2;
  del_lens = del_lens2; del_angs = del_angs2;

  del_idn = 0;			/* first free space */
  del_idmax = MAX_DEL_TRIS * 3 * npts;

  poly_idn = 0; poly_idmax = *poly_npts;
  poly_pts = poly_pts2;
  
  /*************************************************************/
  /* initialise the data. */
  freeinit(&sfl, sizeof *sites);

  sje_readsites(xpts, ypts, npts);
  next = nextone;


  siteidx = 0;			/* defined globally in defs.h */

  geominit();
  if(plot) plotinit();		/* may need to initialise pxmin etc */

  voronoi(triangulate, next);

  /*************************************************************/
  if (sje_debug) {
    printf("read %d vertices %d lines %d edges\n", vnum, lnum, ednum);
  }


  /* Post-process the memory structures. */

  /* Now that we have finished with Fortune code, remove the pointers.
   * See memory.c for details. */

  /*printf("%d fortune pointers allocated\n", num_fortune_pointers);*/
  for(i=0; i< num_fortune_pointers; i++) {
    free(fortune_pointers[i]);
  }


  /* Initialise the info array. */
  for (i=0; i<npts;i++) {
    info[RIND(i,0,npts)] =  (Sfloat)(i + first_index); /* index num */
    info[RIND(i,1,npts)] = -1.00;	/* id of nearest neigh */
    info[RIND(i,2,npts)] = -1.00;	/* distance to nearest neigh */
    info[RIND(i,3,npts)] = -1.00;	/* area of polygon of this unit */
  }


    
  /* find the rejects. */
  find_rejects(npts);


  init_neighs(npts); find_neighs();
  /*write_neighs(npts);*/

  find_nnd(xpts, ypts, npts, info, sneighs);


  /* Find the areas of the polygons. */
  if (need_areas) {
    find_vertices(npts);
    find_areas(npts, info);
    find_internal_angles(xpts, ypts, npts, ias);

    /* When returning vertice numbers, increase numbers by 1 to account
     * for 0/1 array problem.  This needed only if we want to return
     * the vertices.
     */
    for(i=0; i< npts; i++) {
      for(j=0; j< numvertices[i]; j++)
	verticeso[VOIND(i,j)] += first_index;
    }
  }

  del_ids2[del_idn] = -1;	/* Add -1 terminator so we know
				 * how many Delaunay triangles were found.
				 */

  *poly_npts = poly_idn;

  
  /* Clear-up memory. */
#ifdef alloc_vx
  myfree(vx); myfree(vy);
  myfree(verticeso);
#endif
  
  myfree(la); myfree(lb); myfree(lc); myfree(lb1); myfree(lb2);
  myfree(el); myfree(ev1); myfree(ev2);
  myfree(reject);

  myfree(numneighs); myfree(neighs);
  myfree(numvertices); 
}


void myfree(void *ptr)
{
  /* free the memory, but check first that the ptr is not NULL. */
  if (ptr) free(ptr);
}



void find_rejects(int npts)
{
  int i, edge;
  int line, v1, v2, p1, p2;
  for(i=0; i < npts; i++) {
    reject[i] = 0;
  }

  if (sje_debug) {
    printf("find_rejects: boundary %f %f %f %f\n", sje_minx, sje_maxx,
	   sje_miny, sje_maxy);
  }

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

#ifdef print_rejects
  /* print out how many rejects found. */
  /*int num_rejects = 0;*/
    
  for (i=0; i< npts; i++) {
    if (reject[i] ) {
      printf("%d ", (i + first_index));
      num_rejects++;
    }
  }
  printf("; total rejects: %d\n", num_rejects);
#endif
  
}



int out_of_bounds(int v)
{
  /* Return 1 if vertice v is out of bounds. */
  Sfloat x, y;
  x = vx[v]; y = vy[v];


  if ( ( x < sje_minx) || (x > sje_maxx) ||
       (y <  sje_miny) || (y > sje_maxy)) {
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
    if (numneighs[i] >= MAX_NUM_NEIGHS) {
      printf("%s:%d maximum number of neighbours (%d) exceeded\n",
	     __FILE__, __LINE__, numneighs[i]);
      exit(-1);
    }
      
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




void find_nnd(Sfloat *xpts, Sfloat *ypts, int npts,
	      Sfloat *info, int *sneighs)
{
  /* Find the NND for each datapoint.
   * We can also optionally sort the neighbours according to distance.
   */

  FILE	*nndfp = NULL, *snfp = NULL;
  char  *file = "nnds";
  int   p;
  Sfloat mindist, dist;
  int   minidx, v;
  Sfloat px, py,  dx, dy, dist2,  x, y;
  int n, i;
  int   num_sneighs;
  int	first_check;


  if (0 && !(nndfp = fopen( file, "w"))) {
    printf("%s: %s could not be opened for writing",
	   "find_nnd", file);
    exit(-1);
  }

  if (sort_neighs) {
    if (0 && !(snfp = fopen("sneighs", "w"))) {
      printf("Could not open sneighs for writing\n");
      exit(-1);
    }
  }

  /* initialise the sneighs array to store -1. */
  /* TODO: This loop can be optimised! */
  /* If we are not sorting the neighbours, need to return -1 anyway. */
  num_sneighs = npts*MAX_NUM_NEIGHS;
  for(i=0; i< num_sneighs; i++) {
    sneighs[i] = -1;
  }

  mindist = 99; minidx = 888;	/* prevent compiler thinking these might
				 * be uninitialised. */
  for (p=0; p<npts; p++) {

    /* printf("finding nearest to point %d\n", p); */
    if (reject[p] && ignore_rejects) { 
      if (nndfp) fputs("-1 -1\n", nndfp);
      if (snfp) fputs("\n", snfp);
      
      continue;		/* jump to next p in the for loop. */
    }

    first_check = 1;
    px = xpts[p];	py = ypts[p];


    /* For each cell, look amongst its neighbours
     * to find its closest neighbour. */
	
    for(n=0; n < numneighs[p]; n++) {
      i = neighs[NIND(p,n)];
      x = xpts[i]; y = ypts[i];
      dx = (x - px); dy = (y - py);
      dist2 = (dx*dx) + (dy*dy);
	    
      if (first_check ||(dist2 < mindist)) {
	mindist = dist2;
	minidx = i;
	first_check = 0;
      }

      if (sort_neighs) {
	if (n >= MAX_DISTS) {
	  printf("%s:%d: too many distances %d\n",
		 __FILE__, __LINE__, n);
	  exit(-1);
	}
	dists[n].key = i + first_index;
	dists[n].dist = dist2;
      }
    }

    dist = (Sfloat)sqrt((double)mindist);
    v = minidx + first_index;
    if (nndfp) fprintf(nndfp, "%.4f %d\n", dist, v);

    info[RIND(p,1,npts)] = v;	/* i.d. of nearest neigh */
    info[RIND(p,2,npts)] = dist; /* distance to nearest neigh */


    if (sort_neighs) {
      /* Sort the neighbours according to distance from data point. */

      qsort (dists, numneighs[p], sizeof (Keydist), keydist_cmp);
      /*  printf("storing %d sneighs of cell %d\n", numneighs[p], p); */
      for(i=0; i<numneighs[p]; i++) {
	if (snfp) fprintf(snfp,"%d ", dists[i].key);
	sneighs[SNIND(p,i,npts)] = dists[i].key;
      }
	  
      if (snfp) fputs("\n", snfp);
    }

  }

  
  /*      ##print "nearest to px py is pts[minidx][0] pts[minidx][1]\n"; */

/*      if (sort_neighs) { */
/*  	## output the sorted neighbour list too. */
/*  	&write_neighs(sorted_neighs_file); */
/*      } */

  if (nndfp) fclose(nndfp);
	
  if (sort_neighs) {
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
#ifdef alloc_vx
  verticeso = (int*)calloc(max_numvertices_o*npts, sizeof(int));
#endif


  if (!vertices1 || !vertices2 || !verticeso) {
    printf("couldn't allocate vertices1/2/o\n");
    exit(-1);
  }

  
  if (! numvertices) { 
    printf("%s: could not allocate space for numvertices\n", "sjevorxx");
    exit(-1);
  }
  
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
      printf("%s:%d exceeded max_numvertices: %d\n",
	     __FILE__, __LINE__, max_numvertices);
      exit(-1);
    }
    vertices1[VIND(p,numv)] = v1;
    vertices2[VIND(p,numv)] = v2;

    p = p2;
    numv = numvertices[p];
    if ((numvertices[p] ++) > max_numvertices) {
      printf("%s:%d exceeded max_numvertices: %d\n",
	     __FILE__, __LINE__, max_numvertices);
      exit(-1);
    }
    vertices1[VIND(p,numv)] = v1;
    vertices2[VIND(p,numv)] = v2;
  }    

  /* test code to print out the vertices belonging to point 0. */
  /*
    for(i=0; i< numvertices[0]; i++) {
    printf("%d %d  ", vertices1[VIND(0,i)], vertices2[VIND(0,i)]);
    }
    printf("\n");
  */
    
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
    /* printf("finding vertices for site %d\n", p);*/
    /* This won't work if the data point has any vertices 
     * out at infinity, so we should ignore them for now. */
    if (!reject[p]) {
      /* Order the vertices for data point p.
       * Go through the list of vertices, finding the next
       * one until we return to the first point. */
	
	  
      first = vertices1[VIND(p,0)]; vertices1[VIND(p,0)] = done;
      next = vertices2[VIND(p,0)]; vertices2[VIND(p,0)] = done;
	
      norder = 0;
      verticeso[VOIND(p,norder)] = first; norder++;
      if (norder > max_numvertices_o) {
	printf("%s:%d max_numvertices_o (%d) reached\n",
	       __FILE__, __LINE__, norder);
	exit(-1);
      }

      looking = 1;
      nverts = numvertices[p];
      while (looking) {
	/*printf("looking for vertice %d\n", next);*/
	verticeso[VOIND(p,norder)] = next;
	norder++;		/* todo: sje: check this doesn't exceed max. */
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
  }
    
  /* After sorting, no longer need the unordered vertices. */
  myfree(vertices1); myfree(vertices2);
    
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



void find_areas(int npts, Sfloat *info)
{
  /* Find the area of each valid polygon.
   * Print out the ordered vertices for each data point.
   *
   * To get the areas, the vertices must first have been ordered,
   * and stored in verticeso.  This is done by find_vertices().
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
  Sfloat xi, yi, xi1, yi1;
  Sfloat sum, dsum;

  firstv = 0;			/* keep compiler quiet. */
  
  for (p=0; p < npts; p++) {
    if (reject[p]) {
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
    info[RIND(p,3,npts)] = sum;	/* area of polygon */
  }
  
}



void find_internal_angles(Sfloat *xpts, Sfloat *ypts, int npts, Sfloat *ias)
{
  /* Compute the internal angles from each site to its vertices. */

  /* If a site has a voronoi polygon with N vertices, we can make N
   * triangles by taking any two contiguous vertices with the site.
   * For each triangle, we measure the internal angle, theta, opposite
   * the polygon edge.  The sum of these angles for a site is 360 degrees.
   * Theta is an angle from a triangle, so is bounded in range (0,180).
   *
   * Since we know the coordinates of the three vertices of each
   * triangle, we compute theta using the cosine rule. */
   
  int nextfree = 0;

  int p, n, i, v1, v2;
  Sfloat xi2, yi2, xi1, yi1, x, y, a,b2,c;
  Sfloat theta;

  for (p=0; p < npts; p++) {
    /* Loop over each site. */

    if (reject[p]) {
      continue;
    }

    /* Find coordinates of central site. */
    x = xpts[p]; y = ypts[p];
    
    n = numvertices[p];

    for (i=0; i < n; i++) {
      v1 = verticeso[VOIND(p,i)];
      if (i == (n-1)) {
	v2 = verticeso[VOIND(p,0)];
      } else {
	v2 = verticeso[VOIND(p,i+1)];
      }

      xi1 = vx[v1]; yi1 = vy[v1];
      xi2 = vx[v2]; yi2 = vy[v2];


      /* Find the lengths of each side of the triangle.  Note we only
       * need b^2 rather than b for the side of the triangle opposite theta.
       */
	 
      a  = sqrt( ((x-xi1)*(x-xi1)) +     ((y-yi1)*(y-yi1)) );
      b2 =     ( ((xi2-xi1)*(xi2-xi1)) + ((yi2-yi1)*(yi2-yi1)) );
      c  = sqrt( ((x-xi2)*(x-xi2)) +     ((y-yi2)*(y-yi2)) );
      theta = acos( ((a*a) + (c*c) - b2) / (2*a*c)) * RAD_TO_DEG;

      ias[nextfree++] = theta;
    }
  }
  ias[nextfree++] = -1;		/* terminator mark. */

}

void sjevoradd(Sfloat *xpts, Sfloat *ypts, Sfloat *temp, int npts)
{
  /* Test function to see that matwrap is working okay.
   * temp[i] = xpts[i] + ypts[i] */
  
  int i;
  for (i=0; i< npts; i++)
    temp[i] = xpts[i] + ypts[i];

}
