######################################################################
## Test functions for the 
## Example: test on Fortune's data set.

library(sjevor)
data(fortune.eg)

##tdata <- scan('/home/stephen/mosaics/mouse-horizontal/data/r4')
tdata <- as.matrix(fortune.eg)

npts <- length(tdata)/2                 #assume data is 2d.
dim(tdata) <- c(2,npts); tdata <- t(tdata)

par(mfrow = c(1,1), pty="s")
plot(tdata)


xmin <- 0; xmax <- max(tdata[,1])
ymin <- 0; ymax <- max(tdata[,2])
v <- vorcr(tdata[,1], tdata[,2], xmin, xmax,  ymin, ymax, opts='snag')

del.plot(tdata, v)
vor.plot(tdata, v)                      #new Thu 27 Dec 2001...
title("vorcr_test_fort.ps: Fortune data set using del.plot()")

v$cr


## Plot histograms of lengths and angles.  These distributions are
## biased by the rejected segments lengths and angles...  So these are
## worth ignoring.

par(mfrow=c(2,3));

hist(v$dellens,            main="all delaunay seg lengths");
hist(v$dellens[v$delacc,], main="accepted delaunay seg lengths");
hist(v$dellens[v$delrej,], main="rejected delaunay seg lengths");

hist(v$delangs,            main="all delaunay tri angles");
hist(v$delangs[v$delacc,], main="accepted delaunay tri angles");
hist(v$delangs[v$delrej,], main="rejected delaunay tri angles");





## To test the Delaunay triangulation, we should be able to make the
## list of nearest neighbours from the delaunay triangulation.

ntris <- length(v$delids)/3             #number of triangles.
neighs <- matrix(data = 0,nrow = length(tdata[,1]), ncol = 30)
nns <- vector(mode = "numeric", length = length(tdata[,1]))
ds <- v$delids; ##dim(ds) <- c(3,ntris); ds <- t(ds);

## Neighs[i,] will store the neighbours of site i;
## nns[i] stores the number of neighbours of site i
## ds stores the site id numbers.

for (i in 1:ntris) {

  ## For each triangle, we find the three sites, a,b,c.  We then note
  ## that b and c are neighbours of a, and likewise for b,c.
  a <- ds[i,1]
  b <- ds[i,2]
  c <- ds[i,3]

  neighs[a,nns[a]+1] <- b
  neighs[a,nns[a]+2] <- c
  nns[a] <- nns[a] + 2

  neighs[b,nns[b]+1] <- a
  neighs[b,nns[b]+2] <- c
  nns[b] <- nns[b] + 2

  neighs[c,nns[c]+1] <- a
  neighs[c,nns[c]+2] <- b
  nns[c] <- nns[c] + 2
}


## Now compare the neighbour list from the voronoi representation with
## the Delaunay representation.
for (i in 1:npts) {
  if (!(v$rejects[i]) ) {
    a <- sort(unique(v$neighs[i,]))
    b <- sort(unique(neighs[i,]))
    a <-a[2:length(a)]                  #ignore first element; will be -1
    b <-b[2:length(b)]                  #ignore first element; will be 0
    s <- sum(abs(a-b))
    if (s > 0) 
      warning(cat(paste(i, s, "\n")))
  }
}
    
  
  
## To test the internal angles, we can use a distorted triangular lattice
## since we expect most of the angles to be around 60 degrees.  As we
## add more noise, the spread of angles should get larger and larger.

trigrid <- function (nx, ny, d)
  ## Generate a triangular grid, with nx pts in the x direction,
  ## and ny in the y direction; d is the separation distance
  ## between points.
  {
    odd <- F
    y <- 0
    pts <- array(0, dim=c(nx*ny,2))
    n <- 1
    for (j in 1:ny) {
      odd <- !odd;
      if (odd)
        x <- 0
      else
        x <- d/2
      for (i in 1:nx) {
        x <- x + d;
        pts[n,1] <- x
        pts[n,2] <- y
        n <- n + 1
      }
      y <- y + (d*sin(pi/3))
    }
    pts
  }


p <- trigrid(10,10,20)
p <- p + (6.0 * runif(nrow(p)))         #Add noise, between 1.0 and 30.0
plot(p)

v <- vorcr(p[,1], p[,2], 0,250,  0,200, opts='snag')
v$cr
r <- ianglesplot(v$iangles)


######################################################################
## Compare the Fortuen code with tripack().

library(tripack)
wid<- 500; ht <- 300; npts <- 200

par(mfrow=c(4,3))
for (i in 1:10) {
  xs <- runif(npts)*wid; ys <- runif(npts)*ht
  v <- vorcr(xs, ys, 0, wid, 0, ht)
  print(v$cr)
  v.areas <- v$info[,4]
  rejects <- which(v$info[,2]==-1)

  ##plot(tripack.vm)
  vor.plot(cbind(xs, ys), v, show.rejects=T)
  tripack.vm <- voronoi.mosaic(xs, ys)
  tripack.vm.areas <- voronoi.area(tripack.vm)
  ##rejects <- voronoi.findrejectsites(tripack.vm, 0, wid, 0, ht)

  ## Check that the areas computed by Voronoi and Fortune are the same.
  v.areas.acc <- v.areas[which(!rejects)]
  tripack.areas.acc <- tripack.vm.areas[which(!rejects)]
  stopifnot(all.equal(tripack.areas.acc, v.areas.acc))
}



## Compare the area of Fortune data set between my code and Matlab.
## matlab areas computed by matlab_vorarea_fortune.m.

data(fortune.eg)
v <- vorcr(fortune.eg[,1], fortune.eg[,2], 0,1, 0,1)
v$cr
v.areas <- v$info[,4]
accepts <- (v$info[,2]>0)
m.areas <- scan(file = "~/langs/R/sjevor/extras/matlab_fortune_areas.dat")
diffs <- v.areas[accepts] - m.areas[accepts]
stopifnot( (diffs < 1e-9))              #all diffs must be less than 1e-9
