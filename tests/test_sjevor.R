## Extra tests for the Voronoi library.
## Tue 12 Mar 2002
library(sjevor)
invisible(.libPaths( Sys.getenv("MYRLIB")))
library(deldir)
print(1)
##library(sjedmin)
print(2)
## This big data file caused problems -- there are 4500 cells, which
## were originally generated from a toroidal sample of 500 cells.
## Until Tue 12 Mar 2002, it was crashing with more than 4000 cells.
## This was because of realloc problems in sje_readsites().

# d <- matrix(data = scan("./nadata_bad.txt"), ncol =2, byrow = T)
# plot(d[,1], d[,2])

# w <- 1:4500
# v <- vorcr(d[w,1], d[w,2], -1000, 2000, -1000, 2000, 'n')
# plot(v)
# vorcr.polygons(v)

# v <- vorcrb(d[w,1], d[w,2], -1000, 2000, -1000, 2000, 'n')
# plot(v)
# vorcr.polygons(v)


## Now let's compare the output from deldir package with the output from
## my package.  Create several dmin mosaics, and compare:
## - number of neighborus
## - Voronoi area.

scale <- 2.2
wid <- 1000*scale; ht <- 1000*scale;
dmin <- 8*scale
dminsd <- 2

##npts <- 1000*(scale^2)
npts <- 100*(scale^2)
for (i in 1:20) {
  cat(paste("i",i,"\n"));
  ##d <- dminsd(wid, ht, npts, dmin, dminsd)
  ## selecting pts at random doesn't seem to work at high densities.
  x <- runif(npts)*wid; y <- runif(npts)*ht; d <- list(x=x, y=y);
  ##plot(d)
  v <- vorcr(d$x, d$y, 0, wid, 0, ht)
  try <- deldir(d$x, d$y, plot=FALSE, wl='tess')
  
  ## check the number of neighs
  dd.numneighs.valid <- try$summary[-which(v$rejects),3]
  my.numneighs.valid <- v$numneighs[-which(v$rejects)]

  ##stopifnot(identical(all.equal(dd.numneighs.valid, my.numneighs.valid),T))
  diffs <- abs(dd.numneighs.valid - my.numneighs.valid)
  print(quantile(diffs))
  
  ## check the areas.
  dd.areas.valid <- try$summary[-which(v$rejects),8]
  my.areas.valid <- v$info[-which(v$rejects),4]
  ##stopifnot(identical(all.equal(dd.areas.valid, my.areas.valid),T))
  diffs <- abs(dd.areas.valid - my.areas.valid)
  print(quantile(diffs))

  ##plot(v)
  table(v$numneighs)
}
