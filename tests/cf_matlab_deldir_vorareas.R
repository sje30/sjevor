## Check my sjevor library against the alternatives (matlab, deldir in R)
## Wed 26 Feb 2003

width <- 3000; height <- 3000;
npts <- 1000


## If this is TRUE, we regenerate some dmins and run matlab again...
need.newdmins <- FALSE

dmins <- c(6, 10, 20, 40, 45, 50)
sds <- c(1, 4, 4, 4, 2, 4)
nreps <- length(dmins)

if (need.newdmins) {
  ## Generate some random cell distributions, using dmin algorithm.
  library(sjedmin)
  
  for (i in 1:nreps) {
    dat <- dminlul(width, height, npts, dmin=dmins[i], dminsd=sds[i])
    write.table(signif(cbind(dat$x, dat$y),7),
                paste("dmin",i,".txt",sep=""),
                quote=F,col.names=F, row.names=F)
  }

  ## Now run matlab if need be.
  system("matlab -nojvm < matlab_vorareas.m")
}

## Otherwise, don't worry about making new dmins...
library(sjevor)
library(deldir)                         #used for comparison.



compare <- function(infile, matfile) {

  my.quantiles <- function(x, y) {
    ## Compare two distributions
    print(signif(quantile(abs(v.areas - mat.areas), na.rm=T,
                          probs=c(0, .05, .25, .50, 0.75, 0.95, 1)),4))
  }

  
  mat.areas <- scan(matfile, quiet=T)
  dat <- read.table(infile)
  v <- vorcr(dat[,1], dat[,2], 0, width, 0, height, fuzz=0)
  v.areas <- v$info[,4]
  v.areas[which(v.areas<0)] <- NaN
  
  cat("sjevor vs matlab\n")
  my.quantiles(v.areas, mat.areas)
  ##
  ## compare with deldir.  deldir has no NaN numbers, so should wipe them out.
  try <- deldir(dat[,1], dat[,2], plot=F, wl='tess')
  deldir.areas <- try$summary[,8]
  deldir.areas[which(is.nan(v.areas))] <- NaN
  cat("sjevor vs deldir\n")
  my.quantiles(v.areas, deldir.areas)
  ##
  ##
  cat("deldir vs matlab\n")
  my.quantiles( deldir.areas, mat.areas)
}

##par(mfcol=c(2,nreps))
for (i in 1:nreps) {
  print(i)
  infile <- gsub("XX", as.character(i), "dminXX.txt")
  matfile <- gsub("XX", as.character(i), "dminXX_matareas.txt")
  compare(infile, matfile)
}

compare("w81s.on.d", "w81s.on.matareas.txt")
compare("triarray.dat", "triarray_matareas.txt")
## most areas seem about on order of 1e-5.  Typical errors range
## between 10-8 and 10-4.  Small differences may be due to numerical
## fuzz added (e.g. "frac" used by deldir; I think matlab 6 no longer
## uses fuzz.)
