## Do vorcr in R.
## Wed 07 Jun 2000
## Can be loaded as a source file, as in:
## source("/home/stephen/langs/R/vorr/vorcr.r")


dyn.load("/home/stephen/langs/R/vorr/libvor.so")


vorcr <- function(x, y, xl, xh, yl, yh, fuzz = 0, opts = 'nags') {
    ## Do Voronoi analysis and then return various useful bits of info.
    
    dms <- c(xl, xh, yl, yh)
    npts <- length(x)
    z <- .C("sjevor",
            as.double(x), as.double(y),
            as.double(dms),
            as.character(opts),
            info = double(npts*4),
            sneighs = integer(npts*18),
            iangles = double(npts*10),

            ## make space for 5n del tris; each triangle needs 3 items.
            delids  = integer(npts*15),
            dellens  = double(npts*15),
            delangs = double(npts*15),
            as.integer(npts)
            )
    info <- z$info; dim(info) <- c(npts,4)


    sneighs <- z$sneighs;
    maxnumneighs <- length(sneighs) / npts;
    dim(sneighs) <- c(npts,maxnumneighs)
    ## find the number of neighbours of each site.
    numneighs <- apply(sneighs,1, function (x) length(x[x>0]))

    ## can now shorten the list of internal angles, since their should
    ## be the same as sum(numneighs).
    iangles <- z$iangles[1:sum(numneighs)]
    rejects <- (info[,3] < 0)
    validdists <- info[!rejects,3]
    meannnd <- mean(validdists)
    sdnnd   <- sqrt(var(validdists))
    cr <-  meannnd / sdnnd


    ## Process Delaunay triangle information.
    delidmax <- which(z$delids == -1) -1

    ## Truncate Delaunay info to right length.
    delids  <- z$delids[1:delidmax]
    dim(delids) <- c(3, delidmax/3); delids <- t(delids);

    dellens <- z$dellens[1:delidmax]
    dim(dellens) <- c(3, delidmax/3); dellens <- t(dellens);

    delangs <- z$delangs[1:delidmax]
    dim(delangs) <- c(3, delidmax/3); delangs <- t(delangs);

    list(info = info, neighs = sneighs, cr = cr, meannnd = meannnd,
         sdnnd = sdnnd, rejects = rejects, iangles = iangles,
         delids = delids, dellens = dellens, delangs = delangs,
         numneighs = numneighs)
}




ianglesplot <- function(angles, show=TRUE)  {
  ## Produce a histogram of the internal angles of each Voronoi polygon.
  abins <- seq(0,180, by=5)
  ah <- hist(angles, breaks=abins, plot=F)
  
  cdf <- cumsum(ah$counts)/sum(ah$counts)

  if (show) {
    par(lab = c(10,10,0))
    par(xaxp = c(0,18,1))   ## Can I get xaxp to work --can I heck!?!?
    plot(ah$mids,cumsum(ah$counts)/sum(ah$counts), type="l",
         xlab="angle (deg)", ylab="cumulative probability")
    title("interior angles")
  }

  list(x=ah$mids, y=cdf)
}

## dyn.unload("/home/stephen/langs/R/vorr/libvor.so")
