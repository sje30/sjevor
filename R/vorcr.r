## Do vorcr in R.
## Wed 07 Jun 2000
## Can be loaded as a source file, as in:
## source("/home/stephen/langs/R/vorr/vorcr.r")

vorcr <- function(x, y, xl, xh, yl, yh, fuzz = 0, opts = 'nags') {
    ## Do Voronoi analysis and then return various useful bits of info.
    
    dms <- c(xl, xh, yl, yh)
    npts <- length(x)

    ## maximum number of polygon npts to allow.  for each line, there
    ## will be four values (x1,y1 x2, y2), so we have five times the
    ## number of lines to be safe.
    max.num.neighs <- 18
    max.del.tris <- 5           #see also sjevor.h -- these values should agree

    polynpts1 <- integer(1); polynpts1[1] <- npts * 4 * 5
    z <- .C("sjevor",
            as.double(x), as.double(y),
            as.double(dms),
            as.character(opts),
            info = double(npts*4),
            sneighs = integer(npts*max.num.neighs),
            iangles = double(npts*10),
            ## make space for 5n del tris; each triangle needs 3 items.
            delids  = integer(npts*max.del.tris*3),
            dellens = double(npts*max.del.tris*3),
            delangs = double(npts*max.del.tris*3),
            polypts = double(polynpts1[1]), #4 for each line: (x1,y1, x2,y2)
            polynpts = as.integer(polynpts1),
            vertices.x = double(max.del.tris*npts),
            vertices.y = double(max.del.tris*npts),
            vertices = integer(npts*max.num.neighs),
            as.integer(npts)
            )

    info <- z$info; dim(info) <- c(npts,4)
    colnames(info) <- c("id", "nn id", "dist", "area")


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

    polypts <- z$polypts[1:z$polynpts]
    vertices.xy <- cbind(z$vertices.x, z$vertices.y)


    ## Normally ignore.rejects is true so that we reject triangles
    ## that involve reject sites.
    ignore.rejects <- TRUE;
    if (ignore.rejects) {
      anyrej <- rejects[delids]; dim(anyrej) <- c(length(anyrej)/3,3);
      ## rejtri[i] is true if the ith triangle should be rejected.
      rejtri <- apply(anyrej, 1, any)
    } else {
      ## accept all delauanay triangles.
      ## check that delidmax below is the right thing to do...
      rejtri <- logical(length = (delidmax/3))# all elements FALSE.
    }

    delacc <- which(!rejtri)            #ids of accepted triangles.
    delrej <- which(rejtri)             #ids of rejected triangles.
    
    res <- list(info = info,
                pts = cbind(x, y),
                neighs = sneighs, cr = cr, meannnd = meannnd,
                sdnnd = sdnnd, rejects = rejects, iangles = iangles,
                delids = delids, dellens = dellens, delangs = delangs,
                delacc = delacc, delrej = delrej, polypts = polypts,
                numneighs = numneighs,
                vertices.xy = vertices.xy,
                vertices= matrix(z$vertices, nrow=npts, byrow=T))

    class(res) <- "sjevor"
    res
}

vorcrb <- function(x, y, xl, xh, yl, yh, fuzz = 0, f=0.3, opts = 'nags') {
  ## Toroidal boundary conditions...
  ## This format was copied from vorb_c.m
  ## TODO: del.plot doesn't work..., but on the other hand, the vorcr.polygons
  ## routine looks very nice...!
  
  ht <- yh - yl; wid <- xh - xl;
  ncells <- length(x)

  ##pdy = positive dy, ndy = negative dy.
  pdy  <- ht + numeric(ncells); ndy <- -pdy;

  pdx  <- wid + numeric(ncells); ndx  <- -pdx;
  emp  <- numeric(ncells)
  data <- cbind(x,y);			# original data points.
  data1 <- data + cbind(emp, pdy);
  data2 <- data + cbind(emp, ndy);
  
  data3 <- data + cbind(ndx, emp);
  data4 <- data + cbind(pdx, emp);
  
  data5 <- data + cbind(ndx, pdy);
  data6 <- data + cbind(pdx, pdy);
  
  data7 <- data + cbind(ndx, ndy);
  data8 <- data + cbind(pdx, ndy);
  
  ndata <- rbind(data, data1, data2, data3, data4, data5, data6, data7, data8)

  fxl <- xl - (f*wid); fyl <- yl - (f*ht);
  fxh <- xh + (f*wid); fyh <- yh + (f*ht);

  ids <- c(1:ncells, 1:ncells, 1:ncells, 1:ncells, 1:ncells, 1:ncells, 1:ncells,
           1:ncells, 1:ncells)

  pass <- ((ndata[,1] >= fxl) & (ndata[,1] <= fxh) & 
           (ndata[,2] >= fyl) & (ndata[,2] <= fyh));
  ndata <- ndata[which(pass),];
  ids   <- ids[which(pass)];	# keep only cell IDs that are within border.

  v <- vorcr(ndata[,1], ndata[,2], fxl, fxh, fyl, fyh, fuzz, opts)

  ## now post-process v a bi, by removing the information for cells that we
  ## don't care about.
  v$info <- v$info[1:ncells,]
  v$info[,2] <- ids[v$info[,2]]
  
  v$pts <- v$pts[1:ncells,]

  v$neighs <- v$neighs[1:ncells,]
  excess <- which(v$neighs>ncells, arr.ind=T)
  v$neighs[excess] <- ids[v$neighs[excess]]

  v$rejects <- v$rejects[1:ncells]
  v$numneighs <- v$numneighs[1:ncells]
  ## little sanity check -- nearest neighbours should be the first neighbour.
  stopifnot(v$info[,2] == v$neighs[,1])

  ## TODO: delaunay triangles do not yet plot.
  v
}


vor.showpts <- function(v, n=10, max.neighs=NULL) {
  ## View the Voronoi tesselation, and click on cells to show their neighbours.
  ## N is the number of cells to select.
  ## MAX.NEIGHS, if a number, is the number of neighbours to show.  e.g. set to 1
  ## to show just the nearest neighbours.
  first <- TRUE
  for (i in 1:n) {
    if (first) {
      plot(v); first <- FALSE
    }
    id <- identify(v$pts[,1], v$pts[,2], n=1)
    neighs <- v$neighs[id,]
    neighs <- neighs[which(neighs>0)]
    if (is.numeric(max.neighs) && (length(neighs)> max.neighs))
      neighs <- neighs[1:max.neighs]
    points(v$pts[neighs,1], v$pts[neighs,2], pch=19)
  }
}

vorcr.dellens <- function(vor, idxs=NULL) {
  ## Helper function to get the Delaunay Segment lengths.
  ##
  ## IDXS is a vector of triangle indexes for which we want to compute
  ## segment lengths.  If this is NULL, we should compute the
  ## segment lengths for all triangles.

  if (length(idxs) == 0)
    idxs <- 1:dim(vor$delids)[1]
  nsites <- dim(vor$info)[1]

  ## ds will be a sparse, upper triangular matrix.
  ds <- matrix(0, nrow=nsites, ncol = nsites)

  ## Use sort() to ensure that first index is always smaller than the
  ## second index, to make an upper triangular matrix.
  ## e.g. t(apply(cbind( c(1,7,8), c(5,2,6)),1,sort))

  ds[ t(apply(cbind(vor$delids[idxs,1], vor$delids[idxs,2]),1,sort)) ] <-
    vor$dellens[idxs,1]
  ds[ t(apply(cbind(vor$delids[idxs,2], vor$delids[idxs,3]),1,sort)) ] <-
    vor$dellens[idxs,2]
  ds[ t(apply(cbind(vor$delids[idxs,3], vor$delids[idxs,1]),1,sort)) ] <-
    vor$dellens[idxs,3]
  ds[which(ds>0)]
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

del.plot <- function(v) {
  ## Plot the Delaunay triangulation.
  ## pts is the 2d set of data points; v is the voronoi information from
  ## that data set.

  pts <- v$pts
  ## First draw the sites.
  par(col ="black")
  par(mfrow=c(1,1), pty="s")
  plot(pts, type="n")                   #don't plot points, just set ranges.
  text(pts[,1], pts[,2], seq(1:length(pts))) #label the points.

  
  ## draw rejected tris.
  par(col ="red")
  for (i in v$delrej) {
    t <- c(v$delids[i,], v$delids[i,1]);
    lines(pts[t,])
  }
  
  ## draw accepted tris.
  par(col ="blue")
  for (i in v$delacc) {
    t <- c(v$delids[i,], v$delids[i,1]);
    lines(pts[t,])
  }

  par(col = "black")                    #reset to usual colour.
}


plot.sjevor <- function(v, show.pts=T, show.areas=F, show.rejects=F, ...) {

  ## line-based approach to doing the plot.  We take the vector
  ## v$polypts; each set of 4 consecutive values is (x0, y0, x1, y1)
  ## -- for a line to be drawn from (x0,y0) to (x1,y1).
  
  np <- length(v$polypts)
  x0 <- v$polypts[seq(from=1, to=np, by=4)]
  y0 <- v$polypts[seq(from=2, to=np, by=4)]
  x1 <- v$polypts[seq(from=3, to=np, by=4)]
  y1 <- v$polypts[seq(from=4, to=np, by=4)]

  if (show.areas) {
    plot(v$pts[,1], v$pts[,2], type="n", asp=1, xlab="", ylab="")
    text(v$pts[,1], v$pts[,2], signif(v$info[,4],3))
  } else {
    if (show.pts) {
      if (show.rejects) {
        plot(v$pts[,1], pch=1, v$pts[,2], asp=1, xlab="", ylab="")
        rejects <- which(v$info[,4] < 0)
        points(v$pts[rejects,1], v$pts[rejects,2], pch=19, asp=1)
      } else {
        plot(v$pts[,1], v$pts[,2], asp=1, xlab="", ylab="", ...)
      }
    } else {
      ## don't want to see points
      plot(v$pts[,1], v$pts[,2], type="n",asp=1, xlab="", ylab="")
    }
  }
  segments(x0, y0, x1, y1)
}

vorcr.polygons <- function(v) {
  ## For each site, find its surrounding polygon and plot it.  This is
  ## a slower method than vor.plot(), but it shows how to use the vertice
  ## information associated with each site.

  ## Check that the vertices were created.
  if (!any(v$vertices > 0)) {
    stop("no vertice information for each site was stored. Use option a")
  }
  pts <- v$pts
  npts <- dim(pts)[1]
  plot(pts)
  for (i in 1:npts) {
    nverts <- v$numneighs[i]
    if (nverts >0) {
      polygon(v$vertices.xy[v$vertices[i,1:nverts],])
    }
  }
}
