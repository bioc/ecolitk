cPlotCircle <- function(radius=1, xlim=c(-2, 2), ylim=xlim, edges=300, main=NULL, main.inside=NULL, ...) {
  plot.new()
  plot.window(xlim, ylim, ...)
  linesCircle(radius, edges=edges, ...)
  title(main)
  text(0, 0, main.inside)
}

linesCircle <- function(radius, center.x = 0, center.y = 0, edges=300, ...) {
  xy <- polar2xy(radius, seq(0, 2*pi, length=edges))
  ##x <- (radius * cos(seq(0, 2*pi, length=edges))) + center.x
  ##y <- (radius * sin(seq(0, 2*pi, length=edges))) + center.y
  lines(xy$x + center.x, xy$y + center.y, ...)
}

chromPos2angle <- function(pos, len.chrom, rot=pi/2, clockwise=TRUE) {
  if (any(abs(pos) > len.chrom, na.rm=TRUE))
    warning(paste(pos, ">", len.chrom, ": abs(pos) > len.chrom !!!"))
  
  theta <- pos * 2 * pi / len.chrom
  
  if (clockwise)
    theta <- - theta

  theta <- theta + rot

  return(theta)
}

polygonDisk <- function(radius, center.x=0, center.y=0, edges=300, ...) {
  ##x <- (radius * cos(seq(0, 2*pi, length=edges))) + center.x
  ##y <- (radius * sin(seq(0, 2*pi, length=edges))) + center.y
  
  xy <- polar2xy(radius, seq(0, 2*pi, length=edges))
  polygon(xy$x + center.x, xy$y + center.y, ...)
  
}

# linesPolar <- function(theta, radius, center.x = 0, center.y = 0, ...) {
#   xy <- polar2xy(radius, theta0, theta1, length=edges))
#   xy$x <- xy$x + center.x
#   xy$y <- xy$y + center.y
#   lines(xy$x, xy$y, ...)
# }

pointsArc <- function(theta0, theta1, radius, center.x = 0, center.y = 0, ...) {
  xy <- polar2xy(radius, seq(theta0, theta1, length=length(radius)))
  ##x <- (radius * cos(seq(0, 2*pi, length=edges))) + center.x
  ##y <- (radius * sin(seq(0, 2*pi, length=edges))) + center.y
  points(xy$x + center.x, xy$y + center.y, ...)
}

linesArc <- function(theta0, theta1, radius, center.x = 0, center.y = 0, ...) {
  xy <- polar2xy(radius, seq(theta0, theta1, length=length(radius)))
  ##x <- (radius * cos(seq(0, 2*pi, length=edges))) + center.x
  ##y <- (radius * sin(seq(0, 2*pi, length=edges))) + center.y
  lines(xy$x + center.x, xy$y + center.y, ...)
}

arrowsArc <- function(theta0, theta1, radius, center.x = 0, center.y = 0, edges = 10,
                      length = 0.25, angle = 30, code = 2, ...) {
  xy <- polar2xy(radius, seq(theta0, theta1, length=edges))
  xy$x <- xy$x + center.x
  xy$y <- xy$y + center.y
  lines(xy$x, xy$y, ...)
  n <- length(xy$x)
  if (code == 2 | code == 3)
    arrows(xy$x[n-1], xy$y[n-1], xy$x[n], xy$y[n], length=length, angle=angle, ...)
  if (code == 1 | code == 3)
    arrows(xy$x[1], xy$y[1], xy$x[2], xy$y[2], length=length, angle=angle, ...)
  
}

polygonArc <- function(theta0, theta1, radius.in, radius.out,
                       center.x = 0, center.y = 0,
                       edges=10,
                       col="black",
                       border = NA,
                       ...) {
  
  
  if (length(edges) == 1)
    edges <- rep(edges, length=length(theta0))

  col <- rep(col, length = length(theta0))

  ok <- ! (is.na(theta0) | is.na(theta1))
  
  for (i in seq(along=theta0[ok])) {
    theta.seq <- seq(theta0[ok][i], theta1[ok][i], length=edges[ok][i])
    x <- c(radius.in * cos(theta.seq), radius.out * cos(rev(theta.seq))) + center.x 
    y <- c(radius.in * sin(theta.seq), radius.out * sin(rev(theta.seq))) + center.y
    polygon(x, y, col=col[ok][i], border=border, ...)
  }
}

polygonChrom <- function(begin, end, len.chrom,
                         radius.in, radius.out,
                         total.edges=300,
                         edges=max(round(abs(end-begin)/len.chrom* total.edges), 2, na.rm=TRUE),
                         rot=pi/2, clockwise=TRUE,
                         ...) {

  theta0 <- chromPos2angle(begin, len.chrom, rot=rot, clockwise=clockwise)
  theta1 <- chromPos2angle(end, len.chrom, rot=rot, clockwise=clockwise)

  if (any(theta0 == theta1, na.rm=TRUE))
    warning(paste("identical angles for: ", which(theta0 == theta1), collapse=TRUE))
  
  polygonArc(theta0, theta1, radius.in, radius.out, edges=edges, ...)
}

linesChrom <- function(begin, end, len.chrom, radius,
                       total.edges=300,
                       edges=max(round(abs(end-begin)/len.chrom* total.edges), 2, na.rm=TRUE),
                       rot=pi/2, clockwise=TRUE,
                       ...) {

  theta0 <- chromPos2angle(begin, len.chrom, rot=rot, clockwise=clockwise)
  theta1 <- chromPos2angle(end, len.chrom, rot=rot, clockwise=clockwise)

  if (any(theta0 == theta1, na.rm=TRUE))
    warning(paste("identical angles for: ", which(theta0 == theta1), collapse=TRUE))

  if (length(edges) == 1)
    edges <- rep(edges, length=length(theta0))

  ok <- ! (is.na(theta0) | is.na(theta1))
  
  for (i in seq(along=theta0[ok])) {
    linesArc(theta0[ok][i], theta1[ok][i], rep(radius, edges[ok][i]), ...)
  }
}

ecoli.len <- 4639221
