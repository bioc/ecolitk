\name{polygonChrom}
\alias{polygonChrom}
\alias{linesChrom}
\alias{cPlotCircle}
\alias{chromPos2angle}
\title{ Functions to plot circular chromosomes informations }
\description{
  Functions to plot circular chromosomes informations
}
\usage{
cPlotCircle(radius=1, xlim=c(-2, 2), ylim=xlim, edges=300, main=NULL,
            main.inside, ...)

chromPos2angle(pos, len.chrom, rot=pi/2, clockwise=TRUE)

polygonChrom(begin, end, len.chrom, radius.in, radius.out,
             total.edges = 300,
             edges = max(round(abs(end - begin)/len.chrom *
                     total.edges), 2, na.rm = TRUE),
             rot = pi/2, clockwise = TRUE, ...)

linesChrom(begin, end, len.chrom, radius,
             total.edges = 300,
             edges = max(round(abs(end - begin)/len.chrom *
                     total.edges), 2, na.rm = TRUE),
             rot = pi/2, clockwise = TRUE, ...)

ecoli.len
}

\arguments{
  \item{radius}{radius}
  \item{xlim, ylim}{ range for the plot. Can be used to zoom-in a
    particular region.}
  \item{pos}{position (nucleic base coordinate)}
  \item{begin}{ begining of the segment (nucleic base number). }
  \item{end}{ end of the segment (nucleic base number). }
  \item{len.chrom}{ length of the chromosome in base pairs }
  \item{radius.in}{ inner radius }
  \item{radius.out}{ outer radius }
  \item{total.edges}{ total number of edges for the chromosome}
  \item{edges}{ number of edges for the specific segment(s) }
  \item{rot}{ rotation (default is \code{pi / 2}, bringing the angle
    zero at 12 o'clock)}
  \item{clockwise}{ rotate clockwise. Default to \code{TRUE}.  }
  \item{main, main.inside}{main titles for the plot}
  \item{\dots}{ optional graphical parameters }
}
\details{
  The function \code{chromPos2angle} is a convenience function.
  The variable ecoli.len contains the size of the Escheria coli genome
  considered (K12).
}
\value{
  Except \code{chromPos2angle}, the function are solely used for their
  border effects.
}
\author{ laurent <laurent@cbs.dtu.dk>}
\seealso{  }
\examples{

data(ecoligenomeSYMBOL2AFFY)
data(ecoligenomeCHRLOC)

## find the operon lactose ("lac*" genes)
lac.i <- grep("^lac", ls(ecoligenomeSYMBOL2AFFY))
lac.symbol <- ls(ecoligenomeSYMBOL2AFFY)[lac.i]
lac.affy <- unlist(lapply(lac.symbol, get, envir=ecoligenomeSYMBOL2AFFY))

beg.end <- lapply(lac.affy, get, envir=ecoligenomeCHRLOC)
beg.end <- matrix(unlist(beg.end), nc=2, byrow=TRUE)

lac.o <- order(beg.end[, 1])

lac.i <- lac.i[lac.o]
lac.symbol <- lac.symbol[lac.o]
lac.affy <- lac.affy[lac.o]
beg.end <- beg.end[lac.o, ]

lac.col <- rainbow(length(lac.affy))

par(mfrow=c(2,2))

## plot

cPlotCircle(main="lac genes")
polygonChrom(beg.end[, 1], beg.end[, 2], ecoli.len, 1, 1.2, col=lac.col)
rect(0, 0, 1.1, 1.1, border="red")

cPlotCircle(xlim=c(0, 1.2), ylim=c(0, 1.1))
polygonChrom(beg.end[, 1], beg.end[, 2], ecoli.len, 1, 1.1, col=lac.col)
rect(0.4, 0.8, 0.7, 1.1, border="red")

cPlotCircle(xlim=c(.45, .5), ylim=c(.85, 1.0))
polygonChrom(beg.end[, 1], beg.end[, 2], ecoli.len, 1, 1.03, col=lac.col)

mid.genes <- apply(beg.end, 1, mean)
mid.angles <- chromPos2angle(mid.genes, ecoli.len)
xy <- polar2xy(1.03, mid.angles)
xy.labels <- data.frame(x = seq(0.45, 0.5, length=4), y = seq(0.95, 1.0, length=4))
segments(xy$x, xy$y, xy.labels$x, xy.labels$y, col=lac.col)
text(xy.labels$x, xy.labels$y, lac.symbol, col=lac.col)

}
\keyword{ hplot}
\keyword{ dplot}
