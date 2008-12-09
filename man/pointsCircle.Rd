\name{pointsCircle}
\alias{linesCircle}
\alias{polygonDisk}
\alias{arrowsArc}
\alias{pointsArc}
\alias{linesArc}
\alias{polygonArc}
\title{ Functions to plot circular related figures }
\description{
  Functions to plot circular related figures
}
\usage{
linesCircle(radius, center.x = 0, center.y = 0, edges = 300, ...)
polygonDisk(radius, center.x = 0, center.y = 0, edges=300,
...)
arrowsArc(theta0, theta1, radius, center.x = 0, center.y = 0, edges = 10,
          length = 0.25, angle = 30, code = 2, ...)
pointsArc(theta0, theta1, radius, center.x = 0, center.y = 0, ...)
linesArc(theta0, theta1, radius, center.x = 0, center.y = 0, ...)
polygonArc(theta0, theta1, radius.in, radius.out,
           center.x = 0, center.y = 0,
           edges = 10,
           col = "black",
           border = NA,
           ...)
}
\arguments{
  \item{theta0, theta1}{ start and end angles for the arc}
  \item{radius}{ radius of the circle }
  \item{radius.in}{ inner radius }
  \item{radius.out}{ outer radius }
  \item{center.x, center.y}{Coordinates for the center of the circle
    (default to (0, 0))}
  \item{edges}{ number of edges the shape is made of }
  \item{col}{color}
  \item{border}{border (see \code{\link[base]{polygon}})}
  \item{length, angle, code}{see the corresponding parameters for the
    function \code{\link{arrows}}}
  \item{\dots}{ optional graphical paramaters}
}
\details{
  Details to come\ldots for now the best to run the examples and experiment by
  yourself\ldots
%%  The functions \code{linesArc}, \code{pointsArc}, \code{arrowsArc} 
}
\value{
  Function only used for their border effects.
}
\author{ laurent }
\examples{
par(mfrow=c(2,2))
n <- 10
thetas <- rev(seq(0, 2 * pi, length=n))

rhos <- rev(seq(1, n) / n)

xy <- polar2xy(rhos, thetas)
colo <- heat.colors(n)

plot(0, 0, xlim=c(-2, 2), ylim=c(-2, 2), type="n")
for (i in 1:n)
  linesCircle(rhos[i]/2, xy$x[i], xy$y[i])

plot(0, 0, xlim=c(-2, 2), ylim=c(-2, 2), type="n")
for (i in 1:n)
  polygonDisk(rhos[i]/2, xy$x[i], xy$y[i], col=colo[i])

plot(0, 0, xlim=c(-2, 2), ylim=c(-2, 2), type="n", xlab="", ylab="")
for (i in 1:n)
  polygonArc(0, thetas[i],
             rhos[i]/2, rhos[i],
             center.x = xy$x[i], center.y = xy$y[i], col=colo[i])

plot(0, 0, xlim=c(-2, 2), ylim=c(-2, 2), type="n", xlab="", ylab="")
for (i in (1:n)[-1]) {
  linesCircle(rhos[i-1], col="gray", lty=2)
  polygonArc(thetas[i-1], thetas[i],
             rhos[i-1], rhos[i], col=colo[i],
             edges=20)
  arrowsArc(thetas[i-1], thetas[i],
             rhos[i] + 1, col=colo[i],
             edges=20)
}
  
}
\keyword{ aplot }
