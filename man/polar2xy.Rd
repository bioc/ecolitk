\name{polar2xy}
\alias{polar2xy}
\alias{xy2polar}
\alias{rotate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Functions to perform polar coordinate related functions}
\description{
  Functions to perform polar coordinate related functions
}
\usage{
polar2xy(rho, theta)
xy2polar(x, y)
rotate(x, y, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ cartesian coordinate }
  \item{y}{ cartesian coordinate  }
  \item{rho}{ polar radius \code{rho} }
  \item{theta}{ polar angle \code{theta}}
  \item{alpha}{ angle to perform rotation }
}
\details{
  \code{y} and \code{theta} can be respectively missing. In this case,
  \code{x} and \code{rho} are expected to be lists with entries
  \code{x, y}, \code{rho, theta} respectively.
}
\value{
  %Lists of
  %\item x,y
  %or
  %\item rho, theta.
}
\examples{


n <- 40
nn <- 2
thetas <- seq(0, nn * 2 * pi, length=n)

rhos <- seq(1, n) / n

plot(c(-1, 1), c(-1, 1), type="n")
abline(h=0, col="grey")
abline(v=0, col="grey")

xy <- polar2xy(rhos, thetas)

points(xy$x, xy$y, col=rainbow(n))

}
\keyword{ manip}
