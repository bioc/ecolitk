\name{wstringapply}
\alias{wstringapply}
\title{ Apply a function on a window sliding on a string }
\description{
  Apply a function on a window sliding on a string.
}
\usage{
wstringapply(x, SIZE, SLIDE, FUN, ...)
}
\arguments{
  \item{x}{ The string }
  \item{SIZE}{ The size of the window (number of characters). }
  \item{SLIDE}{ offset to move at each slide}
  \item{FUN}{ The function to be applied }
  \item{\dots}{ optional parameter for the function \code{FUN} }
}
\details{
  Apply the function \code{FUN} to substrings of \code{x} of length \code{SIZE}.
}
\value{
  A list of size \code{nchar(x) - SIZE}.
}
\author{ L, Gautier }
\examples{
}
\keyword{ manip }

