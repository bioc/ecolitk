\name{gccontent}
\alias{gccontent}
\title{ function to compute gccontent }
\description{
  A simple R function to compute the GC content of a sequence
}
\usage{
gccontent(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a vector of mode \code{character} }
}
\details{
  This a simple (and not particularly fast) function to compute the GC
  content of sequence. When speed is an issue, one should use the
  function in the package \code{matchprobes}. This function only exists
  to avoid dependency on this package.
}
\value{
  The GC content (\code{numeric})
}
%\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\keyword{ manip }% at least one, from doc/KEYWORDS
