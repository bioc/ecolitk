\name{ecoligenomeBNUM}
\alias{ecoligenomeBNUM}
\alias{ecoligenomeBNUM2SYMBOL}
\alias{ecoligenomeBNUM2ENZYME}
\alias{ecoligenomeBNUM2GENETYPE}
\alias{ecoligenomeBNUM2GENBANK}
\alias{ecoligenomeBNUM2GENEPRODUCT}
\docType{data}
\title{ Environment for 'bnum' identifiers}
\description{
  An environment to associate Affymetrix probe set IDs with 'bnum' IDs
}
\usage{
data(ecoligenomeBNUM)
data(ecoligenomeBNUM2SYMBOL)
data(ecoligenomeBNUM2ENZYME)
data(ecoligenomeBNUM2GENETYPE)
data(ecoligenomeBNUM2GENBANK)
data(ecoligenomeBNUM2GENEPRODUCT)
}
\format{
  These are \code{environment} objects.
}
\details{
  Escherichia coli genes are sometimes identified by 'bnum's. This
  identfier is typically a 'b' followed by digits.  
}
\source{
  BNUM numbers were parsed out of the Affymetrix identifiers.
  BNUM2* were obtained from the GenProtEC website.
}
\keyword{datasets}