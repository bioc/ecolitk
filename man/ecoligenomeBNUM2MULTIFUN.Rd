\name{ecoligenomeBNUM2MULTIFUN}
\alias{ecoligenomeBNUM2MULTIFUN}
\alias{ecoligenomeBNUM2STRAND}
\docType{data}
\title{ Environment }
\description{
  An environment to store associtations between 'bnum' identifiers (key)
  and 'MultiFun' identifiers (or strand information).
}
\usage{data(ecoligenomeBNUM2MULTIFUN)}
\format{
  The format is:
length 0 <environment>
 - attr(*, "comments")= chr "GenProtEC: MultiFun assignments for E. coli modules
September 17th, 2003"
}
\details{
  'MultiFun' is a classification scheme. The structure is
  'approximately tree-like'.
  Several 'MultiFun' numbers can be assigned to one 'bnum'.
}
\source{
  "http://genprotec.mbl.edu/files/MultiFun.txt"
}
\keyword{datasets}
