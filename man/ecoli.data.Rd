\name{ecoli.m52.genome}
\alias{ecoli.m52.genome}
\alias{ecoligenomeCHRLOC}
\alias{ecoligenomeSYMBOL2AFFY}
\alias{ecoligenomeSYMBOL}
\alias{ecoligenomeSTRAND}
\alias{ecoli.len}
%\docType{data}
\title{ Escherichia coli data}
\description{
  Meta-data related to Escherichia coli
}
\usage{
data(ecoli.m52.genome)
data(ecoligenomeCHRLOC)
data(ecoligenomeSYMBOL2AFFY)
data(ecoligenomeSYMBOL)
data(ecoligenomeSTRAND)
ecoli.len
}
\format{
  The format for \code{ecoli.m52.genome} is \code{character}
  with genome sequence.
  The format for \code{ecoligenomeCHRLOC} is
  an environment (as a hash table). Each key
  is an Affyemtrix probe set ID, and each value is vector of
  two integers (begining and end - see the details
  below)
  The format for \code{ecoligenomeSYMBOL2AFFFY}
  is an environment (as a hash table). Each key is
  a gene symbol name.
  The format for \code{ecoligenomeSYMBOL}
  is an environment (as a hash table). Each key
  is an Affymetrix probe set id
  \code{ecoli.len} is a variable containing the size of the genome in
  \code{ecoli.m52.genome}.
}
\details{
  The environments \code{ecoligenomeSYMBOL2AFFFY} and
  \code{ecoligenomeSYMBOL} are like the ones in the
  data packages built by \code{annBuilder}.
  
  The environment
  \code{ecoligenomeCHRLOC} differs: two integers are associated
  with each key, one corresponds to the begining of the segment
  the other to the end.

  The environment ecoligenomeSTRAND returns a \code{logical}.
  \code{TRUE} means that the orientation is `+', \code{FALSE} means
  that the orientation is '-' (and \code{NA} is used when irrelevant
  for the key).
  
}
\source{
  \url{http://www.genome.wisc.edu/sequencing/k12.htm} and
  \url{http://www.biostat.harvard.edu/complab/dchip/info_file.htm}
}

\examples{
data(ecoli.m52.genome)


}
\keyword{datasets}
