\name{ecoligenome.operon}
\alias{ecoligenome.operon}
\docType{data}
\title{ Known operon in E.coli - data.frame}
\description{
  The known operon in the Escherichia coli genome.
}
\usage{data(ecoligenome.operon)}
\format{
  A data frame with 932 observations (genes ) on the following 4 variables.
  \describe{
    \item{gene.name}{a character vector}
    \item{gene.annotation}{a character vector}
    \item{operon.name}{a factor with levels the names of the operons}
    \item{operon.comments}{a factor with levels the comments for the operons}
  }
}
\details{
  For some operons, the source of information specifies the existence of
  regulating elements such as promoter, terminator, box, etc\ldots
  In those cases, the \code{gene.name} is set to \code{"Regulation"},
  and the \code{gene.annotation} gives what kind of regulating element
  it is. If volonteers, it would be neat to map those on the genome\ldots
  Besides that, not much to add. The data structure is fairly straightforward.
}
\source{
  Built from the webpage:
  \url{http://www.cib.nig.ac.jp/dda/backup/taitoh/ecoli.operon.html}
}
% \references{
%   ~~ possibly secondary sources and usages ~~
% }
\examples{
## something that might be useful when working with Affymetrix data
library(Biobase)
data(ecoligenome.operon)
data(ecoligenomeSYMBOL2AFFY)
ecoligenome.operon$affyid <- unname(unlist(multiget(ecoligenome.operon$gene.name, ecoligenomeSYMBOL2AFFY)))
}
\keyword{datasets}
