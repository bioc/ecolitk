\name{linkedmultiget}
\alias{linkedmultiget}
\title{ A function to look for values across linked environments}
\description{
  A function to look for values across linked environments}
}
\usage{
linkedmultiget(x, envir.list = list(), unique = TRUE)
}

\arguments{
  \item{x}{The keys in the first environment in the list.}
  \item{envir.list}{A list of environments.}
  \item{unique}{Simplify the list returned by ensuring that the values
    for each key are unique.}
}
\details{
  Environments can be considered as hashtables. The keys are obviously
  strings, but in some cases the associated values are also
  strings. This is the case for annotation environments (as built with
  the package \code{AnnBuilder}). This function helps to look for values
  across several environments: the keys have associated values in a
  first environment, these values are used as keys in the second
  environments, etc...
}
\value{
  A list of length the length of \code{x}.
}
\author{ Laurent Gautier }
\seealso{ \code{\link{mget}} }
\examples{
data(ecoligenomeBNUM)
data(ecoligenomeBNUM2MULTIFUN)
data(multiFun)

## get 5 Affymetrix IDs
set.seed(456)
my.affyids <- sample(ls(ecoligenomeBNUM), 5)

## get the MULTIFUN annotations for them
r <- linkedmultiget(my.affyids, list(ecoligenomeBNUM,
                    ecoligenomeBNUM2MULTIFUN, multiFun))

print(r)
}
\keyword{ manip }
