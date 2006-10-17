.buildBNUM2GENBANK <- function() {
  return(.buildGENEPRODFUN(what="genbank"))
}

.buildBNUM2GENEPRODUCT <- function() {
  return(.buildGENEPRODFUN(what="gene product"))
}

.buildGENEPRODFUN <- function(my.file=NULL, my.url="http://genprotec.mbl.edu/files/geneproductfunctions.txt", what="") {
  choices <- c("bnum", "genbank", "gene", "gene type", "gene product")
  what <- match.arg(what, choices)
  what.i <- match(what, choices)
  env.x2y <- new.env(hash=TRUE)
  ##env.genbank <- new.env(hash=TRUE)
  if (is.null(my.file))
    con <- url(my.url, open="r")
  else
    con <- file(my.file, open="r")

  ## skip comments
  mycomments <- ""
  line <- readLines(con, n=1)
  while (length(grep("^bnum", line)) == 0) {
    line <- readLines(con, n=1)
    mycomments <- paste(mycomments, line, sep="\n")
  }
  mycomments <- paste(mycomments, "\n", sep="")
  line <- readLines(con, n=1)

  ##
  while (length(line) != 0) {
    if (line != "") {
      m <- strsplit(line, "\t", extended=TRUE)[[1]]
      bnum <- m[1]
      genbank <- m[what.i]
      if (genbank == "0")
        genbank <- NA
      if (exists(bnum, envir=env.x2y))
        assign(bnum, unique(c(genbank, get(bnum, envir=env.x2y))), envir=env.x2y)
      else
        assign(bnum, genbank, envir=env.x2y)
    }
    line <- readLines(con, n=1)
  }

  attr(env.x2y, "comments") <- mycomments
  return(env.x2y)

}


.buildMultiFun <- function(my.file=NULL, my.url="http://genprotec.mbl.edu/files/MultiFun.txt") {
  env.multiFun <- new.env(hash=TRUE)

  if (is.null(my.file))
    con <- url(my.url, open="r")
  else
    con <- file(my.file, open="r")

  ## skip comments
  mycomments <- ""
  line <- readLines(con, n=1)
  while (length(grep("^[0-9]", line)) == 0) {
    line <- readLines(con, n=1)
    mycomments <- paste(mycomments, line, sep="\n")
  }
  mycomments <- paste(mycomments, "\n", sep="")

  ##
  while (length(line) != 0) {
    if (line != "") {
      ## trim leading white spaces
      m <- regexpr("(\\\w\\\.?)+", line, perl=TRUE)
      multiFun <- substr(line, m, m+attr(m, "match.length")-1)
      multiFun <- sub("\\\.$", "", multiFun)
      multiFunAnnotation <- substr(line, m + attr(m, "match.length") + 1, nchar(line))
      assign(multiFun, multiFunAnnotation, envir=env.multiFun)
    }
    line <- readLines(con, n=1)
  }
  attr(env.multiFun, "comments") <- mycomments
  return(env.multiFun)
}

.buildBNUM2MULTIFUN <- function(my.file=NULL, my.url="http://genprotec.mbl.edu/files/multifunassignments.txt") {
  env.bnum2multifun <- new.env(hash=TRUE)

  if (is.null(my.file))
    con <- url(my.url, open="r")
  else
    con <- file(my.file, open="r")

  ## skip comments
  mycomments <- ""
  line <- readLines(con, n=1)
  while (length(grep("^Bnum", line)) == 0) {
    line <- readLines(con, n=1)
    mycomments <- paste(mycomments, line, sep="\n")
  }
  mycomments <- paste(mycomments, "\n", sep="")
  line <- readLines(con, n=1)

  ##
  while (length(line) != 0) {
    if (line != "") {
      ## trim leading white spaces
      m <- regexpr("b\\\d+", line, perl=TRUE)
      bnum <- substr(line, m, m+attr(m, "match.length")-1)
      m <- regexpr("\\\s\\\S+", line, perl=TRUE)
      multiFun <- substr(line, m+1, m+attr(m, "match.length")-1)
      if (exists(bnum, envir=env.bnum2multifun))
        assign(bnum, unique(c(multiFun, get(bnum, envir=env.bnum2multifun))), envir=env.bnum2multifun)
      else
        assign(bnum, multiFun, envir=env.bnum2multifun)
    }
    line <- readLines(con, n=1)
  }

  return(env.bnum2multifun)
}

.buildBNUM2SYMBOL <- function(my.url="http://genprotec.mbl.edu/files/MultiFun.txt") {

}

.buildMULTIFUN2GO <- function(filename) {
  env.multiFun2GO <- new.env(hash=TRUE)
  con <- file(filename, open="r")

  ## skip comments
  mycomments <- ""
  line <- readLines(con, n=1)
  while (length(grep("^!", line)) > 0) {
      line <- readLines(con, n=1)
      mycomments <- paste(mycomments, line, sep="\n")
    }
  mycomments <- paste(mycomments, "\n", sep="")
  ##

  while (length(line) != 0) {
    ## get MultiFun number
    m <- regexpr("^MultiFun:(\\.+?) ", line, perl=TRUE)
    multiFun <- substr(line, 10, attr(m, "match.length") - 1)
    ## get the GOs
    m <- regexpr("> \\.+$", line, perl=TRUE)
    GOs <- substr(line, m+2, m+attr(m, "match.length"))
    if (GOs == "GO:.") {
      ## no go ;)
      GOs <- NA
    } else {
      GOs <- strsplit(GOs, " > ")[[1]]
      GOs <- strsplit(GOs, " ; ")
      GOs <- unlist(lapply(GOs, function(x) x[2]))
    }
    assign(multiFun, GOs, envir=env.multiFun2GO)

    line <- readLines(con, n=1)
  }
  attr(env.multiFun2GO, "comments") <- mycomments
  return(env.multiFun2GO)
}


.buildMultiFunGraphNEL <- function(filename) {
  require("graph") || stop("The graph packages is needed for this operation")

  nid.i <- 1
  nname.i <- 2

  linesInFile <- readLines(filename)

  tmp <- strsplit(linesInFile, ";")

  r <- lapply(tmp, function(x, y) strsplit(x[nid.i], y), "\\\.")
  r.names <- unlist(lapply(tmp, function(x, y) x[nname.i]))

  mFunNodenames <- unique(unlist(lapply(tmp, function(x) x[nname.i])))
  nodename2i <- new.env(hash=TRUE)
  nodename2nodeid <- new.env(hash=TRUE)
  nodeid2nodename <- new.env(hash=TRUE)

  multiassign(mFunNodenames, seq(along=mFunNodenames), nodename2i)
  mFunEdges <- vector("list", length=length(mFunNodenames))
  names(mFunEdges) <- mFunNodenames
  multiassign(unlist(lapply(tmp, function(x) x[nid.i])), unlist(lapply(tmp, function(x) x[nname.i])),
              envir=nodeid2nodename)

  for (i in seq(along=r)) {
    n <- length(r[[i]][[1]])
    if (n == 1)
      next
    parent <- paste(r[[i]][[1]][seq(1, n-1, length=n-1)], collapse=".")
    if (! exists(parent, nodeid2nodename))
      next
    parent.i <- get(get(parent, nodeid2nodename), nodename2i)
    child.i <- get(tmp[[i]][nname.i], nodename2i)
    ##parent.i <- get(get(parent, nodeid2nodename), nodename2i)
    ##child.i <- i
    if (is.null(mFunEdges[[parent.i]])) {
      mFunEdges[[parent.i]]$edges = child.i
      mFunEdges[[parent.i]]$weights = 1
    } else {
      mFunEdges[[parent.i]]$edges = c(mFunEdges[[parent.i]]$edges, child.i)
      mFunEdges[[parent.i]]$weights = c(mFunEdges[[parent.i]]$weights, 1)
    }

    if (is.null(mFunEdges[[child.i]])) {
      ##meshedges[[child.i]]$edges = parent.i
###meshedges[[parent.i]]$edges = paste(r[[i]][[1]], collapse=".")
                                        #meshedges[[child.i]]$weights = 1
    } else {
###meshedges[[parent.i]]$edges = c(meshedges[[parent.i]]$edges, paste(r[[i]][[1]], collapse="."))
      ##meshedges[[child.i]]$edges = c(meshedges[[child.i]]$edges, parent.i)
      ##meshedges[[child.i]]$weights = c(meshedges[[child.i]]$weights, 1)
    }
  }

  ##meshnodenames <- seq(along=meshnodenames)
  ##names(meshedges) <- meshnodenames

  gmesh <- new("graphNEL", nodes=mFunNodenames, edgeL=mFunEdges, edgemode="directed")

  return(gmesh)
}


linkedmultiget <- function(x, envir.list=list(), unique=TRUE) {

  if (! is.character(x))
    stop("x must be a vector of mode 'character'")

  f <- function(x, y) {
    tmp <- mget(x, envir=y, ifnotfound=NA)
    tmp <- unlist(tmp)
    if (all(is.na(tmp)))
      tmp <- as.character(tmp)
    if (! is.character(tmp))
      stop("Values in environments must be of mode 'character'")
    return(tmp)
  }

  ##r <- vector("list", length=length(x))
  r <- as.list(x)
  for (i in seq(along=envir.list)) {
    r <- lapply(r, f, envir.list[[i]])
    if (unique)
      r <- lapply(r, unique)
  }
  names(r) <- x
  return(r)
}
