.buildMultiFun <- function(my.url="http://genprotec.mbl.edu/files/MultiFun.txt") {
  env.multiFun <- new.env(hash=TRUE)
  con <- url(my.url, open="r")
  
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
  return(env.multiFun)
}

.buildBNUM2MULTIFUN <- function(my.url="http://genprotec.mbl.edu/files/multifunassignments.txt") {
  env.bnum2multifun <- new.env(hash=TRUE)
  con <- url(my.url, open="r")
  
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
    m <- regexpr("^MultiFun:(\.+?) ", line, perl=TRUE)
    multiFun <- substr(line, 10, attr(m, "match.length") - 1)
    ## get the GOs
    m <- regexpr("> \.+$", line, perl=TRUE)
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
