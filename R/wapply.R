wstringapply <- function(x, SIZE, SLIDE, FUN, ...) {
  n <- (nchar(x) - SIZE) %/% SLIDE
  res <- vector("list", length = n)
  for (i in seq(0, n-1, length=n)) {
    ##x.start <- (i-1) * SIZE + 1
    ##x.stop <- min(nchar(x), i * SIZE)
    w.x <- substring(x, i * SLIDE + 1, i * SLIDE + SIZE)
    res[[i+1]] <- FUN(w.x, ...)
    ##res[[i]] <- FUN(w.x)
  }
  return(res)
}

# wapply <- function(x, SIZE, FUN, ...) {
#   n <- length(x)
#   res <- vector("list", length = n-1)
#   for (i in seq(1, n-1, length=n-1)) {
#     w.x <- x[i, i+SIZE]
#     res[[i]] <- FUN(w.x, ...)
#     ##res[[i]] <- FUN(w.x)
#   }
#   return(res)
# }

gccontent <- function(x) {x <- toupper(x); n <- nchar(x); sum(as.integer(strsplit(x, "")[[1]] %in% c("G", "C")))/n}
