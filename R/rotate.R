xy2polar <- function(x, y) {
  if (missing(y)) {
    y <- x$y
    x <- x$x
  }
  rho <- sqrt(x^2 + y^2)
  theta <- atan(y/x)
  return(list(rho=rho, theta=theta))
}

polar2xy <- function(rho, theta) {
  if (missing(theta)) {
    theta <- rho$theta
    rho <- rho$rho
  }
  x <- rho * cos(theta)
  y <- rho * sin(theta)
  return(list(x=x, y=y))
}

rotate <- function(x, y, alpha) {
  pol <- xy2polar(x, y)
  pol$theta <- pol$theta + alpha
  xy <- polar2xy(pol)
  return(xy)
}
                    
