## RCI script
##  Here is the RCI function I've used 
##  x and y are vectors of counts in a datatable (e.g. counts of poor and non-poor in a zone)
##  sort.var is the vector variable by which zones are sorted (e.g distance)
##  The calculation uses matrices for speed (at one point I needed faster outputs)
##  t() gives a transpose of a matrix or vector
##  %*% is matrix product

#ME: BE AWARE ORDER IS REVERSED I.E. IF Y IS CENTRALISED, RCI WILL BE POSITIVE
rci <- function(x, y, sort.var) {
  ordering <- order(sort.var)
  N <- length(ordering)
  x <- x[ordering]
  y <- y[ordering]
  
  ##  cumulative sums
  csum.x <- cumsum(x) / sum(x)
  csum.y <- cumsum(y) / sum(y)
  
  ## The below has been checked to be right;
  out <- t(csum.y[-N]) %*% csum.x[-1] - t(csum.y[-1]) %*% csum.x[-N]
  return(out)
}
