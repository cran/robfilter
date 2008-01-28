weighted.quantile <- function(dat, weights, p){
  weights <- weights[order(dat)]
  dat     <- sort(dat)
  q       <- sum(weights)*p
  m <- w  <-0
  for (i in 1 : length(weights)){
    w <- w+weights[i]
    if (w >= q) {m<- dat[i] ; break}
  }
  m
}