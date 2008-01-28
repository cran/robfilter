WRMfit <- function(xdat, ydat, x0=xdat[length(xdat)], weights=rep(1,length(xdat))){
  weights     <- weights/sum(weights)
  xdat        <- xdat[weights!=0]
  ydat        <- ydat[weights!=0]
  weights     <- weights[weights!=0]
  sm          <- rep(0,length(xdat))
  for (i in 1:length(xdat)){
     sm[i] <- (weighted.quantile.top((ydat[i]-ydat[-i])/(xdat[i]-xdat[-i]), weights[-i], 0.5) + weighted.quantile((ydat[i]-ydat[-i])/(xdat[i]-xdat[-i]),weights[-i],0.5))/2
  }
  bm <- (weighted.quantile.top(sm, weights, 0.5) + weighted.quantile(sm, weights, 0.5))/2
  am <- (weighted.quantile.top(ydat-bm*(xdat-x0), weights, 0.5) + weighted.quantile(ydat-bm*(xdat-x0), weights, 0.5))/2
  c(am,bm)  # am is fitted value
}
