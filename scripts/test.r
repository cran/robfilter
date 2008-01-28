## Required dependencies:
require(robustbase)
require(MASS)

## Clean start:
rm(list=ls(all=TRUE))

Require <- function(f) 
  source(paste("../R", f, sep="/"))

set.seed(43)
n <- 100
y <- cumsum(runif(n)) - .5*(1:n)
o <- sample(1:n, n/10)
j <- n/2
## Add noise
y[o] <- y[o] + rnorm(n/10)
## Add jump
y[j:n] <- y[j:n] + 2

Require("trend.R")
Require("scale.R")
Require("robust-filter.R")
Require("double-window.R")
Require("hybrid.R")
load("../R/sysdata.rda")

message(":: Testing 'robust.filter' ...") 
y.rf <- robust.filter(y, 
trend.extraction="LMS") print(h.rf) plot(y.rf) 

message(":: Testing 'double.window' ...")
y.dw <- double.window(y, outer.n=11, inner.n=5)
print(y.dw)
plot(y.dw)

message(":: Testing 'hybrid' ...")
y.hy <- hybrid(y, n=11)
print(y.hy)
plot(y.hy, ylim=c(-6, 2))
