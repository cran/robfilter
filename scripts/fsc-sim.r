##
## fsc-sim.r - Finite Sample Correction coef simulation
##

## Repeatable research:
set.seed(42)

## Sample size range:
r <- 1:100

## Simulate coef. for function 'fun'
sim.corr <- function(fun) {
  res <- numeric(length(r))
  cat(":: Simulating .. ")
  for (x in r) {
    d <- replicate(10000, 1/fun(rnorm(x)))
    res[x] <- mean(d)
    cat(round(res[x], 4), ".. ")    
  }
  return(res)
}

## lsh estimator:
lsh <- function(res) {
  n <- length(res)
  m <- ceiling((n - 1)/2)
  x <- sort(res)
  return(min(abs(x[(m+1):n] - x[1:(n-m)])))
}

## Read timecor data:
timecor.dat <- read.table("../data/timecor.dat")
timecorrection <- timecor.dat[,-(17:18)]
colnames(timecorrection)<- c("T-MAD","T-LSH","T-QN","T-SN",
                             "L-MAD","L-LSH","L-QN","L-SN",
                             "M-MAD","M-LSH","M-QN","M-SN",
                             "W-MAD","W-LSH","W-QN","W-SN")


## Simulate sizecor data:
sizecorrection <- data.frame(MAD=sim.corr(function(x) mad(x, constant=1)),
                             QN=sim.corr(function(x) Qn(x, constant=1)),
                             SN=sim.corr(function(x) Sn(x, constant=1)),
                             LSH=sim.corr(function(x) lsh(x, c=1))
                             )

save(timecorrection, sizecorrection, file="sysdata.rda")
