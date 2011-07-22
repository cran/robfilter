scarm.filter <- function(x, 
                         right.width=30, 
                         min.left.width=right.width, 
                         min.width=right.width/3, 
                         max.width=180, 
                         sign.level=0.001,
                         bound.noise.sd=1){
          
    # stopping criteria
    if(missing(x))
        stop("The input data vector is missing with no default.\n")
    if(!is.numeric(x))
        stop("The data vector must be numeric.\n")
    N <- length(x) 
    if(N < min.width)
        stop("The data vector must be longer than min.width.\n")
    if(min.width < 1 | min.width%%1!=0 )
        stop("min.width must be a positive integer >=1")
    if(right.width < 5 | right.width%%1!=0)
        stop("right.width must be a positive integer >=5")
    if(min.left.width < 5 | min.left.width%%1!=0)
        stop("min.left.width must be a positive integer >=5")
    if(right.width > min.left.width)
        stop("min.left.width must not be smaller than right.width")
    if(sign.level <= 0 | sign.level > 0.5)
        stop("sign.level must be a value in (0,0.5)")
    if(max.width <= min.width | max.width%%1!=0)
        stop("max.width must be a positive integer > min.width")
    if(max.width <= right.width + min.left.width | max.width%%1!=0)
        stop("max.width must be a positive integer > right.width+min.left.width")
    if(bound.noise.sd <= 0)
        stop("bound.noise.sd must be > 0")
                                           
    # load required data sets
    data(dfs)
    data(const.Q)
    data(var.n)

    # required functions
    RM <- function(y){
        n <- length(y)
        s <- 1:n
        medquotient.t <- c()
        for (k in s) {
            all.slopes.with.k <- y[k] - y[s[s!=k]]
            medquotient.t <- c(medquotient.t, median(all.slopes.with.k / (k - s[s!=k]) ,na.rm=TRUE))
        }
        beta.RM <- median(medquotient.t, na.rm=TRUE)
        mu.RM <- median(y - beta.RM*(s-n), na.rm=TRUE)
        y.hat <- mu.RM + beta.RM*(s-n)
        f <- which(is.na(y))
        if(length(f) > 0){
            y[f] <- y.hat[f]
        }
        res <- y - y.hat
        list(mu.RM=mu.RM, beta.RM=beta.RM, y.hat=y.hat, res=res)
    }
    Q <- function(x){ 
        n <- length(x)
        if(n > 150){
            c.q <- 1.21 * (n/(n+0.44))
        } else {
            c.q <- const.Q[n]
        }
        h <- numeric(n-2)
        for(i in 1:(n-2)){
            h[i] <- abs(x[i+1]-(x[i]+x[i+2])/2)
        }
        h <- sort(h)
        return(c.q * h[floor(0.5*(n-2))])
    }
    
    # set parameters
    right.width <- round(right.width)
    min.left.width <- round(min.left.width)
    min.width <- round(min.width)
    max.width <- round(max.width)

    # empty result vectors
    signal.est <- adapted.width <- scarm.statistic <- critvals <- noise.sd <- rep(NA,N)  
    
    # degrees of freedom for t-distribution to obtain critical values for the SCARM test statistic
    get.critval <- function(left.width, right.width, sign.level){
        if(left.width <= 100){
            ell <- left.width - left.width%%5
            r <- right.width - right.width%%5
            degree.of.freedom <- dfs[ell/5,r/5]
            critval <- qt(p=1-(sign.level/2), df=degree.of.freedom)
        } else {
            critval <- qnorm(p=1-(sign.level/2))
            
        }
        return(critval)
    }
    
    # this function delivers the scarm test statistic
    scarm.test.statistic <- function(whole.sample, left.width, right.width){
        left.sample <- whole.sample[1:left.width]
        right.sample <- whole.sample[(left.width+1):length(whole.sample)]
        left.RM.fit <- RM(left.sample)
        right.RM.fit <- RM(right.sample)
        d.t <- left.RM.fit$beta.RM - right.RM.fit$beta.RM
        if(length(left.sample) <= 300){
            v.left <- var.n[length(left.sample)]
        } else {
            v.left <- 4.77e-07 + (17.71*(1/length(left.sample)^3))
        }
        if(length(right.sample) <= 300){
            v.right <- var.n[length(right.sample)]
        } else {
            v.right <- 4.77e-07  + (17.71*(1/length(right.sample)^3))
        }
        noise.sd <- Q(whole.sample)
        var.t <- max(bound.noise.sd, noise.sd)^2 * (v.left+v.right)
        test.statistic <- d.t/sqrt(var.t)
        return(list(test.statistic=test.statistic,noise.sd=noise.sd))    
    }
    
    # the main algorithm
    whole.width <- min.width
    for(i in min.width:N){
        # is whole.width large enough for test application?
        if(whole.width >= min.left.width+right.width){
            left.width <- whole.width - right.width
            whole.sample <- x[(i-whole.width+1):i]
            if(length(which(!is.na(whole.sample[1:left.width])))>=min.left.width/2 && length(which(!is.na(whole.sample[(left.width+1):whole.width])))>=right.width/2){
                critval <- get.critval(left.width, right.width, sign.level)
                scarm.test <- scarm.test.statistic(whole.sample, left.width, right.width)
                noise.sd[i] <- scarm.test$noise.sd
                # apply test
                if(abs(scarm.test$test.statistic) > critval){
                    whole.width <- min.width
                    whole.sample <- x[(i-whole.width+1):i]
                }
                mu.est <- RM(whole.sample)$mu.RM     
                adapted.width[i] <- whole.width    
                whole.width <- min(whole.width+1,max.width)  
                scarm.statistic[i] <- scarm.test$test.statistic
                critvals[i] <- critval
            } else {
                mu.est <- NA
                whole.width <- min.width
                adapted.width[i] <- whole.width
            }
        } else {
            whole.sample <- x[(i-whole.width+1):i]
            if(length(which(!is.na(whole.sample)))>=min.width){
                mu.est <- RM(whole.sample)$mu.RM
                adapted.width[i] <- whole.width
                whole.width <- min(whole.width+1,max.width)
            } else {
                whole.width <- min.width
                adapted.width[i] <- whole.width
                mu.est <- NA
            }
        }
        # store RM signal extraction
        signal.est[i] <- mu.est      
    }
    result <- list(signal=signal.est, adapted.width=adapted.width, test.statistic=scarm.statistic, critvals=critvals, noise.sd=noise.sd,
                x=x, right.width=right.width, min.left.width=min.left.width, min.width=min.width, max.width=max.width, sign.level=sign.level, bound.noise.sd=bound.noise.sd)
    return(structure(result, class = "scarm.filter"))
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Default output
print.scarm.filter <- function(x, ...) {
  N <- length(x$signal)
  cat('$signal \n')  
    if(N <= 100){
    print(x$signal,...)
  } else {
    print(x$signal[1:50])
    cat('Only the first 50 signal estimations are printed.\n')    
  }
}

## # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
## Default plot
plot.scarm.filter <- function(x, ...){
  plot(x$x, type="l", main="Data and SCARM signal extraction", xlab="time", ylab="", ...)
  lines(x$signal, col=2, lwd=2)
  legend(x="topleft", legend=c("Data", "SCARM signal extraction"),lty=c(1, 1),col=c("black", "red"),lwd=c(1, 2))
}
