scarm.filter <- function(x, 
                         right.width=15, 
                         min.left.width=right.width, 
                         min.width=floor(right.width/3), 
                         max.width=200, 
                         sign.level=0.001,
                         bound.slope.diff=0,
                         bound.noise.sd=0,
                         rtr=TRUE
                         ){

################
# Preparations #
################
    right.width <- round(right.width)
    min.left.width <- round(min.left.width)
    min.width <- round(min.width)
    max.width <- round(max.width)
    N <- length(x) 
    tseries <- x

##############################
# Stopping and Warning rules #
##############################

    if(missing(x))
        stop("The input data vector is missing with no default.\n")
        
    if(!is.numeric(x))
        stop("The data vector must be numeric.\n")
        
    if(N < min.width)
        stop("Length of time series must be at least 'min.width' =", min.width, ".\n")

    if(min.width < 5){
        min.width <- 5
        warning("'min.width' must be at least 5; 'min.width' is set to 5.\n")
    }
    
    if(right.width < 5){
        right.width <- 5
        warning("'right.width' must be at least 5; 'right.width' is set to 5.\n")
    }    
        
    if(min.left.width < 5){
        min.left.width <- 5
        warning("'min.left.width' must be at least 5; 'min.left.width' is set to 5.\n")
    }  
    
    if(right.width > min.left.width){
        min.left.width <- right.width
        warning("'min.left.width' must not be smaller than 'right.width'; 'min.left.width' is set to 'right.width' = ", right.width,".\n")
    }
    
    if(sign.level <= 0 | sign.level > 0.5)
        stop("sign.level must be a value in (0,0.5)")
        
    if(max.width < min.width)
        stop("'max.width' must not be smaller than 'min.width'.\n")
        
    if(max.width < right.width + min.left.width | max.width%%1!=0)
        stop("'max.width' must not be smaller than 'right.width'+'min.left.width' = ", right.width+min.left.width,".\n")
    
    if(bound.slope.diff < 0){
        bound.slope.diff <- 0
        warning("'bound.slope.diff' must be at least 0; 'bound.slope.diff' is set to 0.\n")
    }
    
    if(bound.noise.sd < 0){
        bound.noise.sd <- 0
        warning("'bound.noise.sd' must be at least 0; 'bound.noise.sd' is set to 0.\n")
    }
                                               
#############################################
# Internal functions and required data sets #
#############################################

    data(dfs)
    data(const.Q)
    data(var.n)

    get.B.t <- function(y){
        n <- length(y)
        j <- 1:n
        B.t <- matrix(NA,n,n)
        for (i in j) {
            B.t[j[j<i],i] <- B.t[i,j[j<i]] <- (y[i] - y[j[j<i]]) / (i - j[j<i])
        }
        return(B.t)
    }

    get.beta.RM <- function(B.t){
        beta.i <- apply(B.t,1,median, na.rm=T)
        beta.RM <- median(beta.i, na.rm=T)
        return(beta.RM)
    }
    
    get.mu.RM <- function(y,beta.RM){
        n <- length(y)
        j <- 1:n
        mu.RM <- median(y - (beta.RM*(j-n)), na.rm=T)
        return(mu.RM)
    }
    
    get.B1 <- function(y, ell){
        n <- length(y)
        r <- n - ell
        I.left <- 1:ell
        I.right <- (ell+1):n
        B1 <- matrix(NA,ell,r)
        for (i in I.left) {
            B1[i,I.right-ell] <- (y[i] - y[I.right]) / (i - I.right)
        }
        return(B1)
    }   
    
    get.Q <- function(y){ 
        y <- na.omit(y)
        n <- length(y)
        h <- numeric(n-2)
        for(i in 1:(n-2)){
            h[i] <- abs(y[i+1]-(y[i]+y[i+2])/2)
        }
        n <- length(h)+2
        if(n > 150){
            c.q <- 1.21 * (n/(n+0.44))
        } else {
            c.q <- const.Q[n]
        }
        Q.est <- c.q * sort(h)[floor(0.5*(n-2))]
        return(Q.est)
    }    

    # degrees of freedom for t-distribution to obtain critical values for the SCARM test statistic
    get.critval <- function(ell, r, sign.level){
        if(ell <= 100){
            ell <- ell - ell%%5
            r <- r - r%%5
            degree.of.freedom <- dfs[ell/5,r/5]
            critval <- qt(p=1-(sign.level/2), df=degree.of.freedom)
        } else {
            critval <- qnorm(p=1-(sign.level/2))
            
        }
        return(critval)
    }
    
    # this function delivers the scarm test statistic
    scarm.test.statistic <- function(y, B, ell, r, noise.sd){
        n <- ell+r
        left.sample <- y[1:ell]
        right.sample <- y[(ell+1):n]
        
        B.t.left <- B[1:ell,1:ell]
        B.t.right <- B[(ell+1):n,(ell+1):n]

        beta.left <- get.beta.RM(B.t.left)
        beta.right <- get.beta.RM(B.t.right)
        
        d.t <- beta.left - beta.right
        if(ell <= 300){
            v.left <- var.n[ell]
        } else {
            v.left <- 4.77e-07 + (17.71*(1/ell^3))
        }
        if(r <= 300){
            v.right <- var.n[r]
        } else {
            v.right <- 4.77e-07  + (17.71*(1/r^3))
        }
        var.t <- max(bound.noise.sd, noise.sd)^2 * (v.left+v.right)
        test.statistic <- d.t/sqrt(var.t)
        return(list(test.statistic=test.statistic, noise.sd=noise.sd, 
                    beta.left=beta.left, beta.right=beta.right, 
                    B.t.left=B.t.left, B.t.right=B.t.right))    
    }

####################
# Internal objects #
####################
    r <- round(right.width)
    ell.min <- round(min.left.width)
    n.min <- round(min.width)
    n.max <- round(max.width)
    signal.est <- slope.est <- slope.diff <- adapted.width <- scarm.statistic <- critvals <- noise.sd <- rep(NA,N)  
    n <- n.min
    I.win <- 1:n
    y <- tseries[I.win]
    B <- get.B.t(y)

###################
# Main  Algorithm #
###################
    for(i in n.min:N){
        if(all(is.na(y[(n-n.min+1):n])) || length(which(!is.na(y[(n-(min(r,n))+1):n]))) < n.min){
            # there are not enough recent observations!
            mu.est <- NA
            beta.est <- NA
            noise.sd.est <- NA
            # Update step
            if(i != N) {
                b.new     <- (tseries[i+1] - tseries[i+c((-n.min + 2):0)]) / (i+1 - (i+c((-n.min + 2):0)))
                B     <- cbind(rbind(B[(n-n.min+2):n, (n-n.min+2):n], b.new), c(b.new,NA))
                n <- n.min
                adapted.width[i] <- n
                y       <- tseries[(i - n + 2):(i + 1)]
            }   
        } else {
            if(n >= ell.min+r){
            # is n large enough for test application?
                ell <- n - r
                if(length(which(!is.na(y[1:ell])))>=ell.min/2 && length(which(!is.na(y[(ell+1):n])))>=r/2){
                # are there enough observations for testing? 
                    noise.sd.est <- get.Q(y)
                    critval <- get.critval(ell=ell, r=r, sign.level=sign.level)
                    scarm.test <- scarm.test.statistic(y=y, B=B, ell=ell, r=r, noise.sd=noise.sd.est)
                    scarm.statistic[i] <- scarm.test$test.statistic
                    slope.diff[i] <- scarm.test$beta.left - scarm.test$beta.right
                    critvals[i] <- critval
                    if(abs(scarm.test$test.statistic) > critval && abs(slope.diff[i]) > bound.slope.diff){ #!!
                    # does test reject the null hypothesis?
                        B <- B[(n-n.min+1):n,(n-n.min+1):n]
                        noise.sd.est <- get.Q(y)
                        
                        beta.est <- get.beta.RM(B)
                        mu.est <- get.mu.RM(beta.RM=beta.est, y=y[(n-n.min+1):n])
                        n <- n.min
                        if(rtr==TRUE){
                            if(n==n.min){
                                rtr.sample <- tseries[(i-n+1):i]
                            } else {
                                rtr.sample <- tseries[max(1,(i-r+1)):i]
                            }
                            min.rtr.sample <- min(rtr.sample, na.rm=T)
                            max.rtr.sample <- max(rtr.sample, na.rm=T)
                            mu.est <- max(min.rtr.sample, min(mu.est, max.rtr.sample))
                        }  
                        if(i != N) {
                            b.new <- (tseries[i+1] - tseries[i+c((-n + 1):0)]) / (i+1 - (i+c((-n + 1):0)))
                            B <- cbind(rbind(B, b.new), c(b.new,NA))
                            adapted.width[i] <- n
                            y <- tseries[(i - n + 1):(i + 1)]
                            n <- n+1
                        }
                    } else {
                    # difference of RM slopes is NOT too large => do not decrease time window
                        beta.est <- get.beta.RM(B)
                        mu.est <- get.mu.RM(beta.RM=beta.est, y=y)                    
                        adapted.width[i] <- n
                        if(i != N) {
                            if(n == n.max){ 
                            # is window width maximal?
                                b.new     <- (tseries[i+1] - tseries[i+c((-n + 2):0)]) / (i+1 - (i+c((-n + 2):0)))
                                B     <- cbind(rbind(B[-1, -1], b.new), c(b.new,NA))
                                y       <- tseries[(i - n + 2):(i + 1)]
                            } else {
                                b.new     <- (tseries[i+1] - tseries[i+c((-n + 1):0)]) / (i+1 - (i+c((-n + 1):0)))
                                B       <- cbind(rbind(B, b.new), c(b.new, NA))
                                y       <- tseries[(i - n + 1):(i + 1)]
                                n       <- n + 1
                            }
                        }  
                        if(rtr==TRUE){
                            if(n==n.min){
                                rtr.sample <- tseries[(i-n+1):i]
                            } else {
                                rtr.sample <- tseries[max(1,(i-r+1)):i]
                            }
                            min.rtr.sample <- min(rtr.sample, na.rm=T)
                            max.rtr.sample <- max(rtr.sample, na.rm=T)
                            mu.est <- max(min.rtr.sample, min(mu.est, max.rtr.sample))
                        }    
                    }
                } else {
                # there are not enough observations for testing!
                    B <- B[(n-r+1):n,(n-r+1):n]
                    noise.sd.est <- get.Q(y)
                    beta.est <- get.beta.RM(B)
                    mu.est <- get.mu.RM(beta.RM=beta.est, y=y[(n-r+1):n])
                    n <- r
                    if(rtr==TRUE){
                        rtr.sample <- tseries[max(1,(i-r+1)):i]
                        min.rtr.sample <- min(rtr.sample, na.rm=T)
                        max.rtr.sample <- max(rtr.sample, na.rm=T)
                        mu.est <- max(min.rtr.sample, min(mu.est, max.rtr.sample))
                    }
                    if(i != N) {
                        b.new <- (tseries[i+1] - tseries[i+c((-n + 1):0)]) / (i+1 - (i+c((-n + 1):0)))
                        B <- cbind(rbind(B, b.new), c(b.new,NA))
                        adapted.width[i] <- n
                        y <- tseries[(i - n + 1):(i + 1)]
                        n <- n+1
                    }  
               }   
           } else {
           # n is not large enough for test application!
                noise.sd.est <- get.Q(y)
                beta.est <- get.beta.RM(B)
                mu.est <- get.mu.RM(y=y, beta.RM=beta.est)
                adapted.width[i] <- n
                if(rtr==TRUE){
                    if(n==n.min){
                        rtr.sample <- tseries[(i-n+1):i]
                    } else {
                        rtr.sample <- tseries[max(1,(i-r+1)):i]
                    }
                    min.rtr.sample <- min(rtr.sample, na.rm=T)
                    max.rtr.sample <- max(rtr.sample, na.rm=T)
                    mu.est <- max(min.rtr.sample, min(mu.est, max.rtr.sample))
                }
                if(i != N) {
                    if(n == n.max){ 
                    # is window width maximal?
                        b.new     <- (tseries[i+1] - tseries[i+c((-n + 2):0)]) / (i+1 - (i+c((-n + 2):0)))
                        B     <- cbind(rbind(B[-1, -1], b.new), c(b.new,NA))
                        y       <- tseries[(i - n + 2):(i + 1)]
                    } else {
                        b.new     <- (tseries[i+1] - tseries[i+c((-n + 1):0)]) / (i+1 - (i+c((-n + 1):0)))
                        B       <- cbind(rbind(B, b.new), c(b.new, NA))
                        y       <- tseries[(i - n + 1):(i + 1)]
                        n       <- n + 1
                    }
                }       
            }
        }
        signal.est[i] <- mu.est
        slope.est[i] <- beta.est
        noise.sd[i] <- noise.sd.est
    }
    result <- list(signal.est=signal.est, slope.est=slope.est, adapted.width=adapted.width, 
                   test.statistic=scarm.statistic, critvals=critvals, noise.sd=noise.sd, slope.diff=slope.diff, 
                   tseries=tseries, r=r, ell.min=ell.min, n.min=n.min, n.max=n.max, sign.level=sign.level, bound.slope.diff=bound.slope.diff, bound.noise.sd=bound.noise.sd)
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
  par(mfrow=c(2,1))
  plot(x$tseries, type="l", main="Data and SCARM signal xaction", xlab="time", ylab="", ...)
  lines(x$signal, col=2, lwd=1)
  legend(x="topleft", legend=c("Data", "SCARM signal xaction"),lty=c(1, 1),col=c("black", "red"),lwd=c(1, 1), bty="n")
  plot(x$adapted.width, type="l",  main="Adapted window widths", xlab="time", ylab="", ...)
}
