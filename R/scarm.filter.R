data(dfs)
data(const.Q)
data(var.n)

scarm.filter <- function(x, 
                         right.width=15, 
                         min.left.width=right.width, 
                         min.width=floor(right.width/3), 
                         max.width=200, 
                         sign.level=0.001,
                         bound.noise.sd=0.01,
                         rtr=TRUE,
                         noise.sd.est.method="Q.adj", # "MAD","Qn","SD"
                         autocorrelations="no"        # high.positive, moderate.positive, small.positive, small.negative, moderate.negative, high.negative
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
        warning("'min.width' must be an integer >=5; 'min.width' is set to 5.\n")
    }
    
    if(right.width < 5){
        right.width <- 5
        warning("'right.width' must be an integer >=5; 'right.width' is set to 5.\n")
    }    
        
    if(min.left.width < 5){
        min.left.width <- 5
        warning("'min.left.width' must be an integer >=5; 'min.left.width' is set to 5.\n")
    }  
    
    if(right.width > min.left.width){
        min.left.width <- right.width
        warning("'min.left.width' must not be smaller than 'right.width'; 'min.left.width' is set to 'right.width' = ", right.width,".\n")
    }
    
    if(sign.level <= 0 | sign.level > 0.5)
        stop("'sign.level' must be a value in (0,0.5)")
        
    if(max.width < min.width)
        stop("'max.width' must not be smaller than 'min.width'.\n")
        
    if(max.width < right.width + min.left.width | max.width%%1!=0)
        stop("'max.width' must not be smaller than 'right.width'+'min.left.width' = ", right.width+min.left.width,".\n")
    
    if(bound.noise.sd <= 0){
        bound.noise.sd <- 0.01
        warning("'bound.noise.sd' must be a value >0; 'bound.noise.sd' is set to 0.01.\n")
    }
                                               
#############################################
# Internal functions and required data sets #
#############################################

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
        if(Q.est==0){
            if(any(h!=0)){
                Q.est <- h[h!=0][1]
            } else {
                Q.est <- bound.noise.sd
            }
        }
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
    
    if(autocorrelations=="no"){
        emp.var <- var.n[,1]
        v.model <- function(n){
            out <- 4.77e-07 + (17.71*(1/n^3))
            return(out)
        }
    }
    if(autocorrelations=="high.positive"){
        emp.var <- var.n[,2]
        v.model <- function(n){
            out <- -0.0003318919 + (0.0476670483/n) + (2.2422001348/(n^2)) + (-5.8731008014/(n^3))
            return(out)
        }
    }  
    if(autocorrelations=="moderate.positive"){
        emp.var <- var.n[,3]
        v.model <- function(n){
            out <- 0.0007355669 + (-0.1485876684/n) + (5.5086877566/(n^2)) + (-5.3249305456/(n^3))
            return(out)
        }
    }   
    if(autocorrelations=="small.positive"){
        emp.var <- var.n[,4]
        v.model <- function(n){
            out <- 0.000473727 + (-0.086723762/n) + (2.442400804/(n^2)) + (10.578915504/(n^3))
            return(out)
        } 
    }   
    if(autocorrelations=="small.negative"){
        emp.var <- var.n[,5]
        v.model <- function(n){
            out <- 0.0002385982 + (-0.0395730416/n) + (0.7798496595/(n^2)) + (9.7595002645/(n^3))
            return(out)
        } 
    } 
    if(autocorrelations=="moderate.negative"){
        emp.var <- var.n[,6]
        v.model <- function(n){
            out <- 0.0002560433 + (-0.0434806450/n) + (0.9224130380/(n^2)) + (5.4903584072/(n^3))
            return(out)
        }
    }  
    if(autocorrelations=="high.negative"){
        emp.var <- var.n[,7]
        v.model <- function(n){
            out <- 0.0003514981 + (-0.0629949278/n) + (1.6446684583/(n^2)) + (-3.1859306183/(n^3))
            return(out)
        }
    }

    # this function delivers the scarm test statistic
    scarm.test.statistic <- function(y, B, ell, r){
        n <- ell+r
        left.sample <- y[1:ell]
        right.sample <- y[(ell+1):n]
        B.t.left <- B[1:ell,1:ell]
        B.t.right <- B[(ell+1):n,(ell+1):n]
        beta.left <- get.beta.RM(B.t.left)
        beta.right <- get.beta.RM(B.t.right)
        if(noise.sd.est.method!="Q.adj"){
            mu.left <- get.mu.RM(beta.RM=beta.left, y=left.sample)
            mu.right <- get.mu.RM(beta.RM=beta.right, y=right.sample)
            left.line <- mu.left + beta.left*((1-ell):0)
            right.line <-  mu.right + beta.right*((1-r):0)
            res.left <- left.sample - left.line
            res.right <- right.sample - right.line
            res <- c(res.left,res.right)
        }
        if(noise.sd.est.method=="MAD"){
            noise.sd.est <- mad(res, na.rm=T)
        }
        if(noise.sd.est.method=="Qn"){
            noise.sd.est <- Qn(na.omit(res))
        }
        if(noise.sd.est.method=="SD"){
            noise.sd.est <- sd(res, na.rm=T)
        }
        noise.sd.est <- max(noise.sd.est, bound.noise.sd)
        d.t <- beta.left - beta.right
        
        if(ell>300){
            v.left <- v.model(ell)
        } else {
            v.left <- emp.var[ell]
        }
        if(r>300){
            v.right <- v.model(r)
        } else {
            v.right <- emp.var[r]
        }
  
        var.t <- noise.sd.est^2 * (v.left+v.right)
        test.statistic <- d.t/sqrt(var.t)
        return(list(test.statistic=test.statistic, beta.left=beta.left, beta.right=beta.right, noise.sd.est=noise.sd.est))
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
                b.new <- (tseries[i+1] - tseries[i+c((-n.min + 2):0)]) / (i+1 - (i+c((-n.min + 2):0)))
                B <- cbind(rbind(B[(n-n.min+2):n, (n-n.min+2):n], b.new), c(b.new,NA))
                n <- n.min
                adapted.width[i] <- n
                y <- tseries[(i-n+2):(i+1)]
            }   
        } else {
            if(n >= ell.min+r){
            # is n large enough for test application?
                ell <- n - r
                if(length(which(!is.na(y[1:ell])))>=ell.min/2 && length(which(!is.na(y[(ell+1):n])))>=r/2){
                # are there enough observations for testing? 
                    if(noise.sd.est.method=="Q.adj"){
                        noise.sd.est <- get.Q(y)
                    }
                    critval <- get.critval(ell=ell, r=r, sign.level=sign.level)
                    scarm.test <- scarm.test.statistic(y=y, B=B, ell=ell, r=r)
                    scarm.statistic[i] <- scarm.test$test.statistic
                    slope.diff[i] <- scarm.test$beta.left - scarm.test$beta.right
                    critvals[i] <- critval
                    noise.sd[i] <- noise.sd.est <- scarm.test$noise.sd.est
                    if(abs(scarm.test$test.statistic) > critval){ #!!
                    # does test reject the null hypothesis?
                        B <- B[(n-n.min+1):n,(n-n.min+1):n]            
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
                    beta.est <- get.beta.RM(B)
                    mu.est <- get.mu.RM(beta.RM=beta.est, y=y[(n-r+1):n])
                    regression.line <- mu.est + beta.est*((1-r):0)
                    res <- y[(n-r+1):n] - regression.line
                    if(noise.sd.est.method=="Q.adj"){
                        noise.sd.est <- get.Q(y)
                    }
                    if(noise.sd.est.method=="MAD"){
                        noise.sd.est <- mad(res, na.rm=T)
                    }
                    if(noise.sd.est.method=="Qn"){
                        noise.sd.est <- Qn(na.omit(res))
                    }
                    if(noise.sd.est.method=="SD"){
                        noise.sd.est <- sd(res, na.rm=T)
                    }
                    noise.sd.est <- max(bound.noise.sd, noise.sd.est)
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
                beta.est <- get.beta.RM(B)
                mu.est <- get.mu.RM(y=y, beta.RM=beta.est)
                regression.line <- mu.est + beta.est*((1-n):0)
                res <- y - regression.line
                if(noise.sd.est.method=="Q.adj"){
                    noise.sd.est <- get.Q(y)
                }
                if(noise.sd.est.method=="MAD"){
                    noise.sd.est <- mad(res, na.rm=T)
                }
                if(noise.sd.est.method=="Qn"){
                    noise.sd.est <- Qn(na.omit(res))
                }
                if(noise.sd.est.method=="SD"){
                    noise.sd.est <- sd(res, na.rm=T)
                }
                noise.sd.est <- max(bound.noise.sd, noise.sd.est)
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
                   test.statistic=scarm.statistic, critvals=critvals, 
                   noise.sd=noise.sd, slope.diff=slope.diff,
                   tseries=tseries, right.width=r, min.left.width=ell.min, min.width=n.min, max.width=n.max, 
                   sign.level=sign.level, bound.noise.sd=bound.noise.sd, rtr=rtr, 
                   noise.sd.est.method=noise.sd.est.method, autocorrelations=autocorrelations)
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
  plot(x$tseries, type="l", main="Time series and SCARM signal extraction", xlab="time", ylab="", ...)
  lines(x$signal, col=2, lwd=1)
  legend(x="topleft", legend=c("Data", "SCARM signal extraction"),lty=c(1, 1),col=c("black", "red"),lwd=c(1, 1), bty="n")
  plot(x$adapted.width, type="l",  main="Adapted window widths", xlab="time", ylab="", ...)
}
