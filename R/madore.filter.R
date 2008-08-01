madore.filter <- function(Y, byrow=FALSE, min.width=20, max.width=200, 
                             start.width=min.width, test.sample.size=min.width/2,
                             width.search="geometric", 
                             rtr.size=10, extraction.lag=0, 
                             NA.sample.size=10, minNonNAs=5){

################
# Preparations #
################
    min.width <- round(min.width)
    max.width <- round(max.width)
    test.sample.size <- round(test.sample.size)
    rtr.size <- round(rtr.size)
    extraction.lag <- round(extraction.lag)
    NA.sample.size <- round(NA.sample.size)
    minNonNAs <- round(minNonNAs)

##############################
# Stopping and Warning rules #
##############################
    if(missing(Y))
        stop("Input data set is missing with no default.\n")

    if(dim(Y)[1] < min.width)
        stop("Length of data set must be at least 'min.width' =", min.width, ".\n")

    if(!is.logical(byrow))
        stop("'byrow' must be either 'TRUE' or 'FALSE'.\n")

    if(!is.numeric(min.width))
        stop("'min.width' must be numeric.\n")

    if(!is.numeric(max.width))
        stop("'max.width' must be numeric.\n")

    if(max.width <= min.width)
        stop("'max.width' must be greater than 'min.width'.\n")

    if(min.width < 10){
        min.width <- 10
        warning("'min.width' must be at least 10; 'min.width' is set to 10.\n")
        }

    if(!is.numeric(start.width))
        stop("'start.width' must be numeric.\n")

    if(start.width < min.width){
        start.width <- min.width
        warning("'start.width' must be at least 'min.width'=", min.width,"; 'start.width' is set to 'min.width'=", min.width,".\n")
        }

    if(start.width > max.width){
        start.width <- max.width
        warning("'start.width' must be at most 'max.width'=", max.width,"; 'start.width' is set to 'max.width'=", max.width,".\n")
        }

    if(!is.numeric(test.sample.size))
        stop("'test.sample.size' must be numeric.\n")

    if(test.sample.size < 5){
        test.sample.size <- 5
        warning("'test.sample.size' must be at least 5; 'test.sample.size' is set to 5.\n")
        }

    if(test.sample.size > min.width/2)
        warning("'test.sample.size' is greater than 'min.width'/2. If 'test.sample.size' is greater than n(t)/2,
        where n(t) is the currently used window width, 'test.sample.size' is set to floor(n(t)/2).\n")

    if(width.search!="linear" & width.search!="binary" & width.search!="geometric")
        stop("'width.search' must be a character, either 'linear', 'binary' or 'geometric'.\n")

    if(!is.numeric(rtr.size))
        stop("'rtr.size' must be numeric.\n")

    if(rtr.size < 0){
        rtr.size <- 0
        warning("'rtr.size' must be at least 0; 'rtr.size' is set to 0.\n")
        }

    if(rtr.size > min.width)
        warning("'rtr.size' is greater than 'min.width'. If 'rtr.size' is greater than n(t),
        where n(t) is the currently used window width, 'rtr.size' is set to n(t).\n")

    if(!is.numeric(extraction.lag))
        stop("'extraction.lag' must be numeric.\n")

    if(extraction.lag < 0){
        extraction.lag <- 0
        warning("'extraction.lag' must be at least 0; 'extraction.lag' is set to 0.\n")
    }

    if(extraction.lag > min.width/2){
        extraction.lag <- min.width/2
        warning("'extraction.lag' must not be greater than 'min.width'/2; 'extraction.lag' is set to 'min.width'/2.\n")
    }

    if(!is.numeric(NA.sample.size))
        stop("'NA.sample.size' must be numeric.\n")

    if(NA.sample.size < 10){
        NA.sample.size <- 10
        warning("'NA.sample.size' must be at least 10; 'NA.sample.size' is set to 10.\n")
        }

    if(NA.sample.size > min.width){
        NA.sample.size <- min.width
        warning("'NA.sample.size' must not be greater than 'min.width'; 'NA.sample.size' is set to 'min.width'.\n")
    }

    if(!is.numeric(minNonNAs))
        stop("'minNonNAs' must be numeric.\n")

    if(minNonNAs < 5){
        minNonNAs <- 5
        warning("'minNonNAs' must be at least 5; 'minNonNAs' is set to 5.\n")
    }

    if(minNonNAs > NA.sample.size){
        minNonNAs <- NA.sample.size
        warning("'minNonNAs' must not be greater than 'NA.sample.size'; 'minNonNAs' is set to 'NA.sample.size'.\n")
    }
    
    if(max.width > 100)
        warning("For window widths > 100 the correction factor 'cqn' for the Qn scale estimator is chosen approximately.")

######################
# Internal functions #
######################

    data(critvals)
    data(const)

    oRM <- function(y){
        n <- length(y)
        s <- 1:n
        medquotient.t <- c()
        for (k in s) {
            medquotient.t <- c(medquotient.t, median((y[k] - y[s[s!=k]]) / (k - s[s!=k]) ,na.rm=TRUE))
        }
        beta.oRM <- median(medquotient.t, na.rm=TRUE)
        mu.oRM   <- median(y - beta.oRM*(s-n), na.rm=TRUE)
        y.hat <- mu.oRM + beta.oRM*(s-n)
        f <- which(is.na(y))
        if(length(f) > 0){
            y[f] <- y.hat[f]
        }
        if(extraction.lag>0){
            mu.oRM <- y.hat[n-extraction.lag+1]
        }
        res <- y - y.hat
        list(mu.oRM=mu.oRM, y.hat=y.hat, res=res)
    }

    get.test.sample.size <- function(n, test.sample.size){
        if(test.sample.size > n/2) return(floor(n/2))
        if(test.sample.size <= n/2) return(test.sample.size)
    }

    get.critval <- function(n, test.sample.size, critvals){
        if(n <= 600) return(critvals[n, test.sample.size])  else
        return(2 * qhyper(1 - 0.1/2, floor(n/2), floor(n/2), test.sample.size) - test.sample.size)
    }

    get.TS <- function(res){
        return(abs(sum(sign(res), na.rm = TRUE)))
    }

    ww.adaption <- function(x){
        n <- length(x)
        if(length(which(!is.na(x[(n-NA.sample.size+1):n]))) < minNonNAs){
            n.a <- NA
        } else {
            oRM.extr <- oRM(x)
            test.sample.size2  <- get.test.sample.size(n, test.sample.size)
            test.index <- (n - test.sample.size2 + 1):n
            TS <- get.TS(oRM.extr$res[test.index])
            if(TS > get.critval(n, test.sample.size2, critvals)){
                if(n > min.width){
                    direction <- -1
                } else {
                direction <- 0
                }
            } else {
                direction <- 0
            }
            if(direction==-1){
                if(width.search=="linear"){
                    i=1
                    n.a <- n
                    repeat{
                        if(n == min.width) break
                        n.a <- n.a-i
                        x.win <- x[(n-n.a+1):n]
                        oRM.extr <- oRM(x.win)
                        test.sample.size2  <- get.test.sample.size(n.a, test.sample.size)
                        test.index <- (n.a - test.sample.size2 + 1):n.a
                        TS <- get.TS(oRM.extr$res[test.index])
                        if(TS <= get.critval(n.a, test.sample.size2, critvals)) break
                    }
                }
                if(width.search=="binary"){
                    n.l <- min.width
                    n.u <- n
                    n.a <- n.l + floor((n.u-n.l)/2)
                }
                if(width.search=="geometric"){
                    i=1
                    repeat{
                        n.l <- n-(2^i)
                        if(n.l<=min.width){
                            n.l <- min.width
                            n.u <- n-2^(i-1)
                            n.a <- n.l + floor((n.u-n.l)/2)
                            break
                        }
                        x.win <- x[(n-n.l+1):n]
                        oRM.extr <- oRM(x.win)
                        test.sample.size2  <- get.test.sample.size(n.l, test.sample.size)
                        test.index <- (n.l - test.sample.size2 + 1):n.l
                        TS <- get.TS(oRM.extr$res[test.index])
                        if(TS > get.critval(n.l, test.sample.size2, critvals)){
                            i <- i+1
                        } else {
                            n.u <- n-2^(i-1)
                            n.a <- n.l + floor((n.u-n.l)/2)
                            break
                        }
                    }
                }
                if(width.search!="linear"){
                    repeat{
                        x.win <- x[(n-n.a+1):n]
                        oRM.extr <- oRM(x.win)
                        test.sample.size2  <- get.test.sample.size(n.a, test.sample.size)
                        test.index <- (n.a - test.sample.size2 + 1):n.a
                        TS <- get.TS(oRM.extr$res[test.index])
                        if(n.a==n.l){
                            break
                        }
                        if(n.a==n.u & TS <= get.critval(n.a, test.sample.size2, critvals)){
                            break
                        }
                        if(n.a==n.u & TS > get.critval(n.a, test.sample.size2, critvals)){
                            n.a <- n.l
                            break
                        }
                        if(TS > get.critval(n.a, test.sample.size2, critvals)){
                            n.u <- n.a
                            n.a <- n.l + floor((n.u-n.l)/2)
                        } else {
                            n.l <- n.a
                            n.a <- n.l + ceiling((n.u-n.l)/2)
                        }
                    }
                }
            }
            if(direction==0){
                n.a <- n
            }
        }
        return(n.a)
    }

    OGK.qn <- function(x, cqn){
        trimstd <- function(x, alphahalbe){
            n <- length(x)
            k <- floor(n*alphahalbe)
            trim <- sqrt(sum(sort(x^2)[(k+1):(n-k)])/(n-(2*k)-1))
            return(trim)
        }
        gnanakett.qn <- function(x, cqn){
            if(is.vector(x) == TRUE){
                k <- length(x)
            } else {
                k <- length(x[1,])
            }
            umat <- diag(rep(1,k))
            for (i in 1:(k-1)){
                for (j in (i+1):k){
                    umat[i,j] <- umat[j,i] <- 0.25*((Qn((x[,i]+x[,j]),cqn))^2-(Qn((x[,i]-x[,j]),cqn))^2)
                }
            }
            return(umat)
        }
        dia        <- apply(x, 2, Qn, cqn)
        anders     <- which(dia < 0.002)
        dia[anders] <- apply(as.matrix(x[, anders]), 2, trimstd, 0.1)
        extra      <- which(dia==0)
        dia[extra] <- 10
        dmat       <- diag(dia)
        inv.d      <- solve(dmat)
        inv.d[extra,extra] <- 0
        y          <- x %*% inv.d
        rmat       <- gnanakett.qn(y, cqn)
        emat       <- eigen(rmat)$vectors
        dia[extra] <- 0.002
        dmat       <- diag(dia)
        amat       <- dmat %*% emat
        z          <- y %*% emat
        diaz       <- apply(z, 2, Qn, cqn)
        zanders    <- which(diaz < 0.002)
        diaz[zanders] <- apply(as.matrix(z[,zanders]), 2, trimstd, 0.1)
        diaz[diaz<=0.002] <- 0.002
        Gamma      <- diag(diaz^2)
        vmat       <- amat %*% Gamma %*%t(amat)
        return(vmat)
    }

    LS.reg <- function(x, res, dat, cqn){
        di.qn  <- mahalanobis(res, center=rep(0,dim(dat)[2]), cov=OGK.qn(res, cqn), tol=2.220446e-16)
        dnqn   <- qchisq(0.95, dim(dat)[2]) * median(di.qn)/ qchisq(0.5, dim(dat)[2])
        kov    <- cov(cbind(x,dat)[di.qn <= dnqn,])
        beta   <- kov[1,2:(dim(dat)[2]+1)] / kov[1,1]
        alpha  <- apply(dat[di.qn <= dnqn,], 2, mean) - (beta * mean(x))
        hat.y  <- alpha + beta*(length(x)-extraction.lag+1)
        return(hat.y)
    }

####################
# Internal objects #
####################
    library(robustbase)
    if(byrow==TRUE) Y <- t(Y)
    T <- dim(Y)[1]
    k <- dim(Y)[2]
    y <- matrix(NA,T,k)
    for(i in 1:k){
        y[,i] <- as.vector(Y[,i])
    }
    Y <- y
    if(!is.numeric(Y))
        stop("Data set must be numeric.\n")
    n <- round(start.width)
    signals <- matrix(NA, nrow=T, ncol=k)
    widths <- matrix(NA, nrow=T, ncol=k)
    ov.width <- rep(NA,T)

###################
# Main  Algorithm #
###################
    for(t in n:T){
        Y.win <- Y[(t-n+1):t,]
        nk <- apply(Y.win,2,ww.adaption)
        widths[t,] <- nk
        if(any(!is.na(nk))){
            ov.width[t] <- min(nk,na.rm=TRUE)
            position <- which(!is.na(nk))
            Y.win2 <- Y[(t-ov.width[t]+1):t,position]
            kt <- length(position)
            if(kt == 1){
                oRM.extr <- oRM(Y.win2)
                signals[t,position] <- c(oRM.extr$mu.oRM)
            }
            if(kt > 1){
                oRM.extr <- apply(Y.win2, 2, oRM)
                if(any(is.na(Y.win2))){
                    for(i in 1:kt){
                        j <- which(is.na(Y.win2[,i]))
                        if(length(j) > 0){
                            Y.win2[j,i] <- as.vector(oRM.extr[[i]]$y)[j]
                        }
                    }
                }
                res <- matrix(NA, ncol = kt, nrow = ov.width[t])
                for(i in 1:kt){
                    res[,i] <- oRM.extr[[i]]$res
                }
                x <- 1:ov.width[t]
                cqn <- const$const[which.min(abs(const$n - ov.width[t]))]
                signals[t,position] <- LS.reg(x, res, Y.win2, cqn)
                if(rtr.size>0){
                    if(rtr.size > ov.width[t]){
                        for(i in 1:kt){
                            if(signals[t,position[i]] > max(Y.win[,position[i]], na.rm = TRUE)){
                                signals[t,position[i]] <- max(Y.win[,position[i]], na.rm = TRUE)
                            }
                            if(signals[t,position[i]] < min(Y.win[,position[i]], na.rm = TRUE)){
                                signals[t,position[i]] <- min(Y.win[,position[i]], na.rm = TRUE)
                            }
                        }
                    } else {
                        for(i in 1:kt){
                            if(all(is.na(Y.win[(n-rtr.size+1):n,i]))){
                                if(signals[t,position[i]] > max(Y.win[,position[i]], na.rm = TRUE)){
                                    signals[t,position[i]] <- max(Y.win[,position[i]], na.rm = TRUE)
                                }
                                if(signals[t,position[i]] < min(Y.win[,position[i]], na.rm = TRUE)){
                                    signals[t,position[i]] <- min(Y.win[,position[i]], na.rm = TRUE)
                                }
                            } else {
                                if(signals[t,position[i]] > max(Y.win[(n-rtr.size+1):n,position[i]], na.rm = TRUE)){
                                    signals[t,position[i]] <- max(Y.win[(n-rtr.size+1):n,position[i]], na.rm = TRUE)
                                }
                                if(signals[t,position[i]] < min(Y.win[(n-rtr.size+1):n,position[i]], na.rm = TRUE)){
                                    signals[t,position[i]] <- min(Y.win[(n-rtr.size+1):n,position[i]], na.rm = TRUE)
                                }
                            }
                        }
                    }
                }
            }
            if(ov.width[t] < max.width){
                n <- ov.width[t]+1
            }
        }
        if(is.na(ov.width[t])){
            n <- min.width
        }
    }
    if(extraction.lag>0){
        signals <- rbind(matrix(NA,ncol=k,nrow=min.width-extraction.lag-1), signals[min.width:T,],
            matrix(NA,ncol=k,nrow=min.width-extraction.lag))
    }
    result <- list(signals=signals, widths=widths, overall.width=ov.width, Y=Y, byrow=byrow, min.width=min.width,
        max.width=max.width, start.width=start.width, test.sample.size=test.sample.size, width.search=width.search,
        rtr.size=rtr.size, extraction.lag=extraction.lag, NA.sample.size=NA.sample.size, minNonNAs=minNonNAs)
    return(structure(result, class="madore.filter"))
}

##################
# Default output #
##################
print.madore.filter <- function(x, ...){
    N <- dim(x$signals)[1]
    cat('$signals \n')
    if(N <= 100){
        print(x$signals, ...)
    } else {
        print(x$signals[1:(x$min.width + 9),])
        cat('Only the first 10 signal estimations are printed.\n')
    }
}

################
# Default plot #
################
plot.madore.filter <- function(x, ...){
    T <- dim(x$Y)[1]
    k <- dim(x$Y)[2]
    plot(rep(NA,T), main="Multivariate time series and signal extraction by madore.filter",
    type='l', xlab='time', ylab='', ylim=c(min(x$Y, na.rm=TRUE),max(x$Y, na.rm=TRUE)), las=1,...)
    for(i in 1:k){
        lines(x$Y[,i])
        lines(x$signals[,i], col=2, lwd=2)
    }
    legend(x="topleft", bty="n", legend=c("Time Series", "Filtered Signal"), lty=c(1,1), col=c(1,2), lwd=c(1,1))
}
