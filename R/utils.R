#' Right-boundary procedure
#'
#' True null proportion (pi_0) estimator of Liang and Nettleton (2012), JRSSB
#' @param pval vector of p-values
#' @param lambda.vec vector of lambda candidates (excluding 0 and 1)
#' @return the estimate of the overall true null proportion
#' @export
#' @examples
#' pval <- c(runif(900), rbeta(100, 1, 10))
#' EstNullProp_RB(pval)
EstNullProp_RB <- function(pval, lambda.vec = 0.05*seq_len(19)){
    #fixed lambda
    if(length(lambda.vec)==1){
        return(sum(pval>lambda.vec)/(1-lambda.vec))
    }
    if(length(lambda.vec)<1){
        stop("lambda.vec empty")
    }
    if(sum(lambda.vec>0 & lambda.vec<1) < length(lambda.vec)){
        stop("Values in lambda.vec should be between 0 and 1.")
    }

    m <- length(pval)
    B <- length(lambda.vec)+1   # number of bins
    bin <- c(-0.01, lambda.vec, 1)
    bin.counts <- hist(pval, breaks=bin, plot=FALSE)$counts
    # all but the last bin
    bin.counts <- bin.counts[-length(bin.counts)]
    R <- cumsum(bin.counts)
    tail.m0 <- (m-R+1)/(1-lambda.vec)
    temp <- tail.m0[2:(B-1)]-tail.m0[seq_len(B-2)]
    # index for selected lambda
    if(sum(temp >= 0)>0){
        index <- min((2:(B-1))[temp >= 0])
    }else{
        index <- B-1
    }

    return(tail.m0[index]/m)
}


#' FDR estimation
#'
#' False discovery rate (FDR) estimation from local FDR
#' @param fdr vector of local FDR
#' @return the estimate of the FDR
#' @export
#' @examples
#' lfdr <- c(runif(900), rbeta(100, 1, 10))
#' FDR <- EstFDR(lfdr)
#' sum(FDR<0.05)
EstFDR <- function(fdr){
    o <- order(fdr, decreasing = FALSE)
    FDR <- fdr
    FDR[o] <- cumsum(fdr[o])/seq_along(fdr)
    return(FDR)
}

# pi_0 estimator of Jin and Cai (2007), JASA
EstNullProp_JC <- function(x,u=0,sigma=1)
{
    # x is a vector
    # u is the mean
    # sigma is the standard deviation

    z  = (x - u)/sigma
    xi = c(0:100)/100
    tmax=sqrt(log(length(x)))
    tt=seq(0,tmax,0.1)

    epsest=NULL

    for (j in seq_along(tt)) {

        t=tt[j]
        f  = t*xi
        f  = exp(f^2/2)
        w  = (1 - abs(xi))
        co  = 0*xi
        for (i in seq_len(101)) {
            co[i] = mean(cos(t*xi[i]*z));
        }
        epshat = 1 - sum(w*f*co)/sum(w)
        epsest=c(epsest,epshat)
    }
    return(epsest=1-max(epsest))
}

# bandwidth with normal reference rule
get_normal_bw <- function(xdat, y){
    n <- length(y)
    d <- cbind(xdat, y)
    d_sd <- apply(d, MARGIN=2, sd)
    d_iqr <- apply(d, MARGIN=2, IQR)/1.349
    d_mad <- apply(d, MARGIN=2, mad)
    sigma <- pmin(d_sd, d_iqr, d_mad)
    bw <- 1.06*sigma*n^(-1.0/(4+ncol(d)))
    return(bw)
}

# optimal bandwidth for conditional alternative density estimation
get_bandwidth <- function(bw.init, x, y, info, n.subsample, reltol, lfdr, pi0v){
    m <- length(y)
    # initial bandwidth
    if(is.null(bw.init)){
        # normal-reference rule
        bw.init <- get_normal_bw(xdat=x, y=y)
    }
    if(info) cat("Initial bandwidth: ", bw.init, "\n")

    if(is.null(n.subsample) || n.subsample>m){
        index.subsample <- seq_len(m)
    }else{
        index.subsample <- sample.int(m, size=n.subsample, replace=FALSE)
        if(info) cat("Estimate bandwidth with ", n.subsample, " samples.\n")
    }

    res <- optim(log(bw.init), get_local_f1_loocv_ise_optim,
                control=list(reltol=reltol),
                xdat=x, y=y, lfdr=lfdr, pi0v=pi0v,
                index.subsample=index.subsample)

    #bandwidth
    bw <- exp(res$par)
    if(info) cat("Bandwidth: ", bw, date(), "\n")
    return(bw)
}


# estimate f1, weighted by 1-lfdr
get_f1 <- function(y, lfdr){
    W <- pmax(1-lfdr, 0)
    fit <- density(y, n=1000, weights=W/sum(W))
    f1 <- approx(fit$x, fit$y, y, rule=2, ties="ordered")$y
    return(f1)
}


# given xdat (covariate), f0, f1 and initial values of lfdr and pi0v,
# iteratively update pi0 to maximize the likelihood
get_pi0v_f1 <- function(xdat, lfdr, pi0v, f0, f1, check.gam=FALSE, k.gam=NULL){

    if(is.vector(xdat)){
        xd <- 1
        dat <- data.frame(x1=xdat,lfdr=lfdr)
    }else{# matrix
        xd <- ncol(xdat)
        colnames(xdat) <- paste0('x', seq_len(xd))
        dat <- data.frame(xdat,lfdr)
    }
    xVariables <- paste0('x', seq_len(xd), collapse=',')
    if(is.null(k.gam)){
        txtCmd <- paste0('mgcv::gam(lfdr~s(', xVariables, '), data=dat)')
    }else{
        txtCmd <- paste0('mgcv::gam(lfdr~s(', xVariables,
                        ', k=k.gam),  data=dat)')
    }

    res <- get_pi0v_f1_iter(dat, lfdr, pi0v, f0, f1, k.gam, txtCmd)

    if(check.gam) mgcv::gam.check(res$fit)

    return(list(pi0v=res$pi0v, lfdr=res$lfdr, likhood=res$likhood,
                iter=res$iter, fit.gam=res$fit))
}

# iteratively update pi0 to maximize the likelihood
get_pi0v_f1_iter <- function(dat, lfdr, pi0v, f0, f1, k.gam, txtCmd,
                            epsilon=1e-5, max.iter=200) {
    iter <- 1
    likhood <- sum(log(pi0v*f0+(1-pi0v)*f1))
    pi0.target <- mean(pi0v)

    while(1){
        pi0v.old <- pi0v
        lfdr.old <- lfdr
        likhood.old <- likhood

        # E step
        lfdr <- pi0v*f0/(pi0v*f0+(1-pi0v)*f1)
        dat$lfdr <- lfdr
        # M step
        fit <- eval(parse(text=txtCmd))

        pi0v <- fit$fitted.values
        pi0v <- pmax(pmin(pi0v, 1), 0)
        pi0v <- normalize_pi0(pi0v, pi0.target)

        likhood <- sum(log(pi0v*f0+(1-pi0v)*f1))
        dif <- likhood-likhood.old
        if(dif<0){
            pi0v <- pi0v.old
            likhood <- likhood.old
        }
        if(dif < epsilon) break
        if(iter==max.iter) break
        iter <- iter+1
    }
    return(list(pi0v=pi0v, lfdr=pi0v*f0/(pi0v*f0+(1-pi0v)*f1), likhood=likhood,
                iter=iter, fit.gam=fit))
}

# conditional density
get_local_f1_density <- function(xdat, y, bw, lfdr){
    bw.x <- bw[-length(bw)]
    bw.y <- bw[length(bw)]
    f1 <- rep(0, length(y))

    for(i in seq_along(y)){
        W <- get_kweight(xdat, bw.x, i)*(1-lfdr)
        fit <-density(y, bw=bw.y, weights=W/sum(W))

        f1[i] <- approx(fit$x, fit$y, y[i])$y
    }

    return(f1)
}

# leave-one-out cross-validation wrapper for optim function
# bw in log form so no restriction
get_local_f1_loocv_ise_optim <- function(bw, xdat, y, lfdr, pi0v,
                                        index.subsample=NULL){
    get_local_f1_loocv_ise(xdat=xdat, y=y,
                            bw.x=exp(bw[-length(bw)]), bw.y=exp(bw[length(bw)]),
                            lfdr=lfdr, pi0v=pi0v,
                            index.subsample=index.subsample)
}

# xdat: m*p matrix
get_local_f1_loocv_ise <- function(xdat, y, bw.x, bw.y, lfdr, pi0v,
                                    index.subsample=NULL){
    m <- length(y)

    # leave-one-out f1
    I1 <- f1 <- rep(0, length(index.subsample))

    tmp <- density(y, bw=bw.y)
    lb <- tmp$x[1]
    rb <- tmp$x[length(tmp$x)]
    bin.length.y <- mean(diff(tmp$x))

    for(j in seq_along(index.subsample)){
        i <- index.subsample[j]
        W <- get_kweight_loo(xdat, bw.x, i)*(1-lfdr[-i])
        fit <- density(y[-i], bw=bw.y, from=lb, to=rb, weights=W/sum(W))

        I1[j] <- sum(fit$y^2) * bin.length.y

        f1[j] <- approx(fit$x, fit$y, y[i])$y
    }
    I2 <- -2*f1*(1-lfdr[index.subsample])/pmax(1-pi0v[index.subsample],
                                            .Machine$double.xmin)
    return(mean(I1)+mean(I2))
}



# Gaussian kernel
kern_G <- function (x, xi, h)
{
    exp(-((xi-x)/h)^2/2)/sqrt(2*pi)
}

# product of kernel weight
get_kweight <- function(xdat, bw.x, i){
    # 1d
    if(length(bw.x)==1) return(kern_G(xdat[i], xdat, bw.x))
    # >1d
    W <- rep(1, nrow(xdat))
    for(j in seq_len(ncol(xdat))){
        W <- W*kern_G(xdat[i, j], xdat[, j], bw.x[j])
        W <- W/sum(W)
    }
    return(W)
}

# product of kernel weight leave-one-out
get_kweight_loo <- function(xdat, bw.x, i){
    # 1d
    if(length(bw.x)==1) return(kern_G(xdat[i], xdat[-i], bw.x))
    # >1d
    W <- rep(1, nrow(xdat)-1)
    for(j in seq_len(ncol(xdat))){
        W <- W*kern_G(xdat[i, j], xdat[-i, j], bw.x[j])
        W <- W/sum(W)
    }
    return(W)
}

normalize_pi0_logit <- function(delta, pi0v, pi0.target){
    mean(expit(logit(pi0v)+delta))-pi0.target
}

normalize_pi0 <- function(pi0v, pi0.target){
    delta <- pi0.target-mean(pi0v)
    if(delta==0) return(pi0v)
    if(sum(pi0v+delta>1)>0 | sum(pi0v+delta<0)>0){
        bound <- 1
        if(delta > 0){
            while(1){
                if(normalize_pi0_logit(bound, pi0v, pi0.target)>0){
                    delta.logit <- uniroot(normalize_pi0_logit, c(0,bound),
                                            pi0v=pi0v,
                                            pi0.target=pi0.target)$root
                    break
                }else{
                    bound <- bound * 10
                }
            }
        }else{
            while(1){
                if(normalize_pi0_logit(-bound, pi0v, pi0.target)<0){
                    delta.logit <- uniroot(normalize_pi0_logit, c(-bound,0),
                                            pi0v=pi0v,
                                            pi0.target=pi0.target)$root
                    break
                }else{
                    bound <- bound * 10
                }
            }
        }
        return(expit(logit(pi0v)+delta.logit))
    }else{
        return(pi0v+delta)
    }
}

logit <- function(x){
    log(x/(1-x))
}

expit <- function(x, bound=700){
    res <- exp(x)/(1+exp(x))
    res[x>bound] <- 1
    return(res)
}
