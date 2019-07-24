#' Conditional local FDR (CLfdr)
#'
#' CLfdr returns the local false discovery rate (FDR) conditional on auxiliary
#' covariate information
#'
#' In many multiple testing applications, the auxiliary information is widely
#' available and can be useful. Such information can be summary statistics
#' from a similar experiment or disease, the lengths of gene coding regions,
#' and minor allele frequencies of SNPs.
#'
#' \code{y} is a vector of *m* *z*-values, one of each hypothesis
#' under test. The *z*-values follow N(0,1) if their corresponding
#' null hypotheses are true. Other types of test statistics, such
#' as *t*-statistics and *p*-values can be transformed to
#' *z*-values. In practice, if the distribution of *z*-values is
#' far from N(0,1), recentering and rescaling of the *z*-values
#' may be necessary.
#'
#' \code{x} contains auxiliary covariate information. For a single covariate,
#' \code{x} should be a vector of length *m*. For multiple covariates,
#' \code{x} should be a matrix with *m* rows. The covariates can be either
#' continuous or ordered.
#'
#' \code{pi0.method} specifies the method used to estimate the overall true
#' null proportion. If the *z*-values are generated from the normal
#' means model, the "JC" method from Jin and Cai (2007) JASA can be
#' a good candidate. Otherwise, the right-boundary procedure ("RB",
#' Liang and Nettleton, 2012, JRSSB) is used.
#'
#' \code{bw} are bandwidth values for estimating local alternative density.
#' Suppose there are *p* covariates, then \code{bw} should be a
#' vector of *p*+1 positive numerical values. By default, these
#' bandwidth values are chosen by cross-validation to minimize
#' a certain error measure. However, finding the optimal bandwidth
#' values by cross-validation can be computationally intensive,
#' especially when *p* is not small. If good estimates of bandwidth
#' values are available, for example, from the analysis of a similar
#' dataset, the bandwidth values can be specified explicitly to save time.
#'
#' \code{reltol} specifies the relative convergence tolerance when choosing the
#' bandwidth values (\code{bw}). It will be passed on to
#' [stats::optim()]. For most analyses, the default value
#' of 1e-4 provides reasonably good results. A smaller value such as 1e-5
#' or 1e-6 could be used for further improvement at the cost of more
#' computation time.
#'
#' @author Kun Liang, \email{kun.liang@@uwaterloo.ca}
#' @param x covariates, could be a vector of length *m* or a matrix
#' with *m* rows.
#' @param y a vector of *z*-values of length *m*.
#' @param pval a vector of p-values of length *m*. The p-values are
#' only used to computed the overall true null proportion when
#' \code{pi0.method}="RB".
#' @param pi0.method method to estimate the overall true null proportion (pi0).
#' "RB" for the right-boundary procedure (Liang and Nettleton, 2012, JRSSB) or
#' "JC" (Jin and Cai, 2007, JASA).
#' @param bw.init initial values for bandwidth, optional. If not specified,
#' normal-reference rule will be used.
#' @param bw bandwidth values.
#' @param reltol relative tolerance in optim function.
#' @param n.subsample size of the subsample when esitmating bandwidth.
#' @param check.gam indicator to perform gam.check function on the
#' nonparametric fit.
#' @param k.gam tuning parameter for mgcv::gam.
#' @param info indicator to print out fitting information.
#' @return
#' \item{fdr}{a vector of local FDR estimates. fdr\[i\] is the posteiror
#' probability of the ith null hypothesis is true given all the data.
#' 1-fdr\[i\] is the posterior probability of being a signal (the corresponding
#' null hypothesis is false).}
#' \item{FDR}{a vector of FDR values (q-values), which can be used to control
#' FDR at a certain level by thresholding the FDR values.}
#' \item{pi0}{a vector of true null probability estimates. This contains the
#' prior probabilities of being null.}
#' \item{bw}{a vector of bandwidths for conditional alternative density
#' estimation}
#' \item{fit.gam}{an object of mgcv::gam}
#' @export
#' @md
#' @import stats
#' @importFrom graphics hist
#' @importFrom mgcv gam gam.check
#' @references Liang (2019), Empirical Bayes analysis of RNA sequencing
#' experiments with auxiliary information, to appear in Annals of Applied
#' Statistics
#' @examples
#' data(pso)
#' ind.nm <- is.na(pso$tval_mic)

CLfdr <- function(x, y, pval=NULL, pi0.method="RB", bw.init=NULL, bw=NULL,
                reltol=1e-4, n.subsample=NULL, check.gam=FALSE, k.gam=NULL,
                info=TRUE){

    if(pi0.method == "JC"){
        pi0v <- EstNullProp_JC(y)
    }else if(pi0.method == "RB"){
        if(is.null(pval)) pval <- 2*(1-pnorm(abs(y)))
        pi0v <- EstNullProp_RB(pval)
    }else{
        stop("pi0.method not recognized.")
    }

    if(info) cat("Overall pi0: ", pi0v, date(),"\n")

    m <- length(y)
    if(is.vector(x)){
        if(length(x)!=m) stop("length of x should match y")
    }
    if(is.matrix(x)){
        if(nrow(x)!=m) stop("length of x should match y")
    }
    f0 <- dnorm(y)

    # estimate local FDR without covariate
    fit<-density(y, from=min(y), to=max(y), n=1024)
    f <- approx(fit$x, fit$y, y, rule=2, ties="ordered")$y
    lfdr <- pi0v*f0/f

    # estimate common f1 by weighting of 1-fdr
    f1 <- get_f1(y, lfdr)
    # update pi0 with common f1
    res.up <- get_pi0v_f1(xdat=x, lfdr=lfdr, pi0v=pi0v, f0=f0, f1=f1,
                        check.gam = check.gam, k.gam=k.gam)
    pi0v <- res.up$pi0v
    lfdr <- res.up$lfdr

    # estimate local f1
    if(is.null(bw)) bw <- get_bandwidth(bw.init, x, y, info, n.subsample,
                                        reltol, lfdr, pi0v)
    f1 <- get_local_f1_density(xdat=x, y=y, bw=bw, lfdr=lfdr)
    res.up <- get_pi0v_f1(xdat=x, lfdr=lfdr, pi0v=pi0v, f0=f0, f1=f1,
                        check.gam = check.gam, k.gam=k.gam)

    return(list(fdr=res.up$lfdr, FDR=EstFDR(res.up$lfdr), pi0=res.up$pi0v,
                bw=bw, fit.gam=res.up$fit.gam))
}
