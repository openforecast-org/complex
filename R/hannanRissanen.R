# Hannan-Rissanen estimation of complex ARMA(p,q) (Hannan & Rissanen, 1982), used for the
# starting values in clm(fast=TRUE) and for the screening in auto.carima(fast=TRUE).
# 1. A long AR(m) is estimated via complex least squares, and its residuals approximate the errors.
# 2. The series is regressed on its p lags and q lags of the approximate errors (and the constant),
#    via complex least squares, which gives AR and MA parameters jointly.
# The estimates are consistent, but not efficient.

# Complex least squares: minimises the sum of squared moduli of the residuals
complexLS <- function(X, target){
    XH <- t(Conj(X));
    return(as.vector(solve(XH %*% X, XH %*% target)));
}

# Matrix of lags 1..k of the vector x (rows correspond to the elements of x, NA where unavailable)
lagMatrix <- function(x, k){
    n <- length(x);
    X <- matrix(NA_complex_, n, k);
    for(j in seq_len(k)){
        X[(j+1):n, j] <- x[seq_len(n-j)];
    }
    return(X);
}

# Residuals of ARMA(p,q) with the constant for the series u (all values available).
# The first p residuals are NA, the errors before the start are set to zero.
armaResiduals <- function(u, constant, ar, ma){
    n <- length(u);
    p <- length(ar);
    q <- length(ma);
    errors <- rep(NA_complex_, n);
    errorsLagged <- rep(0+0i, q);
    constantValue <- if(is.null(constant)) 0 else constant;
    for(t in (p+1):n){
        fitted <- constantValue;
        if(p>0){
            fitted <- fitted + sum(ar * u[t-seq_len(p)]);
        }
        if(q>0){
            fitted <- fitted + sum(ma * errorsLagged);
        }
        errors[t] <- u[t] - fitted;
        if(q>0){
            errorsLagged <- c(errors[t], errorsLagged[-q]);
        }
    }
    return(errors);
}

# Approximate errors from a long AR (stage 1). m is the order of the long AR.
hrErrors <- function(u, constant=FALSE, m=NULL){
    n <- length(u);
    if(is.null(m)){
        m <- min(ceiling(10*log10(n)), floor((n-1)/4));
    }
    X <- cbind(if(constant) rep(1+0i, n), lagMatrix(u, m));
    rows <- (m+1):n;
    coefs <- complexLS(X[rows,,drop=FALSE], u[rows]);
    errors <- rep(NA_complex_, n);
    errors[rows] <- u[rows] - X[rows,,drop=FALSE] %*% coefs;
    return(errors);
}

# Stage 2: parameters of ARMA(p,q) with or without constant, given the approximate errors
hannanRissanen <- function(u, p, q, constant=FALSE, errors=NULL){
    n <- length(u);
    if(p==0 && q==0){
        return(list(constant=if(constant) mean(u) else NULL, ar=NULL, ma=NULL));
    }
    if(q>0 && is.null(errors)){
        errors <- hrErrors(u, constant);
    }
    X <- cbind(if(constant) rep(1+0i, n), lagMatrix(u, p), if(q>0) lagMatrix(errors, q));
    rows <- which(apply(!is.na(X), 1, all) & !is.na(u));
    coefs <- complexLS(X[rows,,drop=FALSE], u[rows]);
    k <- as.integer(constant);
    return(list(constant=if(constant) coefs[1] else NULL,
                ar=coefs[k+seq_len(p)],
                ma=coefs[k+p+seq_len(q)]));
}

# Shrink AR / MA parameters towards zero until the roots of the polynomials are outside the unit circle
hrAdmissible <- function(parameters, sign){
    if(length(parameters)==0){
        return(parameters);
    }
    for(i in 1:100){
        if(all(Mod(polyroot(c(1, sign*parameters)))>1)){
            break;
        }
        parameters <- parameters*0.9;
    }
    return(parameters);
}
