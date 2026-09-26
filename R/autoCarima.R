#' Automatic Complex ARIMA
#'
#' Function selects the orders of complex ARIMA (see \link[complex]{carima}) based on an
#' information criterion.
#'
#' The selection follows the principles of \link[smooth]{auto.ssarima}: the AR and MA orders are
#' selected for each order of differencing, and the best model across all of them is chosen.
#' With \code{search="stepwise"}:
#' \enumerate{
#' \item For each \eqn{d} from 0 to \code{orders$i}, the AR and MA orders are selected via a
#' neighbourhood search (similar to Hyndman & Khandakar, 2008): starting from cARIMA(0,d,0),
#' all models with p and q different by at most one are fitted, and the search moves to the best
#' of them, as long as the information criterion improves. If \code{constant=NULL}, the constant
#' is included for \eqn{d=0} and excluded otherwise;
#' \item If \code{constant=NULL}, the best model of each \eqn{d} is re-estimated with the opposite
#' option for the constant and the neighbourhood search continues from it. The best of all the
#' estimated models is returned.
#' }
#' With \code{search="full"}, all combinations of p, d, q and the constant are fitted
#' (except for cARIMA(p,0,q) without constant when \code{constant=NULL}).
#'
#' Models with different orders of differencing and AR use different numbers of initial values,
#' which \link[complex]{clm} constructs from the data. To make the information criteria
#' comparable, all models are estimated on the whole sample, but their criteria are calculated on
#' the common sample without the first \eqn{k} observations, where \eqn{k} is the sum of the
#' maximum orders. The log-likelihood on that sample is calculated at the estimated parameters,
#' with the scale re-estimated.
#'
#' Every model is estimated from the default starting values and from the parameters of the
#' closest model already estimated (with zeros for the new parameters), and the one with the
#' higher likelihood is kept.
#'
#' @template author
#' @template keywords
#'
#' @param y vector or ts object with the complex-valued time series.
#' @param orders list of maximum orders \code{list(ar=p, i=d, ma=q)} (or vector \code{c(p,d,q)}).
#' @param constant if \code{NULL}, the function checks whether the constant is needed.
#' Otherwise \code{TRUE} or \code{FALSE}.
#' @param ic the information criterion used in the selection.
#' @param search \code{"stepwise"} for the neighbourhood search or \code{"full"} for all the
#' combinations of orders.
#' @param distribution the distribution used in the likelihood, see \link[complex]{clm}.
#' @param h the forecast horizon.
#' @param holdout if \code{TRUE}, the last \code{h} observations are held out.
#' @param silent if \code{FALSE}, the progress is printed and the final model is plotted.
#' @param ... other parameters passed to \link[complex]{carima} and \link[complex]{clm}.
#'
#' @return The selected model of the class \code{carima} (see \link[complex]{carima}), with
#' additional elements:
#' \itemize{
#' \item \code{ICs} - data frame with the orders, the constant and the information criterion
#' (on the common sample) of all the models estimated in the process,
#' \item \code{ic} - the information criterion used,
#' \item \code{icObsDropped} - the number of first observations dropped for the criteria.
#' }
#'
#' @references \itemize{
#' \item Svetunkov, I., Boylan, J. E. (2019). State-space ARIMA for supply-chain forecasting.
#' International Journal of Production Research, 58(3), 818-827.
#' \item Hyndman, R. J., Khandakar, Y. (2008). Automatic time series forecasting: the forecast
#' package for R. Journal of Statistical Software, 27(3), 1-22.
#' }
#'
#' @seealso \code{\link[complex]{carima}}, \code{\link[smooth]{auto.ssarima}},
#' \code{\link[smooth]{auto.msarima}}
#'
#' @examples
#' set.seed(41)
#' e <- rcnorm(300, 0, sigma2=1, varsigma2=0.3+0.2i)
#' y <- cumsum(e + c(0, (0.5-0.2i)*e[-300]))
#' model <- auto.carima(y, orders=list(ar=2, i=1, ma=2))
#' model
#' model$ICs
#'
#' @rdname auto.carima
#' @export
auto.carima <- function(y, orders=list(ar=3, i=2, ma=3), constant=NULL,
                        ic=c("AICc","AIC","BIC","BICc"), search=c("stepwise","full"),
                        distribution=c("dcnorm","dcgnorm"),
                        h=0, holdout=FALSE, silent=TRUE, ...){
    cl <- match.call();
    ic <- match.arg(ic);
    search <- match.arg(search);
    distribution <- match.arg(distribution);
    ellipsis <- list(...);
    # The starting values are set by the function
    ellipsis$B <- NULL;

    orders <- carimaOrders(orders);
    arMax <- orders$ar;
    iMax <- orders$i;
    maMax <- orders$ma;

    y <- as.vector(y);
    if(!is.complex(y)){
        y <- as.complex(y);
    }
    obsAll <- length(y);
    holdoutData <- NULL;
    if(holdout){
        if(h<=0){
            stop("holdout=TRUE requires h>0.", call.=FALSE);
        }
        holdoutData <- y[obsAll-h+seq_len(h)];
        y <- y[seq_len(obsAll-h)];
    }

    # Observations dropped for the calculation of the information criteria
    obsDropped <- arMax + iMax + maMax;
    if(length(y)-obsDropped <= 2*(arMax+maMax+1)+3){
        stop("The sample is too small for the maximum orders provided.", call.=FALSE);
    }

    # The information criterion on the common sample
    icTrimmed <- function(model){
        errors <- residuals(model)[-seq_len(obsDropped)];
        obsTrimmed <- length(errors);
        if(distribution=="dcnorm"){
            errorsMatrix <- complex2vec(errors);
            scaleMatrix <- t(errorsMatrix) %*% errorsMatrix / obsTrimmed;
            llikelihood <- -obsTrimmed*(log(2*pi) + 1 + 0.5*log(det(scaleMatrix)));
        }
        else{
            llikelihood <- cgnormConcentrated(errors, model$other$shape)$logLik;
        }
        return(clmIC(llikelihood, nparam(model), obsTrimmed,
                     isTRUE(model$other$shapeEstimated), ic));
    }

    # Starting values from an estimated model for cARIMA(p,d,q) with the given constant
    startFrom <- function(model, p, q, constant){
        coefs <- coef(model);
        pFrom <- model$orders$ar;
        qFrom <- model$orders$ma;
        k <- as.integer(model$constant);
        start <- c(if(constant) {if(model$constant) coefs[1] else 0+0i},
                   coefs[k+seq_len(min(p, pFrom))], rep(0+0i, max(p-pFrom, 0)),
                   coefs[k+pFrom+seq_len(min(q, qFrom))], rep(0+0i, max(q-qFrom, 0)));
        if(distribution=="dcgnorm" && isTRUE(model$other$shapeEstimated)){
            start <- c(start, model$other$shape);
        }
        return(unname(start));
    }

    # Estimated models: list of list(model, IC), named by the orders and the constant
    estimated <- list();
    modelKey <- function(p, d, q, constant){
        paste(p, d, q, constant);
    }
    fitCandidate <- function(p, d, q, constant, startModels=list()){
        key <- modelKey(p, d, q, constant);
        if(!is.null(estimated[[key]])){
            return(estimated[[key]]);
        }
        fit <- function(B=NULL){
            carimaArgs <- c(list(y=y, orders=c(p, d, q), constant=constant,
                                 loss="likelihood", distribution=distribution),
                            if(!is.null(B)) list(B=B), ellipsis);
            tryCatch(do.call(carima, carimaArgs), error=function(e) NULL);
        }
        fits <- c(list(fit()),
                  lapply(startModels, function(model) fit(startFrom(model, p, q, constant))));
        fits <- fits[!sapply(fits, is.null)];
        if(length(fits)==0){
            result <- list(model=NULL, IC=Inf);
        }
        else{
            model <- fits[[which.max(sapply(fits, function(x) as.numeric(logLik(x))))]];
            ICValue <- icTrimmed(model);
            result <- list(model=model, IC=ifelse(is.finite(ICValue), ICValue, Inf));
        }
        estimated[[key]] <<- result;
        if(!silent){
            cat("cARIMA(", p, ",", d, ",", q, ")", ifelse(constant, " with constant", ""),
                ": ", round(result$IC, 3), "\n", sep="");
        }
        return(result);
    }

    constantOptions <- if(is.null(constant)) c(FALSE, TRUE) else constant;

    if(search=="full"){
        for(d in 0:iMax){
            for(constantValue in constantOptions){
                if(d==0 && !constantValue && is.null(constant)){
                    next;
                }
                # By level p+q, so that each model can start from the nested ones
                for(level in 0:(arMax+maMax)){
                    for(p in 0:arMax){
                        q <- level - p;
                        if(q<0 || q>maMax){
                            next;
                        }
                        startModels <- list();
                        if(p>0){
                            startModels$ar <- estimated[[modelKey(p-1, d, q, constantValue)]]$model;
                        }
                        if(q>0){
                            startModels$ma <- estimated[[modelKey(p, d, q-1, constantValue)]]$model;
                        }
                        fitCandidate(p, d, q, constantValue, startModels[!sapply(startModels, is.null)]);
                    }
                }
            }
        }
    }
    else{
        # Neighbourhood search: from the given model, move to the best of the models with p and q
        # different by at most one, as long as the information criterion improves
        neighbourhoodSearch <- function(best, d, constantValue){
            pBest <- best$model$orders$ar;
            qBest <- best$model$orders$ma;
            repeat{
                neighbours <- expand.grid(p=pBest+(-1:1), q=qBest+(-1:1));
                neighbours <- neighbours[neighbours$p>=0 & neighbours$p<=arMax &
                                             neighbours$q>=0 & neighbours$q<=maMax &
                                             !(neighbours$p==pBest & neighbours$q==qBest),,drop=FALSE];
                if(nrow(neighbours)==0){
                    break;
                }
                results <- lapply(seq_len(nrow(neighbours)), function(i){
                    fitCandidate(neighbours$p[i], d, neighbours$q[i], constantValue,
                                 if(!is.null(best$model)) list(best$model) else list());
                });
                neighbourICs <- sapply(results, `[[`, "IC");
                if(min(neighbourICs) < best$IC){
                    best <- results[[which.min(neighbourICs)]];
                    pBest <- neighbours$p[which.min(neighbourICs)];
                    qBest <- neighbours$q[which.min(neighbourICs)];
                }
                else{
                    break;
                }
            }
            return(best);
        }

        # Search for p and q within each d, starting from cARIMA(0,d,0). The constant is included
        # for d=0 and excluded for d>0 (unless it is provided), and is tested at the end.
        bestPerDiff <- vector("list", iMax+1);
        for(d in 0:iMax){
            constantValue <- if(is.null(constant)) (d==0) else constant;
            if(!silent){
                cat("Searching AR and MA orders with d=", d, "...\n", sep="");
            }
            best <- fitCandidate(0, d, 0, constantValue);
            if(!is.null(best$model)){
                best <- neighbourhoodSearch(best, d, constantValue);
            }
            bestPerDiff[[d+1]] <- best;
        }

        # Test the constant on the best model of each d, and continue the search from there:
        # e.g. a random walk with drift is otherwise missed, when the search without constant
        # compensates the drift with MA terms
        if(is.null(constant)){
            if(!silent){
                cat("Testing the constant...\n");
            }
            for(d in 0:iMax){
                bestPerD <- bestPerDiff[[d+1]];
                if(!is.null(bestPerD$model)){
                    constantValue <- !bestPerD$model$constant;
                    flipped <- fitCandidate(bestPerD$model$orders$ar, d, bestPerD$model$orders$ma,
                                            constantValue, list(bestPerD$model));
                    if(!is.null(flipped$model)){
                        neighbourhoodSearch(flipped, d, constantValue);
                    }
                }
            }
        }
    }

    # The best model among all the estimated ones
    ICs <- sapply(estimated, `[[`, "IC");
    if(all(!is.finite(ICs))){
        stop("None of the models could be estimated.", call.=FALSE);
    }
    bestModel <- estimated[[which.min(ICs)]]$model;

    ICsTable <- do.call(rbind, lapply(strsplit(names(estimated), " "), function(x){
        data.frame(p=as.integer(x[1]), d=as.integer(x[2]), q=as.integer(x[3]),
                   constant=as.logical(x[4]));
    }));
    ICsTable[[ic]] <- ICs;
    ICsTable <- ICsTable[order(ICsTable[[ic]]),];
    rownames(ICsTable) <- NULL;

    bestModel$ICs <- ICsTable;
    bestModel$ic <- ic;
    bestModel$icObsDropped <- obsDropped;
    bestModel$autoCall <- cl;
    bestModel <- carimaFinalise(bestModel, h, holdout, holdoutData, silent);

    return(bestModel);
}
