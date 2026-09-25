#' Complex ARIMA
#'
#' Function constructs complex ARIMA(p,d,q) model for a complex-valued time series,
#' without explanatory variables. It is a wrapper of \link[complex]{clm}.
#'
#' The model is
#' \deqn{(1-\phi_1 B-\dots-\phi_p B^p)(1-B)^d y_t = c + (1+\theta_1 B+\dots+\theta_q B^q)\epsilon_t ,}
#' where all the parameters are complex and \eqn{c} is the constant (included only when
#' \code{constant} is not \code{FALSE}). With \code{d>0} the constant acts as a drift.
#'
#' Only the non-seasonal model is supported. Parameters are estimated via \link[complex]{clm},
#' so its restrictions apply: the first lags are extrapolated and the first MA errors are set to
#' zero (there is no \code{initial} and no \code{bounds}). Fixed parameters (\code{arma}, numeric
#' \code{constant}) can only be provided for all of the parameters of the model at once.
#'
#' @template author
#' @template keywords
#'
#' @param y vector or ts object with the complex-valued time series.
#' @param orders orders of the model: either a vector \code{c(p,d,q)} or a list
#' \code{list(ar=p, i=d, ma=q)} (elements that are not provided are set to zero).
#' @param constant if \code{TRUE}, the constant term is included in the model and estimated.
#' Can also be a complex number (the value of the constant), in which case all the other
#' parameters need to be provided via \code{arma}.
#' @param arma list \code{list(ar=..., ma=...)} or a vector (AR parameters followed by MA ones)
#' with the fixed complex AR/MA parameters. All of them need to be provided; the model is then
#' not estimated.
#' @param model previously estimated \code{carima} model. If provided, its orders, constant,
#' loss, distribution and parameters are used, and the model is applied to \code{y} without
#' re-estimation.
#' @param loss the loss function used in the estimation, see \link[complex]{clm}.
#' @param distribution the distribution used in the likelihood, see \link[complex]{clm}.
#' @param h the forecast horizon. If positive, the forecasts are produced and returned in
#' \code{forecast}.
#' @param holdout if \code{TRUE}, the last \code{h} observations are held out: the model is
#' estimated on the rest and its forecasts are compared with the held-out values.
#' @param silent if \code{FALSE}, the function plots the fit and the forecast at the end.
#' @param ... other parameters passed to \link[complex]{clm}, e.g. \code{shape} for
#' \code{distribution="dcgnorm"}, \code{B} or the parameters of the optimiser.
#'
#' @return An object of the class \code{c("carima","clm","greybox")}: the model returned by
#' \link[complex]{clm}, plus
#' \itemize{
#' \item \code{orders} - the list of orders,
#' \item \code{constant} - whether the constant is included,
#' \item \code{h}, \code{holdout} - the horizon and the holdout flag,
#' \item \code{forecast} - the point forecasts (if \code{h>0}),
#' \item \code{holdoutData} - the held-out values (if \code{holdout=TRUE}),
#' \item \code{accuracy} - the error measures on the holdout: real and imaginary parts of the
#' mean error (ME_r, ME_i), MAE (mean modulus of the error) and RMSE (root mean squared modulus
#' of the error).
#' }
#'
#' @seealso \code{\link[complex]{clm}}, \code{\link[smooth]{msarima}}
#'
#' @examples
#' set.seed(41)
#' e <- rcnorm(200, 0, sigma2=1, varsigma2=0.3+0.2i)
#' y <- cumsum(e + c(0, 0.5*e[-200]))
#' model <- carima(y, orders=c(0,1,1), h=10, holdout=TRUE)
#' summary(model)
#' model$accuracy
#' forecast(model, h=5)
#'
#' @rdname carima
#' @export
carima <- function(y, orders=list(ar=0, i=0, ma=0), constant=FALSE, arma=NULL, model=NULL,
                   loss=c("likelihood","OLS","CLS","MSE","MAE","HAM"),
                   distribution=c("dcnorm","dcgnorm"),
                   h=0, holdout=FALSE, silent=TRUE, ...){
    cl <- match.call();
    ellipsis <- list(...);

    # Reuse the previously estimated model
    if(!is.null(model)){
        if(!inherits(model, "carima")){
            stop("model should be estimated via carima().", call.=FALSE);
        }
        orders <- model$orders;
        constant <- model$constant;
        loss <- model$loss;
        distribution <- model$distribution;
        parameters <- coef(model);
        if(distribution=="dcgnorm"){
            ellipsis$shape <- model$other$shape;
        }
    }
    else{
        loss <- match.arg(loss);
        distribution <- match.arg(distribution);
        parameters <- NULL;
    }

    orders <- carimaOrders(orders);
    ordersVector <- c(orders$ar, orders$i, orders$ma);

    # Constant: FALSE, TRUE (estimate) or a value
    if(is.logical(constant)){
        constantValue <- NULL;
    }
    else{
        constantValue <- as.complex(constant);
        constant <- TRUE;
    }

    # Fixed parameters: all or nothing
    if(is.null(model) && (!is.null(arma) || !is.null(constantValue))){
        if(is.list(arma)){
            arma <- c(arma$ar, arma$ma);
        }
        arma <- as.complex(arma);
        if(length(arma)!=orders$ar+orders$ma){
            stop("arma should contain all ", orders$ar+orders$ma, " AR/MA parameters.", call.=FALSE);
        }
        if(constant && is.null(constantValue)){
            stop("The constant needs to be provided as a value when arma is provided, ",
                 "because clm() cannot estimate a part of the parameters.", call.=FALSE);
        }
        parameters <- c(constantValue, arma);
    }

    # Holdout
    y <- as.vector(y);
    if(!is.complex(y)){
        y <- as.complex(y);
    }
    obsAll <- length(y);
    if(holdout){
        if(h<=0){
            stop("holdout=TRUE requires h>0.", call.=FALSE);
        }
        holdoutData <- y[obsAll-h+seq_len(h)];
        y <- y[seq_len(obsAll-h)];
    }

    # Estimate via clm
    formula <- if(constant) y~1 else y~-1;
    clmArgs <- c(list(formula=formula, data=data.frame(y=y), loss=loss, distribution=distribution,
                      orders=ordersVector, parameters=parameters),
                 ellipsis);
    # clm warns about the missing intercept, which is intended in carima()
    ourModel <- suppressWarnings(do.call(clm, clmArgs));
    if(!is.null(parameters)){
        names(ourModel$coefficients) <- names(ourModel$B) <- carimaNames(orders, constant);
    }

    # The clm call is kept in ourModel$call, because vcov.clm() re-estimates the model from it
    ourModel$carimaCall <- cl;
    ourModel$orders <- orders;
    ourModel$constant <- constant;
    ourModel$h <- h;
    ourModel$holdout <- holdout;
    ourModel$model <- paste0("cARIMA(", paste0(ordersVector, collapse=","), ")",
                             ifelse(constant, " with constant", ""));
    class(ourModel) <- c("carima", class(ourModel));

    if(h>0){
        ourModel$forecast <- as.vector(predict(ourModel, newdata=matrix(NA, h, 1))$mean);
        if(holdout){
            ourModel$holdoutData <- holdoutData;
            errors <- holdoutData - ourModel$forecast;
            ourModel$accuracy <- c(ME_r=Re(mean(errors)), ME_i=Im(mean(errors)),
                                   MAE=mean(Mod(errors)), RMSE=sqrt(mean(Mod(errors)^2)));
        }
    }

    if(!silent){
        if(h>0){
            plot(forecast(ourModel, h=h));
        }
        else{
            plot(ourModel, 7);
        }
    }

    return(ourModel);
}

# Orders as a list(ar, i, ma) with single non-negative values
carimaOrders <- function(orders){
    if(is.list(orders)){
        orders <- list(ar=if(is.null(orders$ar)) 0 else orders$ar,
                       i=if(is.null(orders$i)) 0 else orders$i,
                       ma=if(is.null(orders$ma)) 0 else orders$ma);
    }
    else{
        if(length(orders)!=3){
            stop("orders should be a vector c(p,d,q) or a list(ar=p, i=d, ma=q).", call.=FALSE);
        }
        orders <- list(ar=orders[1], i=orders[2], ma=orders[3]);
    }
    if(any(sapply(orders, length)!=1)){
        stop("Only non-seasonal cARIMA is supported: each order should be a single number.",
             call.=FALSE);
    }
    if(any(unlist(orders)<0)){
        stop("Orders should be non-negative.", call.=FALSE);
    }
    return(orders);
}

# Parameter names as clm produces them
carimaNames <- function(orders, constant){
    c(if(constant) "(Intercept)",
      if(orders$ar>0) paste0("yLag", seq_len(orders$ar)),
      if(orders$ma>0) paste0("eLag", seq_len(orders$ma)));
}

#' @param object the model estimated via \code{carima()}.
#' @param level the confidence level for the prediction interval.
#' @param interval the type of interval, see \link[complex]{clm}.
#' @rdname carima
#' @importFrom greybox forecast
#' @importFrom stats predict
#' @export
forecast.carima <- function(object, h=10, interval=c("none", "confidence", "prediction"),
                            level=0.95, ...){
    interval <- match.arg(interval);
    return(predict(object, newdata=matrix(NA, h, 1), interval=interval, level=level, ...));
}

#' @export
print.carima <- function(x, ...){
    cat(x$model, "estimated via carima()\n");
    cat("Loss function used in estimation:", x$loss, "\n");
    if(x$distribution=="dcgnorm"){
        cat("Distribution: complex generalised normal, shape =", round(x$other$shape, 4), "\n");
    }
    cat("Coefficients:\n");
    print(coef(x));
    if(!is.na(logLik(x))){
        ICs <- c(AIC=AIC(x), AICc=AICc(x), BIC=BIC(x), BICc=BICc(x));
        cat("Information criteria:\n");
        print(round(ICs, 4));
    }
    if(!is.null(x$accuracy)){
        cat("Holdout accuracy (h=", x$h, "):\n", sep="");
        print(x$accuracy);
    }
}

# vcov.clm() re-estimates the model via clm(), which warns about the missing intercept.
# In carima() this is intended, so the warning is suppressed.
#' @export
vcov.carima <- function(object, ...){
    return(suppressWarnings(NextMethod()));
}
