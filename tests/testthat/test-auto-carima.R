# auto.carima(): order selection for complex ARIMA

test_that("auto.carima() finds the same model as the full search, and cARIMA(1,1,1) without constant", {
    set.seed(41)
    obs <- 400
    e <- rcnorm(obs, 0, sigma2=1, varsigma2=0.3+0.2i)
    dy <- e
    for(t in 2:obs){
        dy[t] <- (0.5+0.2i)*dy[t-1] + e[t] + (0.3-0.1i)*e[t-1]
    }
    y <- cumsum(dy)
    # With all four constants, cARIMA(2,0,1) with a real constant has a slightly lower AICc than
    # the true model on this sample, so the search is checked against the full search, and the
    # orders against the truth with the constant fixed
    model <- auto.carima(y, orders=list(ar=2, i=1, ma=2))
    expect_s3_class(model, "carima")
    modelFull <- auto.carima(y, orders=list(ar=2, i=1, ma=2), search="full")
    expect_equal(modelFull$orders, model$orders)
    expect_equal(modelFull$constant, model$constant)
    expect_equal(nrow(modelFull$ICs), 3*3*4*2)
    model <- auto.carima(y, orders=list(ar=2, i=1, ma=2), constant="none")
    expect_equal(model$orders, list(ar=1, i=1, ma=1))
    expect_equal(model$constant, "none")
})

test_that("auto.carima() finds the random walk with drift and the white noise", {
    set.seed(7)
    y <- cumsum(rcnorm(400, 0.3+0.1i, sigma2=1, varsigma2=0.2))
    model <- auto.carima(y, orders=list(ar=2, i=1, ma=2))
    expect_equal(model$orders, list(ar=0, i=1, ma=0))
    expect_equal(model$constant, "both")
    y <- 5+2i + rcnorm(400, 0, sigma2=1, varsigma2=0.2)
    model <- auto.carima(y, orders=list(ar=2, i=1, ma=2))
    expect_equal(model$orders, list(ar=0, i=0, ma=0))
    expect_equal(model$constant, "both")
})

test_that("the information criteria are calculated on the common sample", {
    set.seed(41)
    y <- cumsum(rcnorm(300, 0, sigma2=1, varsigma2=0.2))
    model <- auto.carima(y, orders=list(ar=1, i=1, ma=1), h=10, holdout=TRUE)
    expect_equal(model$icObsDropped, 3)
    expect_equal(nobs(model), 290)
    expect_length(model$forecast, 10)
    expect_true(!is.null(model$accuracy))
    # Without dropping observations the criterion would be the one of the model
    modelFit <- carima(y[1:290], orders=c(model$orders$ar, model$orders$i, model$orders$ma),
                       constant=model$constant)
    # The selected model has an AR coefficient close to one, where the two fits from different
    # starting values agree only to about 0.2%
    expect_equal(coef(modelFit), coef(model), tolerance=1e-2)
})

test_that("clm() penalises non-stationary and non-invertible parameters", {
    set.seed(41)
    y <- cumsum(rcnorm(300, 0, sigma2=1, varsigma2=0.2))
    model <- carima(y, orders=c(1,1,0), arma=list(ar=1.2))
    expect_lt(as.numeric(logLik(model)), -1e90)
    model <- carima(y, orders=c(0,1,1), arma=list(ma=0.5-1i))
    expect_lt(as.numeric(logLik(model)), -1e90)
    model <- carima(y, orders=c(1,1,1))
    expect_true(all(Mod(polyroot(c(1, -coef(model)[1])))>1))
    expect_true(all(Mod(polyroot(c(1, coef(model)[2])))>1))
})

test_that("auto.carima(fast=TRUE) finds the same models via Hannan-Rissanen screening", {
    set.seed(41)
    obs <- 400
    e <- rcnorm(obs, 0, sigma2=1, varsigma2=0.3+0.2i)
    dy <- e
    for(t in 2:obs){
        dy[t] <- (0.5+0.2i)*dy[t-1] + e[t] + (0.3-0.1i)*e[t-1]
    }
    y <- cumsum(dy)
    model <- auto.carima(y, orders=list(ar=2, i=1, ma=2), fast=TRUE)
    modelStepwise <- auto.carima(y, orders=list(ar=2, i=1, ma=2))
    expect_equal(model$orders, modelStepwise$orders)
    expect_equal(model$constant, modelStepwise$constant)
    expect_true(!is.null(model$ICsScreening))
    expect_equal(nrow(model$ICsScreening), 3*3*4*2)
    set.seed(7)
    y <- cumsum(rcnorm(400, 0.3+0.1i, sigma2=1, varsigma2=0.2))
    model <- auto.carima(y, orders=list(ar=2, i=1, ma=2), fast=TRUE)
    expect_equal(model$orders, list(ar=0, i=1, ma=0))
    expect_equal(model$constant, "both")
    expect_warning(auto.carima(y, orders=list(ar=1, i=1, ma=1), search="full", fast=TRUE), "fast")
})
