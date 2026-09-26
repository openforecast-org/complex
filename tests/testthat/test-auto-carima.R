# auto.carima(): order selection for complex ARIMA

test_that("auto.carima() finds cARIMA(1,1,1), the same as the full search", {
    set.seed(41)
    obs <- 400
    e <- rcnorm(obs, 0, sigma2=1, varsigma2=0.3+0.2i)
    dy <- e
    for(t in 2:obs){
        dy[t] <- (0.5+0.2i)*dy[t-1] + e[t] + (0.3-0.1i)*e[t-1]
    }
    y <- cumsum(dy)
    model <- auto.carima(y, orders=list(ar=2, i=1, ma=2))
    expect_s3_class(model, "carima")
    expect_equal(model$orders, list(ar=1, i=1, ma=1))
    expect_false(model$constant)
    modelFull <- auto.carima(y, orders=list(ar=2, i=1, ma=2), search="full")
    expect_equal(modelFull$orders, model$orders)
    expect_equal(nrow(modelFull$ICs), 3*3*2 + 3*3)
})

test_that("auto.carima() finds the random walk with drift and the white noise", {
    set.seed(7)
    y <- cumsum(rcnorm(400, 0.3+0.1i, sigma2=1, varsigma2=0.2))
    model <- auto.carima(y, orders=list(ar=2, i=1, ma=2))
    expect_equal(model$orders, list(ar=0, i=1, ma=0))
    expect_true(model$constant)
    y <- 5+2i + rcnorm(400, 0, sigma2=1, varsigma2=0.2)
    model <- auto.carima(y, orders=list(ar=2, i=1, ma=2))
    expect_equal(model$orders, list(ar=0, i=0, ma=0))
    expect_true(model$constant)
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
    expect_equal(coef(modelFit), coef(model), tolerance=1e-4)
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
