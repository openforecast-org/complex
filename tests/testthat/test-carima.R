# carima(): wrapper of clm() for pure complex ARIMA
set.seed(41)
obs <- 300
e <- rcnorm(obs, 0, sigma2=1, varsigma2=0.3+0.2i)
y <- cumsum(e + c(0, (0.5-0.2i)*e[-obs]))

test_that("carima() gives the same model as clm() without intercept", {
    model <- carima(y, orders=c(0,1,1))
    modelClm <- suppressWarnings(clm(y~-1, data.frame(y=y), orders=c(0,1,1)))
    expect_s3_class(model, "carima")
    expect_equal(unname(coef(model)), unname(coef(modelClm)))
    expect_equal(as.numeric(logLik(model)), as.numeric(logLik(modelClm)))
    expect_equal(model$orders, list(ar=0, i=1, ma=1))
})

test_that("orders can be a list and the constant is included when asked", {
    model <- carima(y, orders=list(i=1, ma=1), constant=TRUE)
    modelClm <- clm(y~1, data.frame(y=y), orders=c(0,1,1))
    expect_equal(unname(coef(model)), unname(coef(modelClm)))
    expect_error(carima(y, orders=list(ar=c(1,1))), "non-seasonal")
})

test_that("holdout and forecast work", {
    model <- carima(y, orders=c(0,1,1), h=10, holdout=TRUE)
    expect_equal(nobs(model), obs - 10)
    expect_length(model$forecast, 10)
    expect_equal(model$holdoutData, y[obs-9:0])
    expect_equal(unname(model$accuracy["RMSE"]),
                 sqrt(mean(Mod(model$holdoutData - model$forecast)^2)))
    expect_equal(as.vector(forecast(model, h=10)$mean), model$forecast)
    expect_no_warning(summary(model))
})

test_that("model and arma reuse parameters without estimation", {
    model <- carima(y, orders=c(0,1,1))
    modelReused <- carima(y[1:250], model=model)
    expect_equal(coef(modelReused), coef(model))
    modelFixed <- carima(y, orders=c(0,1,1), arma=list(ma=0.5-0.2i))
    expect_equal(unname(coef(modelFixed)), 0.5-0.2i)
    expect_error(carima(y, orders=c(0,1,1), constant=TRUE, arma=list(ma=0.5)), "constant")
    expect_error(carima(y, orders=c(1,1,1), arma=list(ma=0.5)), "all")
})

test_that("carima() works with the complex generalised normal", {
    model <- carima(y, orders=c(0,1,1), distribution="dcgnorm")
    expect_true(model$other$shapeEstimated)
    modelReused <- carima(y[1:250], model=model)
    expect_equal(modelReused$other$shape, model$other$shape)
})
