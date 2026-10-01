# Partial intercept: y~x-Im(1) (real part only), y~x-Re(1) (imaginary part only)
set.seed(41)
obs <- 400
e <- rcnorm(obs, 0, sigma2=1, varsigma2=0.3+0.2i)
# Drift in the real part only
yReal <- cumsum(0.2 + e + c(0, (0.4-0.1i)*e[-obs]))
partialData <- data.frame(y=yReal)

test_that("the formula defines the type of the intercept", {
    expect_equal(complex:::clmInterceptType(y~x, TRUE), "both")
    expect_equal(complex:::clmInterceptType(y~x-Im(1), TRUE), "real")
    expect_equal(complex:::clmInterceptType(y~x+z-Re(1), TRUE), "imaginary")
    expect_equal(complex:::clmInterceptType(y~x-1, FALSE), "none")
    expect_error(complex:::clmInterceptType(y~x-Im(1)-Re(1), TRUE), "Both parts")
    expect_equal(complex:::clmInterceptStrip(y~x-Im(1)), y~x)
})

test_that("clm() estimates the real part of the intercept only", {
    model <- clm(y~1-Im(1), partialData, orders=c(0,1,1))
    modelFull <- clm(y~1, partialData, orders=c(0,1,1))
    expect_equal(unname(Im(coef(model)[1])), 0)
    expect_equal(model$other$intercept, "real")
    expect_equal(nparam(model), nparam(modelFull)-0.5)
    # The true imaginary drift is zero, so the restriction costs little likelihood
    expect_lt(as.numeric(logLik(modelFull) - logLik(model)), 2)
    expect_equal(deparse(model$call$formula), "y ~ 1 - Im(1)")
    # vcov: zero variance for the fixed part, finite for the others
    covariance <- vcov(model)
    expect_equal(unname(covariance["(Intercept)_i",]), rep(0, 4))
    expect_true(all(diag(covariance)[-2]>0))
    expect_false("(Intercept)_i" %in% rownames(coef(summary(model))))
    expect_true(all(is.finite(predict(model, newdata=matrix(NA, 5, 1))$mean)))
})

test_that("clm() estimates the imaginary part of the intercept only", {
    model <- clm(y~1-Re(1), partialData, orders=c(1,1,1), distribution="dcgnorm")
    expect_equal(unname(Re(coef(model)[1])), 0)
    expect_equal(model$other$intercept, "imaginary")
    expect_false("(Intercept)_r" %in% rownames(coef(summary(model))))
    expect_error(clm(y~1-Im(1), partialData, loss="CLS"), "CLS")
})

test_that("carima() supports the four constants", {
    models <- lapply(c("none","real","imaginary","both"), function(constant){
        carima(yReal, orders=c(0,1,1), constant=constant)
    })
    expect_equal(sapply(models, function(x) x$constant), c("none","real","imaginary","both"))
    expect_equal(sapply(models, nparam), c(2.5, 3, 3, 3.5))
    modelClm <- clm(y~1-Im(1), partialData, orders=c(0,1,1))
    expect_equal(unname(coef(models[[2]])), unname(coef(modelClm)))
    expect_match(models[[2]]$model, "real constant")
    # Reuse keeps the restriction
    modelReused <- carima(yReal[1:300], model=models[[2]])
    expect_equal(coef(modelReused), coef(models[[2]]))
    expect_equal(modelReused$other$intercept, "real")
})

test_that("auto.carima() selects the real constant for a drift in the real part", {
    model <- auto.carima(yReal, orders=list(ar=1, i=1, ma=1))
    expect_equal(model$constant, "real")
    expect_equal(model$orders$i, 1)
    expect_true(all(c("none","real","imaginary","both") %in% model$ICs$constant))
    modelFixed <- auto.carima(yReal, orders=list(ar=1, i=1, ma=1), constant="imaginary")
    expect_true(all(modelFixed$ICs$constant=="imaginary"))
})
