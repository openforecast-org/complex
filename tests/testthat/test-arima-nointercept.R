# cARIMA without intercept (formula y~-1): coefficients, residuals and forecasts
set.seed(3)
obs <- 3000
e <- rcnorm(obs, 0, sigma2=1, varsigma2=0.3)
dy <- e
for(t in 2:obs){
    dy[t] <- (0.5+0.2i)*dy[t-1] + e[t] + (0.3-0.1i)*e[t-1]
}
y <- cumsum(dy)
data <- data.frame(y=y)

test_that("cARIMA(1,1,1) without intercept estimates the parameters with correct signs", {
    model <- suppressWarnings(clm(y~-1, data, orders=c(1,1,1)))
    expect_equal(unname(coef(model)), c(0.5+0.2i, 0.3-0.1i), tolerance=0.05)
    # Residuals follow the model recursion
    b <- coef(model)
    residualsManual <- complex(obs)
    for(t in 3:obs){
        residualsManual[t] <- (y[t]-y[t-1]) - b[1]*(y[t-1]-y[t-2]) - b[2]*residualsManual[t-1]
    }
    expect_equal(unname(residuals(model)[21:obs]), residualsManual[21:obs], tolerance=1e-6)
    # Forecasts follow the model recursion
    yExtended <- y
    residualsExtended <- residuals(model)
    for(i in 1:5){
        n <- length(yExtended)
        yExtended <- c(yExtended, yExtended[n] + b[1]*(yExtended[n]-yExtended[n-1]) +
                           b[2]*residualsExtended[n])
        residualsExtended <- c(residualsExtended, 0)
    }
    expect_equal(as.vector(predict(model, newdata=matrix(NA, 5, 1))$mean), unname(tail(yExtended, 5)))
})

test_that("cARIMA(0,1,2) without intercept is not affected by the sign of MA", {
    model <- suppressWarnings(clm(y~-1, data, orders=c(0,1,2)))
    modelIntercept <- clm(y~1, data, orders=c(0,1,2))
    expect_equal(unname(coef(model)), unname(coef(modelIntercept)[-1]), tolerance=0.02)
})

test_that("cARIMA(0,1,0) without intercept has nothing to estimate and forecasts the last value", {
    for(distribution in c("dcnorm", "dcgnorm")){
        model <- suppressWarnings(clm(y~-1, data, orders=c(0,1,0), distribution=distribution))
        expect_true(is.finite(logLik(model)))
        expect_equal(as.vector(predict(model, newdata=matrix(NA, 3, 1))$mean), rep(y[obs], 3))
    }
})

test_that("cARIMA(0,0,q) without intercept can be estimated", {
    set.seed(4)
    e <- rcnorm(500, 0, sigma2=1, varsigma2=0.3)
    y <- e + c(0, (0.5-0.2i)*e[-500])
    for(loss in c("likelihood", "OLS")){
        model <- suppressWarnings(clm(y~-1, data.frame(y=y), orders=c(0,0,1), loss=loss))
        expect_equal(unname(coef(model)), 0.5-0.2i, tolerance=0.1)
    }
    # CLS is complex-valued and cannot be used in the numeric optimisation needed for MA
    expect_error(clm(y~-1, data.frame(y=y), orders=c(0,0,1), loss="CLS"), "cannot be used with MA")
})
