# Complex generalised normal distribution and its use in clm()

test_that("dcgnorm with shape=1 equals dcnorm", {
    set.seed(41)
    z <- rcnorm(20, 1+1i, sigma2=2, varsigma2=0.6+0.8i)
    expect_equal(dcgnorm(z, 1+1i, scale=2, pseudoscale=0.6+0.8i, shape=1, log=TRUE),
                 dcnorm(z, 1+1i, sigma2=2, varsigma2=0.6+0.8i, log=TRUE),
                 tolerance=1e-12)
})

test_that("dcgnorm integrates to one", {
    grid <- seq(-40, 40, length.out=801)
    step <- diff(grid)[1]
    z <- outer(grid, 1i*grid, "+")
    for(shape in c(0.5, 1, 2)){
        expect_equal(sum(dcgnorm(z, 0, scale=2, pseudoscale=0.6+0.8i, shape=shape))*step^2,
                     1, tolerance=1e-3)
    }
})

test_that("rcgnorm has the right covariance and radial distribution", {
    set.seed(41)
    shape <- 0.5
    scale <- 2
    pseudoscale <- 0.6+0.8i
    z <- rcgnorm(50000, 0, scale, pseudoscale, shape)
    # Variance and pseudo-variance are the scale and pseudo-scale times this factor
    multiplier <- 2^(1/shape)*gamma(2/shape)/(2*gamma(1/shape))
    expect_equal(mean(Mod(z)^2), multiplier*scale, tolerance=0.05)
    expect_equal(mean(z^2), multiplier*pseudoscale, tolerance=0.05)
    # Q^shape follows Gamma(1/shape, scale=2)
    quadraticForm <- 2*(scale*Mod(z)^2 - Re(Conj(pseudoscale)*z^2))/(scale^2-Mod(pseudoscale)^2)
    expect_gt(suppressWarnings(ks.test(quadraticForm^shape, "pgamma", shape=1/shape, scale=2)$p.value),
              0.001)
})

test_that("the concentrated likelihood matches the density at the estimated scales", {
    set.seed(41)
    e <- rcgnorm(2000, 0, 2, 0.6+0.8i, 0.5)
    values <- complex:::cgnormConcentrated(e, 0.5)
    expect_equal(values$logLik,
                 sum(dcgnorm(e, 0, values$scale, values$pseudoscale, 0.5, log=TRUE)))
})

# cAR(1) with heavy-tailed errors
set.seed(41)
obs <- 2100
e <- rcgnorm(obs, 0, scale=4, pseudoscale=1+1i, shape=0.5)
y <- vector("complex", obs)
y[1] <- 10+5i + e[1]
for(i in 2:obs){
    y[i] <- (0.5+0.2i)*y[i-1] + 10+5i + e[i]
}
y <- y[-c(1:100)]
data <- data.frame(y=y)

test_that("clm with dcgnorm and shape=1 reproduces dcnorm", {
    modelNormal <- clm(y~1, data, orders=c(1,0,0))
    modelShapeOne <- clm(y~1, data, orders=c(1,0,0), distribution="dcgnorm", shape=1)
    expect_equal(coef(modelShapeOne), coef(modelNormal), tolerance=1e-4)
    expect_equal(as.numeric(logLik(modelShapeOne)), as.numeric(logLik(modelNormal)),
                 tolerance=1e-6)
    expect_equal(nparam(modelShapeOne), nparam(modelNormal))
})

test_that("clm with dcgnorm estimates the shape and beats dcnorm on heavy tails", {
    modelNormal <- clm(y~1, data, orders=c(1,0,0))
    modelGeneral <- clm(y~1, data, orders=c(1,0,0), distribution="dcgnorm")
    expect_true(modelGeneral$other$shapeEstimated)
    expect_gt(modelGeneral$other$shape, 0.4)
    expect_lt(modelGeneral$other$shape, 0.6)
    expect_equal(nparam(modelGeneral), nparam(modelNormal) + 0.5)
    expect_lt(AICc(modelGeneral), AICc(modelNormal))
    expect_equal(unname(coef(modelGeneral)[2]), 0.5+0.2i, tolerance=0.05)
    # summary and forecasts work
    expect_output(print(summary(modelGeneral)), "shape")
    expect_length(predict(modelGeneral, newdata=matrix(NA, 5, 1))$mean, 5)
})

test_that("dcgnorm is refused with non-likelihood losses", {
    expect_error(clm(y~1, data, loss="OLS", distribution="dcgnorm"))
})
