# Information criteria of clm use the total number of real parameters
# (both parts of the complex variable), see Bedrick & Tsai (1994).
set.seed(41)
obs <- 110
e <- rcnorm(obs, 0, sigma2=25, varsigma2=16+9i)
y <- vector("complex", obs)
y[1] <- 100+50i + e[1]
for(i in 2:obs){
    y[i] <- (0.5+0.2i)*y[i-1] + 100+50i + e[i]
}
y <- y[-c(1:10)]
cAR1Model <- clm(y~1, orders=c(1,0,0), subset=c(1:80))

test_that("nparam counts per series, logLik df counts all real parameters", {
    # Intercept + AR(1): 2 complex = 4 real, covariance matrix: 3 real
    expect_equal(nparam(cAR1Model), 3.5)
    expect_equal(attr(logLik(cAR1Model), "df"), 7)
    expect_equal(cAR1Model$df.residual, 80 - 3.5)
})

test_that("AIC and BIC use the total number of real parameters", {
    ll <- as.numeric(logLik(cAR1Model))
    expect_equal(AIC(cAR1Model), -2*ll + 2*7)
    expect_equal(BIC(cAR1Model), -2*ll + log(80)*7)
})

test_that("AICc and BICc use the multivariate correction with two series", {
    ll <- as.numeric(logLik(cAR1Model))
    # 2 complex coefficients per series, 2 series
    denominator <- 80 - (2 + 2 + 1)
    expect_equal(AICc(cAR1Model), -2*ll + 2*80*7/denominator)
    expect_equal(BICc(cAR1Model), -2*ll + log(80)*80*7/denominator)
})

test_that("ICs are not produced for non-likelihood losses", {
    cAR1OLS <- clm(y~1, orders=c(1,0,0), subset=c(1:80), loss="OLS")
    expect_true(is.na(logLik(cAR1OLS)))
})
