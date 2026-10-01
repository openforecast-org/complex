set.seed(41)
callData <- data.frame(y=cumsum(rcnorm(300, 0, 1, 0.3+0.2i)))

test_that("vcov works when clm is called from a function with variables as arguments", {
    # vcov.clm() re-estimates the model from the stored call; the variables of the wrapper
    # (dist, sc, start) do not exist there, so clm must store their values
    wrapper <- function(dist, sc, start){
        clm(y~1, callData, orders=c(0,1,1), distribution=dist, scaling=sc, B=start);
    }
    fitWrapped <- wrapper("dcnorm", "normalisation", c(0+0i, 0.1+0i))
    fitDirect <- clm(y~1, callData, orders=c(0,1,1), distribution="dcnorm",
                     scaling="normalisation", B=c(0+0i, 0.1+0i))
    expect_equal(fitWrapped$call$distribution, "dcnorm")
    expect_equal(vcov(fitWrapped), vcov(fitDirect))
    expect_s3_class(summary(fitWrapped), "summary.clm")
})

test_that("vcov works for dcgnorm called from a function", {
    wrapper <- function(dist) clm(y~1, callData, orders=c(0,1,1), distribution=dist)
    fit <- wrapper("dcgnorm")
    expect_true(all(is.finite(diag(vcov(fit)))))
})
