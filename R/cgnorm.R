#' Complex Generalised Normal Distribution
#'
#' Density and random generation for the complex generalised normal distribution
#' (also known as the complex generalised Gaussian or the bivariate exponential power
#' distribution) with location \code{mu}, scale \code{scale}, pseudo-scale
#' \code{pseudoscale} and shape \code{shape}.
#'
#' For the error term \eqn{e = z - \mu}, the density is
#' \deqn{f(z) = \frac{2\beta}{\pi \Gamma(1/\beta) 2^{1/\beta} \sqrt{\alpha^4 - |\psi^2|^2}}
#' \exp\left(-\frac{1}{2} Q^\beta\right),}
#' where \eqn{Q = 2(\alpha^2 |e|^2 - \Re(\overline{\psi^2} e^2)) / (\alpha^4 - |\psi^2|^2)},
#' \eqn{\alpha^2} is the scale, \eqn{\psi^2} is the pseudo-scale and \eqn{\beta} is the shape.
#'
#' With \code{shape=1} this is the complex normal distribution with variance \code{scale} and
#' pseudo-variance \code{pseudoscale} (see \link[complex]{dcnorm}). Values of shape below one
#' give heavier tails (\code{shape=0.5} is Laplace-like), values above one give lighter tails
#' (the uniform distribution on an ellipse in the limit). Note that the shape here applies to
#' the quadratic form, so \code{shape=1} corresponds to \code{shape=2} in
#' \link[greybox]{dgnorm}.
#'
#' The variance and the pseudo-variance of the distribution equal \code{scale} and
#' \code{pseudoscale} multiplied by \eqn{2^{1/\beta} \Gamma(2/\beta) / (2 \Gamma(1/\beta))}.
#'
#' Random numbers are generated using the fact that \eqn{Q^\beta} follows the Gamma distribution
#' with shape \eqn{1/\beta} and scale 2, and the angle is uniform in the standardised space.
#'
#' @template author
#' @keywords distribution
#'
#' @param q vector of quantiles (complex).
#' @param n number of observations. Should be a single number.
#' @param mu vector of location parameters (complex).
#' @param scale scale parameter (real, positive).
#' @param pseudoscale pseudo-scale parameter (complex, with the modulus smaller than scale).
#' @param shape shape parameter (real, positive).
#' @param log if \code{TRUE}, then the logarithm of the density is returned.
#'
#' @return \code{dcgnorm} returns the density, \code{rcgnorm} returns a vector of
#' complex random values.
#'
#' @references \itemize{
#' \item Novey, M., Adali, T., Roy, A. (2010). A complex generalized Gaussian distribution -
#' characterization, generation, and estimation. IEEE Transactions on Signal Processing,
#' 58(3), 1427-1433.
#' \item Gomez, E., Gomez-Villegas, M. A., Marin, J. M. (1998). A multivariate generalization
#' of the power exponential family of distributions. Communications in Statistics - Theory
#' and Methods, 27(3), 589-600.
#' }
#'
#' @seealso \code{\link[complex]{dcnorm}}
#'
#' @examples
#' x <- rcgnorm(1000, 1+1i, scale=2, pseudoscale=0.6+0.8i, shape=0.5)
#' plot(x)
#' dcgnorm(x[1:5], 1+1i, scale=2, pseudoscale=0.6+0.8i, shape=0.5)
#'
#' @rdname cgnormal
#' @aliases cgnormal dcgnorm
#' @export dcgnorm
dcgnorm <- function(q, mu=0, scale=1, pseudoscale=0, shape=1, log=FALSE){
    cgnormCheck(scale, pseudoscale, shape);
    errors <- as.vector(q - mu);
    determinant <- scale^2 - Mod(pseudoscale)^2;
    quadraticForm <- 2*(scale*Mod(errors)^2 - Re(Conj(pseudoscale)*errors^2)) / determinant;
    cgnormReturn <- (log(2*shape) - log(pi) - lgamma(1/shape) - log(2)/shape -
                         0.5*log(determinant) - 0.5*quadraticForm^shape);
    if(!log){
        cgnormReturn <- exp(cgnormReturn);
    }
    return(cgnormReturn);
}

#' @rdname cgnormal
#' @export rcgnorm
#' @aliases rcgnorm
#' @importFrom stats rgamma runif
rcgnorm <- function(n=1, mu=0, scale=1, pseudoscale=0, shape=1){
    cgnormCheck(scale, pseudoscale, shape);
    # Covariance-type matrix of the real and imaginary parts (as in dcnorm)
    Sigma <- matrix(c(scale+Re(pseudoscale), Im(pseudoscale),
                      Im(pseudoscale), scale-Re(pseudoscale))/2, 2, 2);
    radius <- sqrt(rgamma(n, shape=1/shape, scale=2)^(1/shape));
    angle <- runif(n, 0, 2*pi);
    cgnormReturn <- t(chol(Sigma)) %*% rbind(radius*cos(angle), radius*sin(angle));
    return(mu + complex(real=cgnormReturn[1,], imaginary=cgnormReturn[2,]));
}

cgnormCheck <- function(scale, pseudoscale, shape){
    if(any(scale<=0)){
        stop("The scale should be positive.", call.=FALSE);
    }
    if(any(Mod(pseudoscale)>=scale)){
        stop("The modulus of the pseudoscale should be smaller than the scale.", call.=FALSE);
    }
    if(any(shape<=0)){
        stop("The shape should be positive.", call.=FALSE);
    }
}

# Concentrated likelihood of the complex generalised normal for given errors and shape.
# The circularity coefficient rho = pseudoscale / scale is estimated from moments
# (exact ML only when shape=1); the scale is the exact ML estimate given rho.
cgnormConcentrated <- function(errors, shape){
    obs <- length(errors);
    rho <- sum(errors^2) / sum(Mod(errors)^2);
    quadraticForm <- 2*(Mod(errors)^2 - Re(Conj(rho)*errors^2)) / (1-Mod(rho)^2);
    scale <- (shape/(2*obs) * sum(quadraticForm^shape))^(1/shape);
    logLik <- obs*(log(2*shape) - log(pi) - lgamma(1/shape) - log(2)/shape -
                       log(scale) - 0.5*log(1-Mod(rho)^2) - 1/shape);
    return(list(scale=scale, pseudoscale=rho*scale, logLik=logLik));
}
