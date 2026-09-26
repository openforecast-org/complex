suppressMessages(library(complex))
# complex GG density, scatter (alpha2 real, psi2 complex)
dcgg <- function(z, alpha2, psi2, beta, log=TRUE){
  D <- alpha2^2 - Mod(psi2)^2
  Q <- 2*(alpha2*Mod(z)^2 - Re(Conj(psi2)*z^2))/D
  l <- log(2*beta) - log(pi) - lgamma(1/beta) - log(2)/beta - 0.5*log(D) - 0.5*Q^beta
  if(log) l else exp(l)
}
rcgg <- function(n, alpha2, psi2, beta){
  S <- matrix(c(alpha2+Re(psi2), Im(psi2), Im(psi2), alpha2-Re(psi2)),2)/2
  L <- t(chol(S)); r <- sqrt(rgamma(n, 1/beta, scale=2)^(1/beta)); ph <- runif(n,0,2*pi)
  x <- L %*% rbind(r*cos(ph), r*sin(ph)); complex(real=x[1,], imaginary=x[2,])
}
# closed-form concentrated logLik: rho from moments, alpha2 exact given rho
concLL <- function(e, beta){
  T <- length(e); rho <- mean(e^2)/mean(Mod(e)^2)
  q <- 2*(Mod(e)^2 - Re(Conj(rho)*e^2))/(1-Mod(rho)^2)
  a2 <- ((beta/(2*T))*sum(q^beta))^(1/beta)
  ll <- T*(log(2*beta) - log(pi) - lgamma(1/beta) - log(2)/beta - log(a2) - 0.5*log(1-Mod(rho)^2) - 1/beta)
  list(ll=ll, alpha2=a2, psi2=rho*a2)
}
set.seed(1)
# 1. beta=1 equals dcnorm
z <- rcnorm(5, 0, sigma2=2, varsigma2=0.6+0.8i)
cat("beta=1 vs dcnorm max diff:", max(abs(dcgg(z,2,0.6+0.8i,1) - dcnorm(z,0,2,0.6+0.8i,log=TRUE))), "\n")
# 2. density integrates to 1 (beta=0.5)
g <- seq(-60,60,length.out=1201); zz <- outer(g, 1i*g, "+")
cat("integral beta=0.5:", sum(dcgg(zz,2,0.6+0.8i,0.5,log=FALSE))*diff(g)[1]^2, "\n")
# 3. simulate beta=0.5, compare closed-form vs full numeric MLE of (alpha2, psi2)
e <- rcgg(20000, 2, 0.6+0.8i, 0.5)
cf <- concLL(e, 0.5)
negLL <- function(p) -sum(dcgg(e, exp(p[1]), complex(real=p[2],imaginary=p[3]), 0.5))
op <- optim(c(log(2),0.5,0.5), negLL, control=list(reltol=1e-12, maxit=5000))
cat("closed-form: LL", round(cf$ll,2), " alpha2", round(cf$alpha2,4), " psi2", format(round(cf$psi2,4)), "\n")
cat("full MLE   : LL", round(-op$value,2), " alpha2", round(exp(op$par[1]),4), " psi2", format(round(complex(real=op$par[2],imaginary=op$par[3]),4)), "\n")
# 4. beta profile from closed form
b <- seq(0.3,1.5,0.01); cat("beta maximising closed-form LL:", b[which.max(sapply(b, function(x) concLL(e,x)$ll))], "(true 0.5)\n")
