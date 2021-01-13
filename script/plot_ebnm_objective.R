library(microbenchmark)
library(ebnm)
library(stats)
library(rgl)
set.seed(123)

############## Functions for simulation and compute loglikelihood #############################
simulate_pn <- function(n, pn, null_case, homosked, est_mu) {
  pi0 <- rbeta(1, shape1 = 10, shape2 = 2)
  scale <- rgamma(1, shape = 4, rate = 1)

  if (est_mu) {
    mu <- runif(1, -20, 20)
  } else {
    mu <- 0
  }

  if (null_case) {
    theta <- mu
  } else {
    if (pn) {
      theta <- mu + rnorm(n, sd = scale)
    } else {
      theta <- mu + rexp(n, 1 / scale) * sample(c(-1, 1), n, replace = TRUE)
    }
    theta[rbinom(n, size = 1, prob = pi0) == 1] <- mu
  }

  if (homosked) {
    s <- 1
  } else {
    s <- sqrt(rexp(n))
  }

  x <- theta + rnorm(n, sd = s)
  out <- list(x = x, s = s, theta = theta)
  return(out)
}

# Functions to compute the log likelihood under the point-normal prior.
loglik_point_normal = function(x, s, w, a, mu) {
  return(sum(vloglik_point_normal(x, s, w, a, mu)))
}

vloglik_point_normal = function(x, s, w, a, mu) {
  if (w <= 0) {
    return(dnorm(x, mu, s, log = TRUE))
  }

  lg <- dnorm(x, mu, sqrt(s^2 + 1/a), log = TRUE)
  if (w >= 1) {
    return(lg)
  }

  lf <- dnorm(x, mu, s, log = TRUE)
  lfac <- pmax(lg, lf)
  result <- lfac + log((1 - w) * exp(lf - lfac) + w * exp(lg - lfac))

  if (any(s == 0)) {
    result[s == 0 & x == mu] <- log(1 - w)
    result[s == 0 & x != mu] <- log(w) + lg[s == 0 & x != mu]
  }

  return(result)
}

## simulation and plot loglikelihood terrain

data = simulate_pn(n = 1000, pn = TRUE, null_case = TRUE, homosked = TRUE, est_mu = TRUE)
I = 100
J = 100

pi0 = seq(0,1,length.out = I)
v = seq(2,5,length.out = J)^2
z = matrix(NA, nrow = I, ncol = J)
for(i in 1:I){
  for(j in 1:J){
    z[i,j] <- loglik_point_normal(x = data$x, s = data$s, w = pi0[i], a = 1/v[j], mu = 0)
  }
}

op <- par(bg = "white")
persp3d(x = pi0, y = v, z = z, theta = 30, phi = 30, expand = 0.5, col = "lightblue")

