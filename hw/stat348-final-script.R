rm(list = ls())
library(extraDistr)
library(geosphere)
library(mvtnorm)
library(clustermq)
library(future.apply)
plan(multiprocess) ## => parallelize on your local computer
set.seed(123)


## load data &  fitted results

load("stat348_finalp2_3.RData")
fitted_par = c(mydata$mu, mydata$s, fit.gp$par)

cond_exp <- function(x_other, mu, S, idx){
  mu[idx] + S[idx,-idx] %*% solve(S[-idx,-idx], x_other - mu[-idx])
}

cv_one <- function(x, mu, S, idx){
  cond_exp(x[-idx], mu, S, idx)
}

Kernel <- function(lam, l, d){
  (lam^2) * exp(-(d/l)^2/2)
}


get_mu_S <- function(par, y, D){
  n = length(y)
  mu = par[1]
  s = par[2]
  lam = par[3]
  l = par[4]
  mu = rep(mu, n)
  S = Kernel(lam, l, D)
  diag(S) = diag(S) + s^2
  return(list(S = S, mu = mu))
}

tmp = get_mu_S(fitted_par, y = mydata$y, D = mydata$D)



lapply(1:3,cv_one, x = mydata$y, S = tmp$S, mu = tmp$mu)
