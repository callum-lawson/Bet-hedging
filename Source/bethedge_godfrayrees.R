##############################################################################
# Old version of threshold germination simulations, kept for its annotations #
##############################################################################

# Define functions --------------------------------------------------------

logistic <- function(x,alpha_G,beta_Gz){
  plogis(alpha_G + beta_Gz*x)
  }

# Set up sims -------------------------------------------------------------

nt <- 100
np <- 100

zmu <- 0
zsd <- 1

iota_mu_min <- -3
iota_mu_max <- 3
iota_sd_min <- 0.01
iota_sd_max <- 2

iota_mu_seq <- seq(iota_mu_min,iota_mu_max,length.out=np) # based on real data
iota_sd_seq <- seq(iota_sd_min,iota_sd_max,length.out=np) # based on real data
iota_mu <- iota_sd <- matrix(NA,nr=np,nc=np) 
iota_mu[] <- rep(iota_mu_seq,times=np)
iota_sd[] <- rep(iota_sd_seq,each=np)
alpha_Y <- 1
beta_Yz <- 2
sig_Y <- 0 # 10
alpha_m <- 1

alpha_G <- beta_Gz <- rbar <- rcns <- matrix(NA,nr=np,nc=np) 
alpha_G <- -iota_mu*sqrt(pi^2/(3*iota_sd^2))
beta_Gz <- sqrt(pi^2/(3*iota_sd^2))
  # calculated from Godfray & Rees 2002

  # if beta_Gz=0 (flat, random prob), then mean and variance undefined
  # but setting beta_Gz to very small shows that results from strategy with
  # very large variance in thresholds
  # (so by one definition, no plasticity is max bet-hedging)
  # is this general, or it tied to underlying idea of thresholds, 
  # which isn't necessarily true?

set.seed(1)
zt <- rnorm(nt,mean=zmu,sd=zsd)
eps_Y <- rnorm(nt,0,sig_Y)

# Density-independent -----------------------------------------------------

Xt <- vector("numeric",length=nt) 
rbar <- matrix(NA,nr=np,nc=np) 

for(i in 1:np^2){
  Xt[1] <- 1
  S <- exp(-alpha_m)  # same for all i
  G <- plogis(alpha_G[i] + beta_Gz[i]*zt)
  Y <- exp(alpha_Y + beta_Yz*zt + eps_Y)  # same for all i
  Gcns <- plogis(alpha_G[i] + beta_Gz[i]*median(zt))
  Ycns <- exp(alpha_Y + beta_Yz*median(zt))  # same for all i
  for(t in 2:nt){
    Xt[t] <- Xt[t-1] + log(G[t]*Y[t] + (1-G[t])*S)
    }
  rbar[i] <- Xt[nt] - Xt[1]
  rcns[i] <- nt * log(Gcns*Ycns + (1-Gcns)*S)
  }

rdiff <- rbar-rcns

library(fields)
image.plot(iota_mu_seq,iota_sd_seq,rbar)
  # there is an optimum where everyone is doing the same thing (no variation)
  # but if the mean is too high (or too low), then better to have variation
  # more evidence that exitnction risk is crucial (since only then would variation help)
image.plot(iota_mu_seq,iota_sd_seq,rdiff)
  # - compared to constant environment, strategies with threshold above mean z 
  # are the ones that benefit most from the increase in variability
  # - even if mean threshold is far from optimum, increasing variance in thresholds forever
  # is never a good thing
  # - suggests that "having variance" and "good with variance" definitions of bet-hedging
  # aren't actually congruent

  # - my definition assumes that INDIVIDUALS are fully-responsive (plastic
  # to the environment, but the POPULATION as a whole can be unresponsive 
  # (if individuals vary a lot)
  # - so assumption is that plasticity is there to start with
  # - but assumption also implicit in Ackermann?
  # - some people would define bet-hedging as having ANY delay,
  # which here would mean a high germination threshold (variation could be small?)
  # - some population-level plasticity (threshold variation < max) always favoured here
  # because information in the signal (rainfall, germination, reproduction correlated)
  # - and minimal bet-hedging -> maximal plasticity favoured because:
  # (a) signal fully reliable
  # (b) no chance of extinction, so popsize variation doesn't matter
  # - random environmental variation in Y is mathematically equivalent to 
  # unreliable signal, becauses messes up relationship between G and Y
  # (or equivalently, between Y and z)

opt <- which(rbar==max(rbar))
iota_mu[opt]; iota_sd[opt]
curve(logistic(x,alpha_G[opt],beta_Gz[opt]),xlim=c(-3,3))

  # Other conclusions:
  # - mean threshold is proportional to intercept of logistic function,
  # so bet-hedgers aren't always the ones with a low germination probability

# Ceiling density-dependence ----------------------------------------------

Xt <- Xtc <- vector("numeric",length=nt) 
Xnmax <- 5 # max log number new seeds

rbar <- Xbar <- Xdiff <- matrix(NA,nr=np,nc=np) 

for(i in 1:np^2){
  Xt[1] <- Xtc <- 1
  S <- exp(-alpha_m)  # same for all i
  G <- plogis(alpha_G[i] + beta_Gz[i]*zt + eps_Y)
  Y <- exp(alpha_Y + beta_Yz*zt)  # same for all i
  Gc <- median(G)
  Yc <- median(Y)
  for(t in 2:nt){
    Xn <- Xt[t-1] + log(G[t]) + log(Y[t])
    Xn <- ifelse(Xn>Xnmax,Xnmax,Xn)
    Xo <- Xt[t-1] + log(1-G[t]) + log(S)
    Xt[t] <- log(exp(Xn)+exp(Xo))

    Xnc <- Xt[t-1] + log(Gc) + log(Yc)
    Xn <- ifelse(Xnc>Xnmax,Xnmax,Xnc)
    Xtc[t] <- Xtc[t-1] + log(exp(Xnc) + (1-Gc)*S)
    }
  rbar[i] <- Xt[nt] - Xt[1]
  Xbar[i] <- mean(Xt)
  Xdiff[i] <- Xbar[i] - mean(Xtc)
  }

library(fields)
image.plot(iota_mu_seq,iota_sd_seq,rbar)
image.plot(iota_mu_seq[iota_mu_seq<0.5],iota_sd_seq,rbar[iota_mu_seq<0.5,]) # zoomed
  # introducing DD creates need for variance in thresholds
  # (so that still have lots germinating in good years, 
  # but in usual good year, not EVERYONE germinates)
  # not sure about this: sensitive to conditions in final year

image.plot(iota_mu_seq,iota_sd_seq,Xbar)
image.plot(iota_mu_seq[iota_mu_seq<0.5],iota_sd_seq,Xbar[iota_mu_seq<0.5,]) # zoomed
  # geometric mean popsize

image.plot(iota_mu_seq,iota_sd_seq,Xdiff)
image.plot(iota_mu_seq[iota_mu_seq>0],iota_sd_seq,Xdiff[iota_mu_seq>0,]) # zoomed
  # diff in geometric mean popsize
  # again, not true to say that bet-hedgers benefit most
  # switch to variable environment can help or HINDER bet-hedgers, depending on mean threshold
  # if threshold too high, adding env var benefits those with fixed threshold than 
  # with variable threshold

opt <- which(rbar==max(rbar))
iota_mu[opt]; iota_sd[opt]
curve(logistic(x,alpha_G[opt],beta_Gz[opt]),xlim=c(-3,3))
  # Under certain DD, Bet-hedging can evolve even when cues reliable
  # (Reed papers, Furness et al 2015 Evolution)
  # (But already found by Ellner?)

# Gompertz density-dependence ---------------------------------------------

Xt <- Xtc <- vector("numeric",length=nt) 
rbar <- Xbar <- Xdiff <- matrix(NA,nr=np,nc=np) 
beta_Yd <- -0.1 # decline in Y with G density

for(i in 1:np^2){
  Xt[1] <- Xtc[1] <- 1
  S <- exp(-alpha_m)  # same for all i
  G <- plogis(alpha_G[i] + beta_Gz[i]*zt) # calculated ahead of time
  Gc <- plogis(alpha_G[i]) # mean(zt)=0
  for(t in 2:nt){
    Y <- exp(alpha_Y + beta_Yz*zt + beta_Yd*(log(G[t])+Xt[t-1]) + eps_Y)  # replaced every i
    Yc <- exp(alpha_Y + beta_Yd*(log(Gc)+Xt[t-1]) + eps_Y) # mean(zt)=0
    Xt[t] <- Xt[t-1] + log(G[t]*Y[t] + (1-G[t])*S)
    Xtc[t] <- Xtc[t-1] + log(Gc*Yc + (1-Gc)*S)
    }
  rbar[i] <- Xt[nt] - Xt[1]
  Xbar[i] <- mean(Xt)
  Xdiff[i] <- Xbar[i] - mean(Xtc)
  }

library(fields)
image.plot(iota_mu_seq,iota_sd_seq,rbar)
image.plot(iota_mu_seq[iota_mu_seq<1],iota_sd_seq,rbar[iota_mu_seq<1,]) # zoomed
  # Whether variance in thresholds benefits depends on form of DD:
  # Ceiling = yes
  # Gompertz = no
  # happens because ceiling = compensating, but Gompertz = over-compensating?
  # But note that adding DD changed location of optimum mean threshold
  # (towards higher threshold -> lower germination) 
  # results for Gompertz more consistent that those for ceiling DD?

image.plot(iota_mu_seq,iota_sd_seq,Xbar)
image.plot(iota_mu_seq[iota_mu_seq<0.5],iota_sd_seq,Xbar[iota_mu_seq<0.5,]) # zoomed
  # geometric mean popsize

image.plot(iota_mu_seq,iota_sd_seq,Xdiff)
image.plot(iota_mu_seq[iota_mu_seq>0],iota_sd_seq,Xdiff[iota_mu_seq>0,]) # zoomed
# diff in geometric mean popsize

# Ricker density-dependence -----------------------------------------------

Xt <- vector("numeric",length=nt) 
rbar <- Xbar <-matrix(NA,nr=np,nc=np) 
beta_Yd <- -0.0001 # decline in Y with G density

for(i in 1:np^2){
  Xt[1] <- 1
  S <- exp(-alpha_m)  # same for all i
  G <- plogis(alpha_G[i] + beta_Gz[i]*zt) # calculated ahead of time
  for(t in 2:nt){
    Y <- exp(alpha_Y + beta_Yz*zt + beta_Yd*G[t]*exp(Xt[t-1]))  # replaced every i
    Xt[t] <- Xt[t-1] + log(G[t]*Y[t] + (1-G[t])*S)
    }
  rbar[i] <- Xt[nt] - Xt[1]
  Xbar[i] <- mean(Xt)
  }

library(fields)
image.plot(iota_mu_seq,iota_sd_seq,rbar)
image.plot(iota_mu_seq[iota_mu_seq<1],iota_sd_seq,rbar[iota_mu_seq<1,]) # zoomed
  # Chaos, coding error, or rounding error?

image.plot(iota_mu_seq,iota_sd_seq,Xbar)
image.plot(iota_mu_seq[iota_mu_seq<0.5],iota_sd_seq,Xbar[iota_mu_seq<0.5,]) # zoomed
  # geometric mean popsize

# Variance efects on logistic ---------------------------------------------

zseq <- seq(-1,1,100)
plot(1,1,type="n",xlim=c(-3,3),ylim=c(0,1))
for(i in 1:np){
  curve(logistic(x,alpha_G[50,i],beta_Gz[50,i]),add=T,col=heat.colors(np)[i])
  }
  # non "bet-hedgers" = threshold (strong plasticity)
  # "bet-hedgers" = flatter (weak plasticity)
  # generalisation of ideas of ackermann etc.?

