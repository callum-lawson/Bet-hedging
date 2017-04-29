############################################################################
# Individual-level p(z) distributions given threshold germination function #
############################################################################

# Load packages -----------------------------------------------------------

library(truncdist)

# Choose parameters -------------------------------------------------------

nt <- 100
ni <- 100
np <- 50
# 20-21 mins

zmu <- 0
zsd <- 1

tau_mu_min <- -1
tau_mu_max <- 1
tau_sd_min <- 0.01
tau_sd_max <- 1

tau_mu_seq <- seq(tau_mu_min,tau_mu_max,length.out=np) # based on real data
tau_sd_seq <- seq(tau_sd_min,tau_sd_max,length.out=np) # based on real data
tau_mu <- tau_sd <- matrix(NA,nr=np,nc=np) 
tau_mu[] <- rep(tau_mu_seq,times=np)
tau_sd[] <- rep(tau_sd_seq,each=np)
alpha_Y <- 1
beta_Yz <- 1
sigma_Y <- 1
alpha_m <- 1

alpha_G <- beta_Gz <- matrix(NA,nr=np,nc=np) 
alpha_G[] <- -tau_mu*sqrt(pi^2/(3*tau_sd^2))
beta_Gz[] <- sqrt(pi^2/(3*tau_sd^2))
# based on Godfray & Rees 2002

set.seed(1)
zt_mat <- eps_Yv <- array(NA,dim=c(nt,ni))
zt_mat[] <- rnorm(nt*ni,mean=zmu,sd=zsd)
zc <- rep(0,nt)
eps_Yc <- rnorm(nt,0,0)
eps_Yv[] <- rnorm(nt*ni,0,sigma_Y)
X0 <- 5

# Calculate individual environment distributions --------------------------

sdmult <- 3
zseq <- seq(zmu-sdmult*zsd,zmu+sdmult*zsd,length.out=ni)

nq <- 100
tau_full_seq <- seq(tau_mu_min-sdmult*tau_sd_max,tau_mu_max+sdmult*tau_sd_max,length.out=nq)
pgivetau <- matrix(nr=ni,nc=nq)

i <- 5
j <- 5
ptau <- dnorm(tau_full_seq,tau_mu[i,j],tau_sd[i,j]) # picking out of matrix
for(k in 1:nq){
  pgivetau[,k] <- dtrunc(x=zseq, spec="norm", mean=zmu, sd=zsd, a=tau_full_seq[k], b=Inf)
  }
pgivetau_w <- rep(ptau,each=nq) * pgivetau # weighted
pz <- rowSums(pgivetau_w)

par(mar=c(5,5,1,1))
plot(pz~zseq,type="l")
library(evd)
curve(dgumbel(x, loc=tau_mu[i,j], scale=tau_sd[i,j]),add=T,col="red")
curve(pnorm(x, mean=tau_mu[i,j], sd=tau_sd[i,j]),add=T,col="red")

# Analytical solution - Gumbel --------------------------------------------

# pnorm gives proportion of pop that will germinate at different z values






