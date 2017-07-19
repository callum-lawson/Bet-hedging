#######################################################################################
# How do the reliability and discrimination ability of a signal relate to each other? #
#######################################################################################

require(MASS)

zmu <- 0
wmu <- 0 # irrelevant
zsd <- 0.1
wsd <- 1
rho <- 0.8
  
nsim <- 10^6
zw_mu <- c(zmu,wmu)
set.seed(1)
zw_sig <- matrix(c(zsd^2,rep(rho*zsd*wsd,2),wsd^2),nr=2,nc=2)
zw <- mvrnorm(n=nsim, mu=zw_mu, Sigma=zw_sig)
colnames(zw) <- c("z","w")
zw <- as.data.frame(zw)

zw$wcat <- cut(zw$w,breaks=100)
byw <- with(zw,split(z,wcat))

plot(density(byw[[50]]))
lines(density(byw[[40]]),col="red")
lines(density(byw[[60]]),col="blue")
  # Increasing zsd only scales z up and down - doesn't change discrimination ability
