#######################################################################################
# How do the reliability and discrimination ability of a signal relate to each other? #
#######################################################################################

require(MASS)

zwdist <- function(zmu,wmu,zsd,wsd,rho=0.8,nsim=10^6,...){
  zw_mu <- c(zmu,wmu)
  zw_sig <- matrix(c(zsd^2,rep(rho*zsd*wsd,2),wsd^2),nr=2,nc=2)
  zw <- mvrnorm(n=nsim, mu=zw_mu, Sigma=zw_sig)
  colnames(zw) <- c("z","w")
  zw <- as.data.frame(zw)
  
  zw$wcat <- cut(zw$w,...)
  byw <- with(zw,split(z,wcat))
  
  list(zw=zw,byw=byw)
}

mybreaks <- c(-1.05,-0.95,-0.05,0.05,0.95,1.05)

# Imperfect information ---------------------------------------------------

set.seed(100)
zw1 <- zwdist(zmu=0,wmu=0,zsd=0.1,wsd=1,breaks=mybreaks)
zw2 <- zwdist(zmu=0,wmu=0,zsd=0.1*2,wsd=1*2,breaks=mybreaks)
zw3 <- zwdist(zmu=-0.2,wmu=-2,zsd=0.1,wsd=1,breaks=mybreaks)
zw4 <- zwdist(zmu=-0.2,wmu=-2,zsd=0.1*2,wsd=1*2,breaks=mybreaks)

par(mfrow=c(1,1))
plot(density(zw1$byw[[1]]),lty=1,col="black",main="")
lines(density(zw1$byw[[3]]),lty=2,col="black")
lines(density(zw1$byw[[5]]),lty=3,col="black")
lines(density(zw2$byw[[1]]),lty=1,col="blue")
lines(density(zw2$byw[[3]]),lty=2,col="blue")
lines(density(zw2$byw[[5]]),lty=3,col="blue")
lines(density(zw3$byw[[1]]),lty=1,col="red")
lines(density(zw3$byw[[3]]),lty=2,col="red")
lines(density(zw3$byw[[5]]),lty=3,col="red")
lines(density(zw4$byw[[1]]),lty=1,col="orange")
lines(density(zw4$byw[[3]]),lty=2,col="orange")
lines(density(zw4$byw[[5]]),lty=3,col="orange")
# Increasing zsd only scales z up and down - doesn't change discrimination ability

plot(z~w,data=zw2$zw[1:1000,],col="blue")
points(z~w,data=zw1$zw[1:1000,],col="black")
points(z~w,data=zw3$zw[1:1000,],col="red")
points(z~w,data=zw4$zw[1:1000,],col="orange")

# Perfect information -----------------------------------------------------

zw1b <- zwdist(zmu=0,wmu=0,zsd=0.1,wsd=1,breaks=mybreaks,rho=1)
zw2b <- zwdist(zmu=0,wmu=0,zsd=0.2,wsd=2,breaks=mybreaks,rho=1)
zw3b <- zwdist(zmu=-0.1,wmu=-1,zsd=0.1,wsd=1,breaks=mybreaks,rho=1)

plot(z~w,data=zw2b$zw[1:1000,],col="blue")
points(z~w,data=zw1b$zw[1:1000,],col="black")
points(z~w,data=zw3b$zw[1:1000,],col="red")
