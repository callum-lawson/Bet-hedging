############################################################################
# Simulations to investigate relationships between environmental variance, #
# stochastic growth rate, and mean time to extinction                      #
############################################################################

nt <- 1000
nsim <- 50

logN0 <- rep(0,nsim)		# N0 = 1
rmean <- rep(-10^-1,nsim)
rsd <- rep(1,nsim)

logN <- r <- matrix(nr=nt,nc=nsim)
for(i in 1:nsim){
	r[,i] <- rnorm(nt,mean=rmean[i],sd=rsd[i]) # change to i if needed
	for(t in 1:nt){
		if(t==1) logN[t,i] <- logN0[i] + r[t,i]
		if(t!=1) logN[t,i] <- logN[t-1,i] + r[t,i]
		}
	}

tau <- 0 			# log scale

matplot(1:nt,logN,type="l",lty=1,col="black")
abline(h=tau,col="red",lty=2)

extimes <- apply(logN,2,function(N) which(N<tau)[1])
