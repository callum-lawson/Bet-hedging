############################################################################
# Is the distribution of log(lambda) altered by demographic stochasticity? #
############################################################################

lambda <- 2^(0:3)
N <- 2^(0:5)

nlambda <- length(lambda)
nN <- length(N)
nrep <- 10^4

i <- 1
j <- 1

ybar <- array(NA,dim=c(nrep,nlambda,nN))

for(i in 1:nlambda){
  for(j in 1:nN){
    repseq <- rep(1:nrep,each=N[j])
    ysim <- rpois(nrep*N[j],lambda[i])
    ybar[,i,j] <- tapply(ysim,repseq,mean)
    }
  }

ybar[ybar==0] <- 0.5

par(mfrow=c(2,2),mar=c(3,3,0.5,0.5))
for(i in 1:nlambda){
  boxplot.matrix(log(ybar[,i,]),labels=N)
  }
