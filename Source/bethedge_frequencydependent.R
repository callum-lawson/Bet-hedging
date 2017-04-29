###############################################################################
# Invasion analyses for bet-hedgers with two phenotypes that have independent #
# carrying capacities                                                         #
###############################################################################

nt <- 100
q <- 0.5
z <- rbinom(nt,prob=q,size=1)
W1 <- c(1/2,3)
W2 <- c(1/2,1/2)

# Separate DD

pmin <- 0
pmax <- 1
np <- 50
pseq <- seq(pmin,pmax,length.out=np)

fcalc <- function(z,f,p,q,W1,W2,v){
  if(p>0 | q>0){
    if(z==0) f1 <- (f*p*W1[1]) / (f*p*W1[1] + (1-f)*q*W1[1])
    if(z==1) f1 <- (f*p*W1[2]) / (f*p*W1[2] + (1-f)*q*W1[2])
    }
  if(p<1 | q<1){
    if(z==0) f2 <- (f*(1-p)*W2[1]) / (f*(1-p)*W2[1] + (1-f)*(1-q)*W2[1])
    if(z==1) f2 <- (f*(1-p)*W2[2]) / (f*(1-p)*W2[2] + (1-f)*(1-q)*W2[2])
    }
  if(p==0 & q==0){
    return(f2)
    }
  if(p==1 & q==1){ 
    return(f1)
    }
  if(!(p==0 & q==0) & !(p==1 & q==1)){
    return(v*f1 + (1-v)*f2)
    }
  }
  # v = carrying capacity for phenotype 1 as proportion of total pop
  # (if K of one type twice as high then 
  # a given increase in freq in that type is worth twice as much)

farr <- array(dim=c(nt,np,np))
farr[1,,] <- 0.001 # starting at 1% frequency

for(i in 1:np){
  for(j in 1:np){
    for(t in 2:nt){
      farr[t,i,j] <- fcalc(z[t],farr[t-1,i,j],pseq[i],pseq[j],W1,W2,v=0.25)
      }
    }
  }

library(fields)
image.plot(pseq,pseq,t(farr[nt,,]),xlab="resident",ylab="invader")

### SIMS

W1 <- c(4,1/2)
W2 <- c(1/2,1/2)

DDgomp <- function(N,NT,alpha=0.5,beta=-0.1){
  exp(log(N) + alpha + beta*log(NT))
  }

alpha1 <- 0.25
alpha2 <- 1
beta <- -0.1 # undercompensating
(K = exp(-alpha/beta))

DDFUN1 <- DDFUN2 <- DDgomp

Np1 <- Nq1 <- Np2 <- Nq2 <- N1 <- N2 <- NpT <- NqT <- array(dim=c(nt,np,np))
NpT[1,,] <- 0.001
NqT[1,,] <- 5 *10^3
p <- 0.5; q <- 0.5

for(i in 1:np){
  for(j in 1:np){
    p <- pseq[i]
    q <- pseq[j]
    for(t in 2:nt){
      if(z[t]==0){
        Np1[t,i,j] <- NpT[t-1,i,j]*p*W1[1]
        Nq1[t,i,j] <- NqT[t-1,i,j]*q*W1[1]
        Np2[t,i,j] <- NpT[t-1,i,j]*(1-p)*W2[1]
        Nq2[t,i,j] <- NqT[t-1,i,j]*(1-q)*W2[1]
        }
      if(z[t]==1){
        Np1[t,i,j] <- NpT[t-1,i,j]*p*W1[2]
        Nq1[t,i,j] <- NqT[t-1,i,j]*q*W1[2]
        Np2[t,i,j] <- NpT[t-1,i,j]*(1-p)*W2[2]
        Nq2[t,i,j] <- NqT[t-1,i,j]*(1-q)*W2[2]
        }
             
      N1[t,i,j] <- Np1[t,i,j] + Nq1[t,i,j]
      N2[t,i,j] <- Np2[t,i,j] + Nq2[t,i,j]
      NTT <-  N1[t,i,j] + N2[t,i,j]
      
      # NpT[t,i,j] <- DDFUN1(Np1[t,i,j],N1[t,i,j],alpha1) + DDFUN2(Np2[t,i,j],N2[t,i,j],alpha2)
      # NqT[t,i,j] <- DDFUN1(Nq1[t,i,j],N1[t,i,j],alpha1) + DDFUN2(Nq2[t,i,j],N2[t,i,j],alpha2)
      
      NpT[t,i,j] <- DDFUN1(Np1[t,i,j],NTT,alpha1) + DDFUN2(Np2[t,i,j],NTT,alpha1)
      NqT[t,i,j] <- DDFUN1(Nq1[t,i,j],NTT,alpha1) + DDFUN2(Nq2[t,i,j],NTT,alpha1)
      }
    }
  }

fp <- NpT/(NpT + NqT)
image.plot(pseq,pseq,t(fp[nt,,]),xlab="resident",ylab="invader")

