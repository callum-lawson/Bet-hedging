#####################################################################################
# Calculation of environmental variance effects, including threshold variance for   #
# bet-hedging, for mixes of flat and convex, linear, and concave response functions #
#####################################################################################

# Set working directories -------------------------------------------------

# maindir <- "D:/Users/calluml/Dropbox/NIOO/"
# maindir <- "C:/Users/Callum/Dropbox/NIOO/"

setwd(paste0(maindir,"Analyses/Bet-hedging"))

# Functions ---------------------------------------------------------------

lamnorm_f <- function(z,mu,sig,gam){
  gam*dnorm(x=z,mean=mu,sd=sig)
  }

lamexp_f <- function(z,b0=0,b1,b2){
  exp(b0 + b1*z + b2*z^2)
  }

mix_f <- function(a,b,p){ 
  (1-p)*a + p*b
  }

pz_f <- function(zseq,zmu,zsd){
  pz <- dnorm(zseq,mean=zmu,sd=zsd)
  if(Inf %in% pz) pz[pz==Inf] <- 1
    # doesn't matter because re-scaled anyway
  return(pz)
  }

rbar_f <- function(pz,r){
  pp <- pz>0
  sum(pz[pp]*r[pp])/sum(pz[pp])
  }

opt <- function(b1,b2){ 
  -b1/(2*b2) 
  }

# Inputs ------------------------------------------------------------------

### General

nz <- 1000 + 1 # +1 -> 0 included
zmin <- -3
zmax <- 3
zseq <- seq(zmin,zmax,length.out=nz)

np <- 100
pseq <- seq(0,1,length.out=np)

zmu <- 0
zsdmin <- 0
zsdmax <- 1
nzsd <- 100
zsdseq <- seq(zsdmin,zsdmax,length.out=nzsd)

### Normal

alpha <- 10 # Sum[A(thin)]
betaseq <- c(0.5,1,2) # c(2,4,8) # A(thin)/A(narrow)
nbeta <- length(betaseq)
gamseq <- alpha*betaseq

mumin <- 0
mumax <- 2
nmu <- 3
museq <- seq(mumin,mumax,length.out=nmu) # 2sd, 1sd, 0sd
mucons <- 0
sigseq <- c(1,10)

### Exponential

zsdmax <- 2
nzsd <- 100
zsdseq <- seq(zsdmin,zsdmax,length.out=nzsd)

b0seq <- c(-0.5,0,0.5)

b1seq <- c(1,1,1)
b2seq <- c(-0.5,0,0.2)

nb <- 3

# Calculations - Exponential ----------------------------------------------

pmax <- array(dim=c(nb,nb,nzsd))
l1 <- l2 <- array(dim=c(nz,nb,nb))

for(i in 1:nb){
  for(j in 1:nb){
    for(k in 1:nzsd){
      
      l1[,i,j] <- lamexp_f(z=zseq,b1=b1seq[j],b2=b2seq[j])
      l2[,i,j] <- rep(exp(b0seq[i]),nz)

      lm <- sapply(pseq,function(p) mix_f(l1[,i,j],l2[,i,j],p))
      rm <- log(lm)
      
      pz <- pz_f(zseq,zmu,zsdseq[k])
      rbarm <- apply(rm,2,function(r) rbar_f(pz,r))
      pmaxpos <- which(rbarm==max(rbarm))
      pmax[i,j,k] <- ifelse(length(pmaxpos)==1, pseq[pmaxpos], NA)
      
      }
    }
  }


# Plots - Exponential -----------------------------------------------------

library(abind)
l12 <- abind(l1,l2,along=4)

pdf(paste0("bethedging_localanalyses_nonlinear_",
  format(Sys.Date(),"%d%b%Y"),".pdf"),
  width=12,height=8
  )

layout(laymat)
par(mar=c(4,4,2,2),las=1,bty="l")
for(i in 1:nbeta){
  for(j in 1:nmu){
    matplot(zseq,log(l12[,i,j,]),type="l",xlim=c(-2,2),
      xlab="z",ylab="r")
    abline(v=0,col="blue",lty=3)
    
    plot(pmax[i,j,]~zsdseq,type="l",ylim=c(0,1),
      xlab=expression(sigma[z]),ylab="p")
    }
  }

dev.off()

# Calculations - Normal ---------------------------------------------------

### !!! Needs to be edited to match Exponential functions, e.g. for pmax !!!

# for each scale value (gamseq),
# for each curve position(museq),
# for each zsd (zsdseq),
# calculate optimum p (pmax)

pmax <- array(dim=c(nbeta,nmu,nzsd))
l1 <- l2 <- array(dim=c(nz,nbeta,nmu))

for(i in 1:nbeta){
  for(j in 1:nmu){
    for(k in 1:nzsd){
      
      l1[,i,j] <- lamnorm_f(z=zseq,mu=museq[j],sig=sigseq[1],alpha)
      # l2[,i,j] <- lamnorm_f(z=zseq,mu=museq[j],sig=sigseq[2],gamseq[i])
      l2[,i,j] <- betaseq[i]*rep(lamnorm_f(z=0,mu=museq[j],sig=sigseq[1],alpha),nz)
        # z=0 -> growth rate in mean environment
      
      lm <- sapply(pseq,function(p) mix_f(l1[,i,j],l2[,i,j],p))
      rm <- log(lm)
      
      pz <- pz_f(zseq,zmu,zsdseq[k])
      rbarm <- apply(rm,2,function(r) rbar_f(pz,r))
      pmax[i,j,k] <- pseq[which(rbarm==max(rbarm))]
      
      }
    }
  }
  

# Plots - Normal ----------------------------------------------------------

laymat <- cbind(
  rep(c(1,7,13),each=2),
  rep(c(1,7,13),each=2),
  c(2,0,8,0,14,0),
  rep(c(3,9,15),each=2),
  rep(c(3,9,15),each=2),
  c(4,0,10,0,16,0), 
  rep(c(5,11,17),each=2),
  rep(c(5,11,17),each=2),
  c(6,0,12,0,18,0)
  )
  
library(abind)
l12 <- abind(l1,l2,along=4)

pdf(paste0("bethedging_localanalyses_normal_",
  format(Sys.Date(),"%d%b%Y"),".pdf"),
  width=12,height=8
  )

layout(laymat)
par(mar=c(4,4,2,2),las=1)
for(i in 1:nbeta){
  for(j in 1:nmu){
    matplot(zseq,log(l12[,i,j,]),type="l",xlim=c(-2,2),
      xlab="z",ylab="r")
    abline(v=0,col="blue",lty=3)
    
    plot(pmax[i,j,]~zsdseq,type="l",ylim=c(0,1),
      xlab=expression(sigma[z]),ylab="p")
    }
  }
  
dev.off()

curve(dnorm(x,zmu,zsdmax,log=T),xlim=c(zmin,zmax))




