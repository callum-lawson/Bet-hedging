####################################################################
# Plots of r~z for pure types and a 50% mix phenotype, for convex, #
# concave, and mixed response types                                #
####################################################################

# maindir <- "D:/Users/calluml/Dropbox/NIOO/"
# maindir <- "C:/Users/Callum/Documents/My Dropbox/NIOO/"

setwd(paste0(maindir,"Analyses/Bet-hedging"))

zmin <- -5
zmax <- 5
nz <- 10^4
zseq <- seq(zmin,zmax,length.out=nz)

zmu <- 0 # -0.5 
zsd <- 1
dp <- 2

lamfunc <- function(b0=0,b1,b2){
  exp(b0 + b1*zseq + b2*zseq^2)
  }

moments <- function(lam){
  pz <- dnorm(zseq,mean=zmu,sd=zsd)
  lam0 <- lam[which(pz==max(pz))[1]]
  lambar <- sum(pz*lam)/sum(pz)
  lamsd <- sqrt(sum(pz*(lam-lambar)^2)/sum(pz))
  r <- log(lam)
  #pz_r <- pz[!is.nan(r)]
  #r <- r[!is.nan(r)]
  rbar <- sum(pz*r)/sum(pz)
  rdelta <- rbar - log(lam0)
  c(lambar=lambar,lamsd=lamsd,rbar=rbar) # lam0=lam0
  }

threematplot <- function(mat,...){
  matplot(zseq,mat,
    type="l",
    lty=c(1,3,2),
    col=c("black","red","black"),
    xlab="",
    ylab="",
    ...
    )
  }

jointrplot <- function(b0a=0,b0b=0,b1a,b1b,b2a,b2b){
  require(gplots)
  y1 <- lamfunc(b0a,b1a,b2a)
  y2 <- lamfunc(b0b,b1b,b2b)
  ym <- apply(cbind(y1,y2),1,mean)
  lammat <- cbind(y1,ym,y2)
  rmat <- log(lammat)
  threematplot(lammat)
  threematplot(rmat)
  stats <- round(apply(lammat,2,moments),dp) 
  textplot(stats)
  }

pseq <- seq(0,1,length.out=10)

levinsplot <- function(b0a=0,b0b=0,b1a,b1b,b2a,b2b){
  require(gplots)
  y1 <- lamfunc(b0a,b1a,b2a)
  y2 <- lamfunc(b0b,b1b,b2b)
  yp1 <- apply(cbind(y1,y2),1,function(x) pseq*x[1])
  yp2 <- apply(cbind(y1,y2),1,function(x) (1-pseq)*x[2])
  matplot(t(yp1[,]),t(yp2[,]),type="l")
  }

levinsplot(b1a=0.1,b1b=-0.1,b2a=0,b2b=0)

pdf(paste0("bethedge_growthcurves_",format(Sys.Date(),"%d%b%Y"),".pdf"),width=7,height=7)
par(mfrow=c(3,3),mar=c(4,4,1,1),bty="l",las=1)

jointrplot(b1a=0.25,b1b=-0.25,b2a=0,b2b=0)
jointrplot(b1a=0.5,b1b=-0.5,b2a=0.05,b2b=0.05)
jointrplot(b1a=0.25,b1b=-0.25,b2a=-0.05,b2b=-0.05)

jointrplot(b1a=0.5,b1b=0,b2a=0,b2b=0)
jointrplot(b1a=1,b1b=0,b2a=-0.1,b2b=0)
jointrplot(b1a=0.5,b1b=0,b2a=0.05,b2b=0)

jointrplot(b1a=0,b1b=0,b2a=0.1,b2b=-0.1)
jointrplot(b1a=0,b1b=0,b2a=0.1,b2b=0)
jointrplot(b1a=0,b1b=0,b2a=-0.1,b2b=0)

jointrplot(b1a=1.5,b1b=-1.5,b2a=-0.25,b2b=-0.25)
jointrplot(b1a=1,b1b=-1,b2a=-0.25,b2b=-0.25)
jointrplot(b1a=0.5,b1b=-0.5,b2a=-0.25,b2b=-0.25)

jointrplot(b0a=0.25,b0b=-0.25,b1a=0,b1b=0,b2a=-0.5,b2b=0)
jointrplot(b0a=-1,b0b=1,b1a=0,b1b=0,b2a=0.1,b2b=-0.1)

jointrplot(b0a=-0.5,b0b=0.5,b1a=0.5,b1b=0,b2a=0,b2b=0)
  # pure sensitive -> mixed
jointrplot(b0a=0.5,b0b=-0.5,b1a=0.5,b1b=0,b2a=0,b2b=0)
  # pure insensitive -> mixed

dev.off()

opt(0.25,-0.05)

opt <- function(b1,b2){ 
  -b1/(2*b2) 
  }

# Ellner & Hairston -------------------------------------------------------

