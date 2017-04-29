###########################################################################
# Plot r~z and lambda~z for range of G values to test congruence between  #
# different bet-hedging definitions                                       #
###########################################################################

# Libraries ---------------------------------------------------------------

library(fields) # for tim.colors

# Functions ---------------------------------------------------------------

Lf <- function(z,G,S,YFUN){
  G*YFUN(z) + (1-G)*S
  }
  # Y functions defined below

estmu <- function(py,y){
  sum(py*y)/sum(py)
  }

estsd <- function(py,y){
  E_Y2 <- estmu(py,y^2)
  EY_2 <- estmu(py,y)^2
  sqrt(E_Y2-EY_2)
  # var(Y) = E[Y^2] - E[Y]^2
  }

Gmatplot <- function(y,...){
  ncol <- ncol(t(y))
  matplot(Gseq,t(y),type="l",lty=1, col=tim.colors(ncol)[1:ncol],...)
  }

# plot responses

Gcurveplot <- function(YFUN,log=T,...){
  gmat <- matrix(nr=nz,nc=nG)
  if(log==F) gmat[] <- sapply(Gseq,function(x) Lf(zseq,G=x,S,YFUN))
  if(log==T) gmat[] <- sapply(Gseq,function(x) log(Lf(zseq,G=x,S,YFUN)))
  matplot(zseq,gmat,col=Gcols,type="l",lty=1,
    xlab="z",ylab=ifelse(log==F,expression(lambda),"r")
    )
  }

# calculate growth statistics (for each combo of mu,sd,G)

growcalc <- function(YFUN){
  
  Lz <- array(dim=c(nz,nmu,nsd,nG))
  Lmu <- Lsd <- rmu <- rsd <- array(dim=c(nmu,nsd,nG))
  
  Lz[] <- Lf(zseq_full,Gseq_full,S,YFUN)
  
  for(i in 1:nmu){
    for(j in 1:nsd){
      for(k in 1:nG){
        Lmu[i,j,k] <- estmu(pz[,i,j,k],Lz[,i,j,k])
        Lsd[i,j,k] <- estsd(pz[,i,j,k],Lz[,i,j,k])
        rmu[i,j,k] <- estmu(pz[,i,j,k],log(Lz[,i,j,k]))
        rsd[i,j,k] <- estsd(pz[,i,j,k],log(Lz[,i,j,k]))
        }
      }
    }
  
  d_Lmu <- Lmu[,nsd,] - Lmu[,1,]
  d_Lsd <- Lsd[,nsd,] - Lsd[,1,]
  d_rmu <- rmu[,nsd,] - rmu[,1,] 
  d_rsd <- rsd[,nsd,] - rsd[,1,] 
  # additive
  
  m_Lmu <- Lmu[,nsd,] / Lmu[,1,]
  m_Lsd <- Lsd[,nsd,] / Lsd[,1,]
  m_rsd <- rsd[,nsd,] / rsd[,1,]
  # multiplicative
  
  list(Lz=Lz,
    Lmu=Lmu,Lsd=Lsd,rmu=rmu,rsd=rsd,
    d_Lmu=d_Lmu,d_Lsd=d_Lsd,d_rmu=d_rmu,d_rsd=d_rsd,
    m_Lmu=m_Lmu,m_Lsd=m_Lsd,m_rsd=m_rsd
    )
  
  }

# key plots

betplots <- function(g,YFUN){
  par(mfrow=c(2,3),bty="l",las=1,mar=c(4.5,4.5,1.5,1.5))
  Gcurveplot(YFUN,log=F)
  Gcurveplot(YFUN,log=T)
  
  Gmatplot(g$rmu[mupos,,],xlab="G",ylab="r")

  plot(g$d_rmu[mupos,]~g$Lsd[mupos,sdpos,],col=Gcols,
    xlab=expression(sigma[lambda]^2),
    ylab=expression(delta~bar(r))
    )
  
  plot(g$d_rmu[mupos,]~I(g$Lsd[mupos,sdpos,]/(2*g$Lmu[mupos,sdpos,])),col=Gcols,
    xlab=expression(sigma[lambda]^2/(2*bar(lambda))),
    ylab=expression(delta~bar(r))
    )
  
  plot(g$Lmu[mupos,sdpos,]~g$Lsd[mupos,sdpos,],col=Gcols,
    xlab=expression(sigma[lambda]^2),ylab=expression(bar(lambda))
    )
  }

# Parameters --------------------------------------------------------------

nz <- 10^3
zmin <- -5
zmax <- 5
zseq <- seq(zmin,zmax,length.out=nz)

nmu <- 5
zmu_min <- -1
zmu_max <- 1
zmuseq <- seq(zmu_min,zmu_max,length.out=nmu)

nsd <- 3
zsd_min <- 0.1
zsd_max <- 1
zsdseq <- seq(zsd_min,zsd_max,length.out=nsd)

nG <- 100
Gmin <- 0
Gmax <- 1
Gseq <- seq(Gmin,Gmax,length.out=nG)

S <- 1.25 # !!! CHANGED !!!
ntot <- nz*nmu*nsd*nG

# Generate climates -------------------------------------------------------

zseq_full <- zmuseq_full <- zsdseq_full <- Gseq_full <- array(dim=c(nz,nmu,nsd,nG))

zseq_full[] <- rep(zseq,times=nmu*nsd*nG)
zmuseq_full[] <- rep(rep(zmuseq,each=nz),times=nsd*nG)
zsdseq_full[] <- rep(rep(zsdseq,each=nz*nmu),times=nG)
Gseq_full[] <- rep(Gseq,each=nz*nmu*nsd)

pz <- array(dim=c(nz,nmu,nsd,nG))
pz[] <- dnorm(zseq_full,zmuseq_full,zsdseq_full)

# Calculate growth rates --------------------------------------------------

Y_exp <- function(z,beta=0.5){
  exp(beta*z)
  }

Y_expsq <- function(z,alpha=1,beta=0,gamma=-0.1){
  exp(alpha+beta*z+gamma*z^2)
  }

Y_quad <- function(z,alpha=1,beta=2,gamma=-0.1){
  exp(alpha+beta*z+gamma*z^2)
  }

Y_cohen <- function(z,tau=0,Y=3){
  ifelse(z<tau,0,Y)
  }

Y_lin <- function(z,alpha=3,beta=0.05){
  exp(alpha+beta*z)
  }

Y_sq <- function(z,beta=0.5){
  2^(beta*z)
  }

Y_sat <- function(z,alpha=1,beta=0.5){
  x <- exp(z)
  alpha*x/(beta+x)
  }

g_exp <- growcalc(Y_exp)
g_expsq <- growcalc(Y_expsq)
g_quad <- growcalc(Y_quad)
g_cohen <- growcalc(Y_cohen)
g_lin <- growcalc(Y_lin)
g_sq <- growcalc(Y_sq)
g_sat <- growcalc(Y_sat)

# Plot results ------------------------------------------------------------

mupos <- 4
sdpos <- 2
Gpos <- 1

Gcols <- tim.colors(nG)[1:nG] # red = higher G

# 1. G (trait)
# 2. Lsd (classic)
# 3. d_rmu (change)

Y_quad <- function(z,alpha=1,beta=2,gamma=-0.1){
  exp(alpha+beta*z+gamma*z^2)
  }
g_quad <- growcalc(Y_quad)
betplots(g_quad,Y_quad)
Gmatplot(g_quad$rmu[mupos,,])

betplots(g_exp,Y_exp)
betplots(g_expsq,Y_expsq)
betplots(g_quad,Y_quad)
betplots(g_cohen,Y_cohen)
betplots(g_lin,Y_lin)
betplots(g_sq,Y_sq)
betplots(g_sat,Y_sat)

### G plots

# Raw outputs

par(mfrow=c(1,1))
Gmatplot(g_exp$Lmu[mupos,,])
Gmatplot(Lsd[mupos,,])
Gmatplot(Lsd[mupos,,]/(2*Lmu[mupos,,]))
Gmatplot(g_exp$rmu[mupos,,])
Gmatplot(rsd[mupos,,])

Gmatplot(Lmu[,sdpos,])
Gmatplot(Lsd[,sdpos,])
Gmatplot(Lsd[,sdpos,]/(2*Lmu[,sdpos,]))
Gmatplot(rsd[,sdpos,])

# Changes in variable environment

Gmatplot(d_Lmu[,])
Gmatplot(d_Lsd[,])
Gmatplot(d_rmu[,],ylab=expression(delta~bar(r)))
Gmatplot(d_rsd[,])

Gmatplot(m_Lmu[,])
Gmatplot(m_Lsd[,])
Gmatplot(m_rsd[,])

### Change Vs lambda variance

plot(d_rmu[mupos,]~Lsd[mupos,1,],col=Gcols,
  xlab=expression(sigma[lambda]^2),
  ylab=expression(delta~bar(r))
  )

plot(d_rmu[mupos,]~Lmu[mupos,1,],col=Gcols,
  xlab=expression(bar(lambda)),
  ylab=expression(delta~bar(r))
  )

plot(d_rmu[mupos,]~I(Lsd[mupos,1,]/(2*Lmu[mupos,1,])),col=Gcols,
  xlab=expression(sigma[lambda]^2/(2*bar(lambda))),
  ylab=expression(delta~bar(r))
  )

plot(Lmu[mupos,sdpos,]~Lsd[mupos,sdpos,],col=Gcols)

### Summary

# G != Lmu (analytical; known?)
# G = Lsd (known)
# G != d_rmu (known?)
# Lsd != d_rmu (new)

# Invasion analyses -------------------------------------------------------

rdiffplot <- function(YFUN){
  par(mfrow=c(1,1))
  plot(I(log(Lf(zseq,G=0.5,S,YFUN))-log(Lf(zseq,G=1,S,YFUN)))~zseq,type="l")
  lines(I(log(Lf(zseq,G=0.5,S,YFUN))-log(Lf(zseq,G=0,S,YFUN)))~zseq,lty=2)
  abline(h=0,col="red",lty=3)
  }

rdiffplot(Y_exp)
  # area over 0 is smaller than area under 0
  # (for certain mean z)
  # maximum fitness is gained by tracing along top lines 
  # (full information / plasticity)
betplots(g_exp,Y_exp)

rdiffplot(Y_expsq)

# Intuitive scenarios -----------------------------------------------------

expscen <- rbind(c(1,2,4),rev(c(1,2,4)))
linscen <- rbind(c(1,2,3),rev(c(1,2,3)))
satscen <- rbind(c(1,3,4),rev(c(1,3,4)))

prod(expscen[1,])
prod(apply(expscen,2,mean))

mean(expscen[1,])
mean(apply(expscen,2,mean))

sd(expscen[1,])
sd(apply(expscen,2,mean))


prod(linscen[1,])
prod(apply(linscen,2,mean))

mean(linscen[1,])
mean(apply(linscen,2,mean))

sd(linscen[1,])
sd(apply(linscen,2,mean))


prod(satscen[1,])
prod(apply(satscen,2,mean))

mean(satscen[1,])
mean(apply(satscen,2,mean))

sd(satscen[1,])
sd(apply(satscen,2,mean))



madscen <- rbind(c(1,3,4),c(0.25,3,15))

prod(madscen[1,])
prod(madscen[2,])
prod(apply(madscen,2,mean))

mean(madscen[1,])
mean(madscen[2,])
mean(apply(madscen,2,mean))

sd(madscen[1,])
sd(madscen[2,])
sd(apply(madscen,2,mean))
  
plot(c(2,16))
points(c(4,8),pch=2)
  # (4,8)->(3,12): 0.75 < 1.50
  # (2,16)->(3,12): 1.50 > 0.75
2*16; 4*8; 3*12
  # mixed -> more variance, 
mean(c(2,16))
mean(c(4,8))
mean(c(3,12))

# two convex

a <- 2^(1:3)
b <- 2^(c(-6,2,10))
c <- apply(rbind(a,b),2,mean)
m <- a/b
log(m[3])+log(m[1])
matplot(cbind(a,b,c),pch=1)
prod(a)
prod(b)
prod(c)
(2^2)^3
  # var env -> mixing improves, others don't
log2(c)
plot(a,b)
cor.test(a,b)

p <- seq(0,1,length.out=100)
f1 <- log(p*a[1]+(1-p)*b[1])
f2 <- log(p*a[3]+(1-p)*b[3])
fm <- f1 + f2
matplot(qlogis(p),cbind(f1,f2),type="l",lty=1)
points(f2~qlogis(p),col="red")
plot(exp(f1)~qlogis(p))
plot(exp(f2)~p)
plot(fm~p)

# convex / concave

a <- 2^(1:3)
b <- 2^(c(1.5,3,1.5))
c <- apply(rbind(a,b),2,mean)
matplot(cbind(a,b,c),pch=1)
prod(a)
prod(b)
prod(c)
a[2]^3
b[2]^3
c[2]^3

# convex / flat 

a <- 2^(1:3)
b <- 2^c(2,2,2)
c <- apply(rbind(a,b),2,mean)
matplot(cbind(a,b,c),pch=1)
matplot(log(cbind(a,b,c)),pch=1)
prod(a)
prod(b)
prod(c)
a[2]^3
b[2]^3
c[2]^3

# concave / flat 

a <- 2^c(0,3,3)
b <- 2^c(2,2,2)
c <- apply(rbind(a,b),2,mean)
matplot(cbind(a,b,c),pch=1)
matplot(log(cbind(a,b,c)),pch=1)
prod(a)
prod(b)
prod(c)
a[2]^3
b[2]^3
c[2]^3

# two concave

a <- 2^(c(1.5,3,1.5))
b <- 2^(c(0,6,0))
c <- apply(rbind(a,b),2,mean)
matplot(cbind(a,b,c),pch=1)
matplot(log(cbind(a,b,c)),pch=1)
prod(a)
prod(b)
prod(c)
a[2]^3
b[2]^3
c[2]^3

# two concave 2

a <- 2^(c(3,3,0))
b <- 2^(c(0,3,3))
c <- apply(rbind(a,b),2,mean)
matplot(cbind(a,b,c),pch=1)
matplot(log(cbind(a,b,c)),pch=1)
prod(a)
prod(b)
prod(c)
a[2]^3
b[2]^3
c[2]^3


