#############################################################################
# Figures to help understand reproductive values perspective on bet-hedging #
# (Grafen 1999 Proc Roy Soc)                                                #
#############################################################################

# maindir <- "D:/Users/calluml/Dropbox/NIOO/"
# maindir <- "C:/Users/Callum/Documents/My Dropbox/NIOO/"

setwd(paste0(maindir,"Analyses/Bet-hedging"))

### FUNCTIONS

rrat <- function(p,r1,r2,logis=F){
	ratio <- r1/(p*r1+(1-p)*r2)
	}

invas <- function(p,r1,r2){
	sapply(p,function(P){
		mean(c(rrat(P,r1[1],r2[1]),rrat(P,r1[2],r2[2])))
		})
	}

invasplot <- function(r1,r2){
	curve(invas(p=x,r1,r2))
 	abline(h=1,col="red")
	}

invassim <- function(r1,r2,
	p1start=0.01,p2start=0.99,nyear=500,random=T,regulated=T){
	
	if(random==T) year <- sample(c(1,2),nyear,replace=T) # Random
	if(random==F) year <- rep(c(1,2),times=nyear/2) # Systematic
	lam1 <- r1[year]
	lam2 <- r2[year]
	p1 <- p2 <- vector()
	p1[1] <- p1start
	p2[1] <- p2start

	for(i in 2:nyear){
		if(regulated==T){
			p1[i] <- p1[i-1]*lam1[i]/( p1[i-1]*lam1[i] + p2[i-1]*lam2[i] )
			p2[i] <- p2[i-1]*lam2[i]/( p1[i-1]*lam1[i] + p2[i-1]*lam2[i] )
			}
		if(regulated==F){
			p1[i] <- p1[1]*prod(lam1[2:i])/( p1[1]*prod(lam1[2:i])+p2[1]*prod(lam2[2:i]) )
			p2[i] <- p2[1]*prod(lam2[2:i])/( p1[1]*prod(lam1[2:i])+p2[1]*prod(lam2[2:i]) )
			}
		# multiply up and calculate frequency at the end
		}

	mydata <- data.frame(year=year,lam1=lam1,lam2=lam2,p1=p1,p2=p2)
	return(mydata)
	
	}

### DIFFERENT FREQUENCIES 

r1 <- c(1,0.9)
r2 <- c(2,0.5)
cor.test(rep(r1,2),rep(r2,2))

invasplot(r1,r2)

### INFERIOR STRATEGY

r3 <- c(1.75,0.5)
invasplot(r3,r1)

v31 <- invassim(r3,r1)
plot(v31$p1)

### ANTICOR
delta1 <- 0.2
delta2 <- 0.9
r4 <- c(1-delta1,1+delta1)
r5 <- c(1.25-delta2,1.25+delta2)
invasplot(r4,r5)

set.seed(1)
v45a <- invassim(r4,r5,regulated=T)
set.seed(1)
v45b <- invassim(r4,r5,regulated=F)
par(mfrow=c(1,2))
plot(v45a$p1,type="l")
plot(v45b$p1,type="l",col="red")

### GEOM REPRODUCTIVE VALUE

sqrt(1.5/1 * 0.5/2)
sqrt(1.5/1) * sqrt(0.5/2)
	# SAME

### MIXED STRATEGIES

r1 <- c(1,1)
r2 <- c(1.5,0.5)

pi <- 0.001
mix <- c(pi*r1[1] + (1-pi)*r2[1], pi*r1[2] + (1-pi)*r2[2])

invas(p=0,r1=mix,r2=r2)
	# mixed strategies can't avoid invasion by pure with
	# smaller geometric mean

u=v=w=x=y=z=c(0,1)
envmat <- expand.grid(u,v,w,x,y,z) + 1
endrat <- apply(envmat,1,function(x){
	prod(r1[x])/prod(r2[x])
	})

hist(endrat,breaks=20)

###############################
### LOGISTIC TRANSFORMATION ###
###############################

mu1 <- 1
mu2 <- 1
sig1 <- 1
sig2 <- 2

nt <- 100
N0 <- 1

rseq1 <- rnorm()

### 
### OTHER
###

np <- 10^3

# a <- c(4,1/2)
# b <- c(2,2)
# -> p = 1/6

a <- c(2,10)
b <- c(4,5)

a <- c(2,1/2)
b <- c(1,1)
p <- seq(0,1,length.out=np)
q1 <- 1/2
q <- c(q1,1-q1)
mix <- function(a,b,p){
  p*a + (1-p)*b
  }

m <- matrix(nr=np,nc=2)
for(i in 1:np){
  m[i,] <- mix(a,b,p[i])
  }

rv <- function(a,b){
  sum(q*(a/b))
  }
r <- r1 <- r2 <- matrix(nr=np,nc=np)
for(i in 1:np){
  for(j in 1:np){
    r[i,j] <- rv(m[i,],m[j,])
    r1[i,j] <- rv(m[i,1],m[j,1])
    r2[i,j] <- rv(m[i,2],m[j,2])
    }
  }
library(fields)
image.plot(p,p,log(r))
image.plot(p,p,r)

sd <- apply(r,2,sd)
plot(sd~p)
opt <- which(sd==min(sd))
p[opt]
prod(a^q); prod(b^q)
m[opt,] # doesn't necessarily minimise variance in absolute / log fitnesses
  # but isn't "betting" on one env more than the other - just gets higher payoff
plot(r[,opt]~p)

maxr <- apply(r,2,max)
opt2 <- which(maxr==min(maxr))

plot(r[,opt2]~p)

matplot(p,cbind(r1[500,],r2[500,]),type="l")


yo <- function(mean,sd){
  (exp(sd^2)-1)*exp(2*mean+sd^2)
  }

curve(yo(mean=1,sd=x),xlim=c(0,1))
curve(yo(mean=2,sd=x),col="red",add=T)

invader <- function(a,b){
  (a[1]/b[1]+a[2]/b[2])/2
  }

res <- c(3/2,3/4)
invader(c(1,1),res)
invader(c(2,1/2),res)
three <- c((1/3*2)+(2/3)*1,(1/3)*(1/2)+(2/3)*1)
invader(three,res)

### ELLNER FUNCTION

lamnorm <- function(z,y){
  exp(-(z-y)^2/2)
  }

ellpure <- function(z,y){
  -(z-y)^2/2
  }

ellmix <- function(z,y1,y2){
  l1 <- lamnorm(z,y1)
  l2 <- lamnorm(z,y2)
  l <- apply(cbind(l1,l2),1,mean)
  log(l)
  }

ellrange <- function(z,y1,y2){
  yseq <- seq(y1,y2,length.out=100)
  lm <- outer(z,yseq,lamnorm)
  l1 <- exp(-(z-y1)^2/2)
  l2 <- exp(-(z-y2)^2/2)
  l <- apply(lm,1,mean)
  log(l)
  }
  # assuming uniform mix

pdf(paste0("normal_mixing_",format(Sys.Date(),"%d%b%Y"),".pdf"),width=4,height=4)

par(mfrow=c(1,1),mar=c(4,4,2,2))

curve(ellpure(x,0),xlim=c(-5,5),ylab="r",xlab="z")
# same as squared function tried earlier
curve(ellmix(x,-0.5,0.5),add=T,col="red")
curve(ellrange(x,-0.5,0.5),add=T,col="blue")

curve(ellpure(x,0),xlim=c(-5,5),ylab="r",xlab="z")
# same as squared function tried earlier
curve(ellmix(x,-1,1),add=T,col="red")
curve(ellrange(x,-1,1),add=T,col="blue")

curve(ellpure(x,0),xlim=c(-5,5),ylab="r",xlab="z")
# same as squared function tried earlier
curve(ellmix(x,-2,2),add=T,col="red")
curve(ellrange(x,-2,2),add=T,col="blue")

# Comp

curve(ellpure(x,0),xlim=c(-2.5,2.5),ylab="r",xlab="z")
curve(ellmix(x,-0.25,0.25),add=T,col="red")
curve(ellmix(x,-0.5,0.5),add=T,col="orange")
curve(ellmix(x,-1,1),add=T,col="blue")
curve(ellmix(x,-2,2),add=T,col="purple")

curve(ellpure(x,0),xlim=c(-2.5,2.5),ylab="r",xlab="z")
curve(ellrange(x,-0.5,0.5),add=T,col="red")
curve(ellrange(x,-1,1),add=T,col="orange")
curve(ellrange(x,-2,2),add=T,col="blue")
curve(ellrange(x,-3,3),add=T,col="purple")

dev.off()

### Cohen two-phenotype example

a <- c(1/4,4)
b <- c(1/2,1/2)
q1 <- 0.9

a <- c(1/2,2)
b <- rev(a)
q1 <- 0.49

q <- c(q1,1-q1)
cohen <- function(x) sapply(x,function(p) exp(sum(q*log(p*a+(1-p)*b))))
curve(cohen(x),xlim=c(0,1))

# Cohen Harmonic means

a <- c(1/4,4)
b <- c(1/2,1/2)
q1 <- 0.9
q <- c(q1,1-q1)

hmean <- function(x){
  1/sum(q*(1/x))
  }

ginv <- function(x,y){
  sum(q*(x/y))
  }

hmean(a)
hmean(b)

ginv(a,b)
ginv(b,a)

### Grafen log scale perspective

a <- c(1/2,2)
b <- c(4,1/4)
q1 <- 0.5
q <- c(q1,1-q1)

ginv <- function(x,y){
  sum(q*(x/y))
  }

glinv <- function(x,y){
  sum(log(q)+log(x)-log(y)) # doesn't work because 
  }

glinv(a,b) # compares how good they are in absolute terms
glinv(b,a)

log(ginv(a,b))
log(ginv(b,a))

### "Swap costs"

par(mfrow=c(2,1),mar=c(4,4,2,2))

curve(exp(5*x)/100,xlim=c(0,1))
abline(h=0.5,lty=2)

curve(5*x-log(100),xlim=c(0,1))
abline(h=log(0.5),lty=2)




