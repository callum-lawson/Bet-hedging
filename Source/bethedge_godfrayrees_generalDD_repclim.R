#################################################################################
# Simulations to see how underlying mean and variance in germination thresholds #
# translates to population growth rates and density under variable rainfall     #
#################################################################################

# Load packages -----------------------------------------------------------

library(parallel)

# Define functions --------------------------------------------------------

logistic <- function(x,alpha_G,beta_Gz){
  plogis(alpha_G + beta_Gz*x)
  }

betpopsim <- function(nt,zt,X0,alpha_Y,beta_Yz,alpha_m,
  alpha_G,beta_Gz,eps_Y,FUN){
  Xt <- vector("numeric",length=nt) 
  G <- plogis(alpha_G + beta_Gz*zt) # density-independent
  Xt[1] <- X0
  for(t in 2:nt){
    Xg <- log(G[t]) + Xt[t-1] 
    Xo <- log(1-G[t]) - alpha_m  + Xt[t-1]  # S = exp(-alpha_m)
    Xn <- Xg + alpha_Y + beta_Yz*zt[t] + eps_Y[t]
    Xndash <- FUN(Xn,Xg)
    Xt[t] <- log(exp(Xo) + exp(Xndash))
    }
  return(Xt)
  }

fullpopsim <- function(nt,ni,np,zt,X0,alpha_Y,beta_Yz,alpha_m,
  alpha_G,beta_Gz,eps_Y,FUN,multi){
  Xt <- array(NA,dim=c(nt,ni,np,np))
  if(multi==T){
    for(i in 1:ni){
      for(p in 1:np){
        for(q in 1:np){
          Xt[,i,p,q] <- betpopsim(nt,zt[,i],X0,alpha_Y,beta_Yz,alpha_m,
            alpha_G=alpha_G[p,q],beta_Gz=beta_Gz[p,q],eps_Y,FUN)
          }
        }
      }
    }
  if(multi==F){
    for(p in 1:np){
      for(q in 1:np){
        Xt[,,p,q] <- betpopsim(nt,zt,X0,alpha_Y,beta_Yz,alpha_m,
          alpha_G=alpha_G[p,q],beta_Gz=beta_Gz[p,q],eps_Y,FUN)
        }
      }
    }
  return(Xt)
  }

DDidentity <- function(Xn,Xg){
  Xn
  }

DDceiling <- function(Xn,Xg){
  Xnmax <- 5
  ifelse(Xn>Xnmax,Xnmax,Xn)
  }

DDgompertz <- function(Xn,Xg){
  beta_Gd <- -0.5
  Xn + beta_Gd*Xg
  }

DDricker <- function(Xn,Xg){
  beta_Gd <- -0.001
  Xn + beta_Gd*exp(Xg)
  }

medfun <- function(x,FUN,...){
  z <- apply(x,c(2,3,4),FUN,...)
  a <- apply(z,c(2,3),median)
    # dims changed: 3,4 -> 2,3
  a[a==-Inf] <- NA
  return(a)
  }

medgrow <- function(x) medfun(x,FUN=function(y) (y[nt]-y[1])/nt )
medsize <- function(x) medfun(x,FUN=median)

gillgrow <- function(y,output){
  lambda <- exp(diff(y)) # y = X = already logged
  mu <- mean(lambda)
  if(output=="mu") return(mu)
  if(output %in% c("dis","gill")){
    lamvar <- var(lambda) 
    dis <- lamvar/(2*mu) # discount
    }
  if(output=="dis") return(dis)
  if(output=="gill") return(mu-dis)
  } # Gillespie 1973 Genet Res

medmu <- function(x) medfun(x,FUN=gillgrow,output="mu")
meddis <- function(x) medfun(x,FUN=gillgrow,output="dis")
medgill <- function(x) medfun(x,FUN=gillgrow,output="gill")

persist <- function(x,Xmin=0){
  z <- apply(x,c(2,3,4),function(y) !T %in% (y<Xmin))
  apply(z,c(2,3),function(x) sum(x)/length(x)) 
    # dims changed: 3,4 -> 2,3
  }

germimage <- function(z,mulim=0:np,sdlim=0:np){
  require(fields)
  image.plot(tau_mu_seq[mulim],tau_sd_seq[sdlim],z[mulim,sdlim])
  }

optiplot <- function(x){
  opti_pos <- which(x==max(x,na.rm=T))
  if(length(opti_pos)>1){
    opti_pos <- sample(opti_pos,1)
    warning("multiple optima")
    }
  curve(logistic(x,alpha_G[opti_pos],beta_Gz[opti_pos]),xlim=c(-3,3))
  list(alpha_G=alpha_G[opti_pos],beta_Gz=beta_Gz[opti_pos])
  }

# Choose parameters -------------------------------------------------------

nt <- 100
ni <- 100
np <- 50
# 20-21 mins

zmu <- 0
zsd <- 1

tau_mu_min <- -5
tau_mu_max <- 5
tau_sd_min <- 0.01
tau_sd_max <- 5

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

# Plot parameter effects --------------------------------------------------

zseq <- seq(zmu-2*zsd,zmu+2*zsd,length.out=ni)
G <- array(NA,dim=c(ni,np,np))
for(i in 1:np){
  for(j in 1:np){
    G[,i,j] <- plogis(alpha_G[i,j] + beta_Gz[i,j]*zseq)
    }
  }

library(fields)
par(mfrow=c(1,1))
matplot(zseq,G[,,25],type="l",col=tim.colors(np),lty=1,main=expression(tau[mu]~change), ylab="G")
matplot(zseq,G[,25,],type="l",col=tim.colors(np),lty=1,main=expression(tau[sigma]~change), ylab="G")
matplot(zseq,qlogis(G[,,25]),type="l",col=tim.colors(np),lty=1,main=expression(tau[sigma]~change), ylab="inv_logit(G)")
matplot(zseq,qlogis(G[,25,5:25]),type="l",col=tim.colors(np),lty=1,main=expression(tau[mu]~change), ylab="inv_logit(G)")
  # change in mean germ threshold -> separate effects of DD
  # change in var in germ threshold -> interactive effects of DD on z

# Run sims ----------------------------------------------------------------

zt_vals <- list(c=zc,v=zt_mat,v=zt_mat)
eps_Y_vals <- list(c=eps_Yc,c=eps_Yc,v=eps_Yv)
DD_FUN_vals <- list(i=DDidentity,c=DDceiling,g=DDgompertz,r=DDricker)

ns <- length(zt_vals)     # number of scenarios
nd <- length(DD_FUN_vals) # number of density values

zt_l <- rep(zt_vals,times=nd)
eps_Y_l <- rep(eps_Y_vals,times=nd)
DD_FUN_l <- rep(DD_FUN_vals,each=ns)
multi_l <- as.list(sapply(zt_l,is.matrix)) # zt must match multi

nc <- length(zt_l) # number of cores

system.time({
CL <- makeCluster(nc)
clusterExport(cl=CL,c(
  "fullpopsim","betpopsim",
  "nc","nt","ni","np",
  "X0","alpha_Y","beta_Yz","alpha_m",
  "zt_l","alpha_G","beta_Gz",
  "eps_Y_l","DD_FUN_l","multi_l"
  ))
Xt_l <- parLapply(CL, 1:nc, function(n){
  fullpopsim(nt=nt,ni=ni,np=np,zt=zt_l[[n]],X0=X0,
    alpha_Y=alpha_Y,beta_Yz=beta_Yz,alpha_m=alpha_m,
    alpha_G=alpha_G,beta_Gz=beta_Gz,eps_Y=eps_Y_l[[n]],
    FUN=DD_FUN_l[[n]],multi=multi_l[[n]])
  })
stopCluster(CL)
})

names(Xt_l) <- paste0("z",names(zt_l),"_e",names(eps_Y_l),"_d",names(DD_FUN_l))[1:nc]

# Calculate statistics ----------------------------------------------------

Xbar_l <- lapply(Xt_l,medsize)
per_l <- lapply(Xt_l,persist)

Xbar_diff_l <- mapply('-',
  Xbar_l[grep("zv_ev",names(Xt_l))],
  Xbar_l[grep("zc",names(Xt_l))],
  SIMPLIFY=F
  )

per_fac_l <- mapply('/',
  per_l[grep("zv_ev",names(Xt_l))],
  per_l[grep("zc",names(Xt_l))],
  SIMPLIFY=F
  )
per_fac_l <- lapply(per_fac_l,function(x){
  x[!is.finite(x)] <- NA
  x
  })

# Constant environment ----------------------------------------------------

germimage(Xbar_l$zc_ec_di)
germimage(Xbar_l$zc_ec_dc)
germimage(Xbar_l$zc_ec_dg)
germimage(Xbar_l$zc_ec_dr)  
germimage(Xbar_l$zc_ec_di,mulim=1:25)  # optimum = full germination
germimage(Xbar_l$zc_ec_dc,mulim=25:50)  # optimum = approx 50% germination?
germimage(Xbar_l$zc_ec_dg,mulim=1:25)  # optimum = full or partial germination
germimage(Xbar_l$zc_ec_dr,mulim=1:25)  # optimum = partial germination

optiplot(Xbar_rick_c)
matplot(1:nt,Xt_l$zc_ec_dr[1:100,1,25,],type="l",col=tim.colors(np),lty=1)
  # With DD, higher popsize achieved by having intermediate G value,
  # i.e. intermediate mean threshold, higher threhsold variance
  # happens because at max Y anyway, so some non-germination helps to save seeds in soil

germimage(per_l$zc_ec_di,mulim=25:50)
germimage(per_l$zc_ec_dc,mulim=25:50)
germimage(per_l$zc_ec_dg,mulim=25:50)
germimage(per_l$zc_ec_dr,mulim=25:50) 
  # there's a minimum optimal germination fraction that leads to extinction
  # within those populations that persist, some grew faster than others
  # too high germination -> smaller popsize, but too low germination -> extinction

optiplot(per_l$zc_ec_di)

# Reliable variability ----------------------------------------------------

germimage(Xbar_l$zv_ec_di)
germimage(Xbar_l$zv_ec_dc)
germimage(Xbar_l$zv_ec_dg)
germimage(Xbar_l$zv_ec_dr) 
  # DI -> no bet-hedging; DD -> bet-hedging helps (for all, but very little for ceiling)

optiplot(Xbar_l$zv_ec_dc)
matplot(1:nt,Xt_l$zv_ec_dg[,1,1,],type="l",col=tim.colors(np),lty=1)

rbar_iden <- medgrow(Xt_l$zv_ec_di)
germimage(rbar_iden)

mubar_iden <- medmu(Xt_l$zv_ec_di)
disbar_iden <- meddis(Xt_l$zv_ec_di)
gillbar_iden <- medgill(Xt_l$zv_ec_di)

germimage(mubar_iden)
germimage(disbar_iden)
germimage(gillbar_iden)
  # Gillespie approximation can be really bad

germimage(per_l$zv_ec_di)
germimage(per_l$zv_ec_dc)
germimage(per_l$zv_ec_dg)
germimage(per_l$zv_ec_dr) 
  # iden -> risk wins
  # ceil -> risk wins (surprising?)
  # gomp and rick -> bet-hedging wins
  # so with some DD forms, bet-hedging can help even if signal is reliable

optiplot(per_l$zv_ec_dg)

# Unreliable variability --------------------------------------------------

germimage(Xbar_l$zv_ev_di)
germimage(Xbar_l$zv_ev_dc)
germimage(Xbar_l$zv_ev_dg)
germimage(Xbar_l$zv_ev_dr)
  # bet-hedging good for everyone

matplot(1:nt,Xt_l$zv_ev_dg[,1,,1],type="l",col=tim.colors(np),lty=1)

germimage(per_l$zv_ev_di)
germimage(per_l$zv_ev_dc)
germimage(per_l$zv_ev_dg)
germimage(per_l$zv_ev_dr)
  # unreliable signal -> bet-hedgers always have lower extinction risk 
  # (except for DI, but perhaps noise not large enough relative to env var?)

optiplot(per_l$zv_ev_di)

# Env var increase --------------------------------------------------------

germimage(Xbar_diff_l$zv_ev_di)
germimage(Xbar_diff_l$zv_ev_dc)
germimage(Xbar_diff_l$zv_ev_dg)
germimage(Xbar_diff_l$zv_ev_dr)
  # those that benefit most are bet-hedgers, but not max bet-hedgers
germimage(Xbar_diff_l$zv_ev_dr,mulim=40:50,sdlim=10:50)

germimage(per_fac_l$zv_ev_di)
germimage(per_fac_l$zv_ev_dc)
germimage(per_fac_l$zv_ev_dg)
germimage(per_fac_l$zv_ev_dr)
  # those that benefit most are bet-hedgers - possibly max bet-hedgers (?)
  # white space -> infinite persistence increase
