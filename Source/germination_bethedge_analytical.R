#####################################################################################
# Derivatives of lambda and r to calculate effects of mean and variance in rainfall #
#####################################################################################

l_z <- function(z,G,S){
  G*exp(z) + (1-G)*S  
  }

# assuming that Y is a linear exponential function of z

dl_dz <- d2l_dz2 <- function(z,G,S){
  G*exp(z)
  }

r_z <- function(z,G,S){
  log(l_z(z,G,S))
  }

dr_dz <- function(z,G,S){
  G*exp(z) / ( G*exp(z) + S*(1-G) )
  }

d2r_dz2 <- function(z,G,S){
  G*(1-G)*S*exp(z) / ( G*exp(z) + S*(1-G) )^2
  }

# Env response

G <- 0.5
S <- 0.8
zlim <- c(-5,5)

par(mfrow=c(1,1))
curve(l_z(x,G,S),xlim=zlim)
curve(dl_dz(x,G,S),xlim=zlim)
curve(d2l_dz2(x,G,S),xlim=zlim)

curve(r_z(x,G,S),xlim=zlim)
curve(dr_dz(x,G,S),xlim=zlim)
curve(d2r_dz2(x,G,S),xlim=zlim)
  # logistic?

# Effect of G
  
z <- 2
Glim <- c(0,1)

par(mfrow=c(1,1))
curve(l_z(z,G=x,S),xlim=Glim)
  # z=0 -> G increases mean lambda
curve(dl_dz(z,G=x,S),xlim=Glim)
  # z=0 -> G increases variability in lambda
curve(d2l_dz2(z,G=x,S),xlim=Glim)
  # z=0 -> G increases positive effects of variability on mean lambda
  # (so G doubly positive for mean lambda in variable environment)

curve(r_z(z,G=x,S),xlim=Glim)
curve(dr_dz(z,G=x,S),xlim=Glim)
  # higher G shows stronger response to change in mean env
  # and year-to-year fluctuations larger
curve(d2r_dz2(z,G=x,S),xlim=Glim)
  # intermediate G shows strongest response to variance


