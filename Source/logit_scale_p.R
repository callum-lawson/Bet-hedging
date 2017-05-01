###########################################################################
# How does stochastic growth rate relate to p when p is on a logit scale? #
###########################################################################

q <- c(1/3,2/3)
lambda1 <- c(1/2,2)
lambda2 <- c(2,1/2)
p <- seq(0,1,length.out=100)
r <- q[1]*log(p*lambda1[1]+(1-p)*lambda2[1]) + (1-q[2])*log(p*lambda1[2]+(1-p)*lambda2[2])

par(mar=c(4,4,1,1))
plot(r~p)
plot(r~qlogis(p),type="l")
  # r non-linear in alpha