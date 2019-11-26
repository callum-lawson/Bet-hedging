### Population growth curve analysis of the storage effect using Grafen's reproductive values ### 

n_eps <- 100 
seq_eps <- seq(-1,1,length.out=100) # climate range

seq_a <- c(1,-1) # climate sensitivity
S <- 0 # adult survival

f_r <- function(S,a,eps){
  S + exp(a * eps)
  # S + a * eps + 1
} 
  # N Births = a * eps

r_A <- f_r(S,seq_a[1],seq_eps)
r_B <- f_r(S,seq_a[2],seq_eps)

r_A[r_A<0] <- NA
r_B[r_B<0] <- NA

w_A_inv <- r_A/r_B
w_B_inv <- r_B/r_A

par(mfrow=c(2,2))

matplot(seq_eps,cbind(w_A_inv,w_B_inv),type="l")
abline(h=1,lty=3,col="gray")
matplot(seq_eps,log(cbind(w_A_inv,w_B_inv)),type="l")

# Is an environmentally-insensitive type uninvadable? ---------------------

r_C <- f_r(S,a=0,seq_eps)
w_A_invC <- r_A/r_C
w_B_invC <- r_B/r_C
matplot(seq_eps,cbind(w_A_invC,w_B_invC),type="l")
abline(h=1,lty=3,col="gray")
  # not with exponential reproduction - C can still be invaded
  # but C seems to be *un*invadable with artificial "survival + linear reproduction" model
  # (even though the extreme types in that model have very convex invasion functions of each other)
  # (Remember that to calculate whether the strategy can be invaded, we're looking at an average over
  # *several* environments (here centered at 1), rather than whether the resident can be beaten in a 
  # *given* environment)

# Is mixed A-B strategy uninvadable? --------------------------------------

r_D <- apply(cbind(r_A,r_B),1,mean)
w_A_invD <- r_A/r_D
w_B_invD <- r_B/r_D
matplot(seq_eps,cbind(w_A_invD,w_B_invD),type="l")
abline(h=1,lty=3,col="gray")
  # appears to be uninvadable, at least for A and B
  # raises question again: why doesn't everything bet-hedge?
  # also suggests that Chesson's storage effect coexistence can regularly by beaten by bet-hedging
