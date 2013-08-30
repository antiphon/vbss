## test the covariance approximation
library(vbglmss)
set.seed(123456)
K <- 500
N <- 425
source("../R/mixed_gaussian_logistic.R")
source("../R/sparseSolve.R")



X <- Matrix(matrix(rpois(K*N, 0.4), nrow=N))
w <- rnorm(K, 0, 4)
logit <- function(x) 1/(1+exp(-x))
y <- rbinom(N, 1, as.numeric(logit(X%*%w)))
y2 <- rnorm(N, as.numeric(X%*%w), 1)



#if(!exists('res2')){
{  T2 <- Sys.time()
  res2<- vbglmss.mixedSS(h=y2, Zh=X, max.iter=100, prior=list(beta=list(diagonal=F), 
                                               gaussian=list(tau=c(100,100))),
                       verbose=T, eps=1e-5)
  T2 <- format(Sys.time()-T2)
}
T1 <- Sys.time()
res1<- vbglmss.mixedSS(h=y2, Zh=X, 
                       prior=list(#beta=list(diagonal=F, Si=res2$beta$Si, S=res2$beta$S, off_from_prior=TRUE),
                         beta=list(diagonal=!F, S=Diagonal(n=K)*0.0001),
                                  gaussian=list(tau=c(100,100))
                                  ),
                                          verbose=T, eps=1e-5)
T1 <- format(Sys.time()-T1)

print(data.frame(T1,T2))

plot(w, res1$beta$m, ylim=range(res1$beta$m, res2$beta$m))
points(w, res2$beta$m, col=2, pch=19, cex=.4)
abline(0,1)
