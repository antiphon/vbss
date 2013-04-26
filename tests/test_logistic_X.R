### test logistic fitting when only SS around
library(vbglmss)
library(devtools)
reload("../")
#source("../R/logistic.R")
#source("../R/vbglmss.R")
N <- 100
P1 <- 50
P0 <- 20
K <- P1+P0
set.seed(123456)
## covariates: SS
X <- matrix(rnorm(K*N, 0, 2), nrow=N)
s <- c(rep(1, P1), rep(0, P0))
w <- s*rnorm(K)
## response:
logit <- function(x)1/(1+exp(-x))
#
y <- rbinom( N, 1, logit(c(X%*%w) ))
## observed covariates:
colnames(X) <- paste0("X.", 1:K)
##
#res<- vbglmss(y=y, X=X, family=binomial, 
#              prior=list(pi=0.05, tau=list(a=10,b=10)), verbose=T)
#plot(1:K, res$SS$g, col=s+2, ylim=0:1, pch=19)

## sparse
X2 <- Matrix(matrix(rpois(K*N, 0.4), nrow=N))
w2 <- s*rnorm(K, 0, 4)

y2 <- rbinom( N, 1, logit(as.numeric(X2%*%w2) ))
## observed covariates:
colnames(X2) <- paste0("X.", 1:K)

res2<- vbglmss(y=y2, X=X2, family=binomial, 
              prior=list(pi=0.5, tau=list(a=100,b=100)), verbose=T, eps=1e-5)
plot(res2$SS$g, col=s+2, ylim=c(0,1))


