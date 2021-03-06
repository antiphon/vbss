### test logistic fitting when only SS around
source("../R/logistic.R")
N <- 10
P <- 5
P0 <- 2
set.seed(1234)
## covariates: SS
X <- matrix(rnorm((P0+P)*N, 0, 2), nrow=N)
s <- c(rep(1,P), rep(0, P0))
w <- s*rnorm((P+P0))
## response:
logit <- function(x)1/(1+exp(-x))
#
y <- rbinom( N, 1, logit(c(X%*%w) ))
## observed covariates:
colnames(X) <- paste0("X.", 1:(P+P0))
##
#res<- vbglmss(y=y, X=X, family=binomial, verbose=TRUE, prior=list(pi=0.01))
res<- vbglmss.fit.logistic.X(y=y, X=X, verbose=TRUE, prior=list(pi=0.01))

plot(1:(P0+P), res$w$pi, col=s+2)