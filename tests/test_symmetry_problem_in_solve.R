## debug symmetry bug
## test the covariance approximation
library(vbglmss)
source("../R/mixed_gaussian_logistic.R")
source("../R/sparseSolve.R")

load("trouble_for_vbglmms.rda")

res1<- with(dat, vbglmss.mixedSS(h=h, y=y, Zy=Zy, Zh=Zh,
                       prior=prior,
                       verbose=T))
