## code to prepare `ClusterScales_one` dataset goes here

num_subj <- 2E2
Z <- rbinom(num_subj,1,.5)
high_risk <- 1:floor((num_subj/2))
hrmat <- diag(num_subj)
hrmat[high_risk,] <- rep(0,num_subj)
dists <- rexp(n = num_subj, rate = .5)
y <- 26 + Z * -2.2 + hrmat %*% dists * -3 + (diag(num_subj)-hrmat) %*% dists *2 + rnorm(num_subj)



ClusterScales_one <- dplyr::tibble(BMI = y,
                            sex = Z,
                            dists = dists)

usethis::use_data(ClusterScales_one, overwrite = TRUE)

# half_way <- floor(num_subj/2)+1
# r <-  y - 26 - Z * -2.2
# denom <- sum((hrmat %*% dists) * dists)
# S <- (hrmat %*% dists) %*% t((hrmat %*% dists)) / denom
# bhat <- (S %*% r)[half_way]/dists[half_way]
# bhat
# 
# lrmat <- diag(num_subj) - hrmat
# r <-  y - 26 - Z * -2.2
# denom <- sum((lrmat %*% dists) * dists)
# S <- (lrmat %*% dists) %*% t((lrmat %*% dists)) / denom
# bhat <- (S %*% r)[1]/dists[1]
# bhat
# 
# 
# ####### 
# 
# head(dnorm(y - 26 - Z * -2.2,mean = dists *.5,sd = 1,log = T))
# 
# head(dnorm(y - 26 - Z * -2.2,mean = dists * 1.5,sd = 1,log = T))
# 
# ib <- rnorm(5)
# pbs <- sapply(ib,function(x) dnorm(y-26-Z *-2.2, mean = dists*x, sd = 1, log = T))
# 
# labs <- apply(pbs,1,function(x) sample(1:5,size = 1,replace=F,prob = exp(x-max(x)) ))
# 
# r <- y - 26 - Z * -2.2 
# mat <- diag((labs==4)*1)
# denom <- sum( (mat %*% dists) * dists)
# S <- (mat %*% dists) %*% t(mat %*% dists) / denom
# bhat <- (S %*% r)[2] / dists[2]
