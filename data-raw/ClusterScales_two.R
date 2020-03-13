## code to prepare `ClusterScales_two.R` dataset goes here

num_subj <- 5E2
Z <- rbinom(num_subj,1,.5)
high_risk <- 1:floor((num_subj/2))
hrmat <- diag(num_subj)
hrmat[high_risk,] <- rep(0,num_subj)
dists <- rgamma(n = num_subj,shape = 15,rate = 15)
dists_sq <- dists^2
dists_mat <- cbind(rpois(num_subj,lambda = 2),dists,dists_sq)
betas_one <- c(2,2,-1)
betas_two <- c(0,0,0)
y <- 26 + Z * -2.2 + hrmat %*% dists_mat %*% betas_one  + (diag(num_subj)-hrmat) %*% dists_mat %*%  betas_two + rnorm(num_subj)

ClusterScales_two <- dplyr::tibble(BMI = y,
                                   sex = Z,
                                   dists_mat = dists_mat)

usethis::use_data(ClusterScales_two, overwrite = TRUE)

