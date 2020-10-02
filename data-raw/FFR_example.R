## code to prepare `FFR_subjects` and `FFR_distances` datasets

set.seed(3431)
num_subj <- 1E3
Z <- rbinom(num_subj,1,.5)
high_risk <- c(rep(0,floor(num_subj/2)),rep(1,num_subj-floor(num_subj/2)))
has_exp <- rbinom(num_subj,size = 1,prob = .95)
cnt <- rpois(num_subj,10)*has_exp
ldists <- lapply(cnt,function(x) runif(x) )
f <- function(x) 2*pweibull(x,shape=5,scale=.5,lower.tail = F)
exposure <- sapply(ldists,function(x) sum(f(x)))


y <- 26 + Z * -2.2 +  high_risk*((exposure))  + (1-high_risk)*(exposure)*0 + rnorm(num_subj)



FFR_subjects <- dplyr::tibble(id=1:num_subj,
                              BMI = y,
                              sex = Z)

FFR_distances <- purrr::map2_dfr(1:length(ldists),ldists,function(x,y) dplyr::tibble(id=x,Distance=y))

FFR_benvo <- rbenvo::benvo(subject_data = FFR_subjects,sub_bef_data = list(FFR=FFR_distances))
usethis::use_data(FFR_benvo,overwrite=TRUE)
