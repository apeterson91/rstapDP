## code to prepare `ClusterScales_two.R` dataset goes here

set.seed(3431)
num_subj <- 5E2
Z <- rbinom(num_subj,1,.5)
high_risk <- c(rep(0,floor(num_subj/2)),rep(1,num_subj-floor(num_subj/2)))
cnt <- rpois(num_subj,10)
ldists <- sapply(cnt,function(x) rgamma(x,shape = 5, rate = 5) )
f <- function(x) (1 -.1*x-.2*x^2 -.1*x^3)*(x<=1.545)
exposure <- cbind(cnt,sapply(ldists,sum),sapply(ldists,function(x) sum(f(x))))
simple <- rexp(num_subj)

y <- 26 + Z * -2.2 +  high_risk*((simple)/max(simple))*10  + (1-high_risk)*(simple)/max(simple)*-10 + rnorm(num_subj)



ClusterScales_two <- dplyr::tibble(id=1:num_subj,
                                   BMI = y,
                                   sex = Z,
                                   simple = simple)

ClusterScales_two_ddata <- purrr::map2_dfr(1:length(ldists),ldists,function(x,y) dplyr::tibble(id=x,Distance=y))

usethis::use_data(ClusterScales_two, overwrite = TRUE)
usethis::use_data(ClusterScales_two_ddata, overwrite = TRUE)

