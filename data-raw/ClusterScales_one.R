## code to prepare `ClusterScales_one` dataset goes here

library(tidyverse)


Z <- rbinom(100,1,.5)
y <- 26 + Z * -2.2 + rnorm(100)


ClusterScales_one <- tibble(BMI = y,
                            sex = Z)

usethis::use_data(ClusterScales_one, overwrite = TRUE)
