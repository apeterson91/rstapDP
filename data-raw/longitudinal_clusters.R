## code to prepare `longitudinal_clusters` dataset goes here

set.seed(3431)
num_subj <- 5E2
riskcat <- sapply(1:num_subj,function(x) sample(1:3,size=1,replace=F,prob = c(.7,.2,.1)))
riskmat <- Matrix::fac2sparse(riskcat)
num_visits <- sample(1:5,size=num_subj,replace=T)
visit_num <- purrr::map_dfr(1:num_subj,function(x) dplyr::tibble(id = x,measurement = 1:num_visits[x]))
num_obs <- sum(num_visits)
Z <- rbinom(num_subj,1,.5)

sjdf <- dplyr::tibble(id=1:num_subj,
                      sex = Z,
                      riskcat = riskcat,
                      subj_effect = rnorm(num_subj,mean = 0,sd = .4),
                      subj_slope = rnorm(num_subj,mean=0,sd=0.3))

has_exp <- rbinom(num_obs,size = 1,prob = .95)
cnt <- rpois(num_obs,10)*has_exp
ldists <- lapply(cnt,function(x) runif(x,0,2) )
majority <- function(x) .5*pweibull(x,shape=5,scale=.6,lower.tail = F) # formerly .3 *
## not going to matter
low_risk <- function(x) return(0)
high_risk <- function(x) 2*pweibull(x,shape=3,scale=1,lower.tail = F)

FFR_effect <- function(riskcat,distance){
  if(all(is.na(distance)))
    return(0)
  if(all(riskcat==1))
    return(sum(majority(distance)))
  else if(all(riskcat==2))
    return(sum(high_risk(distance)))
  else
    return(sum(low_risk(distance)))
}

FFR_distances <- purrr::map_dfr(1:length(ldists),function(x) {dplyr::tibble(ix=x,Distance=ldists[[x]])}) %>% 
  dplyr::right_join(visit_num %>% dplyr::mutate(ix=1:dplyr::n())) %>%
  dplyr::filter(!is.na(Distance)) %>%
  dplyr::select(id,measurement,Distance)

FFR_distances %>% dplyr::right_join(visit_num) %>% dplyr::left_join(sjdf) %>% 
  dplyr::group_by(id,measurement) %>% 
  dplyr::summarise(exposure = FFR_effect(riskcat,Distance)) -> edf

sjdf <- dplyr::left_join(sjdf,visit_num) %>% 
  dplyr::left_join(edf) %>% 
  dplyr::mutate(year = measurement -1)


sjdf$BMI <- 33 +  sjdf$sex* -2.2 + .1*sjdf$year + sjdf$exposure + 
  sjdf$subj_effect + sjdf$subj_slope*sjdf$year + rnorm(num_obs)



longitudinal_clusters <- benvo(subject_data = sjdf,
                              bef_data = list(FFR_distances),
                              bef_names = "FFR",
                              distance_col = "Distance",
                              joining_id = c("id","measurement"))

usethis::use_data(longitudinal_clusters, overwrite = TRUE)
