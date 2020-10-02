# capture_output(
#   xs <-  fdp_staplm(BMI ~ sex + sap(FFR),
#                           benvo = FFR_benvo,
#                           K = 5, ## for speed
#                           iter_max = 30,
#                           burn_in = 20)
# )

# capture_output(
#   lmer1 <- fdp_staplmer(BMI ~ sex + year + sap(FFR) + (year|ID),
#                         benvo = longitudinal_clusters,
#                         iter_max = 30, burn_in = 20)
# )

# test_that("nobs",
# expect_equal(nobs(xs),nrow(FFR_benvo$subject_data)),
  # expect_equal(nobs(lmer1),length(longitudinal_clusters@subject_data$sex))
# )
