capture_output(
  xs <-  fdp_staplm(BMI ~ sex + sap(FFR),
                          benvo = FFR_benvo,
                          K = 5, ## for speed
                          iter_max = 30,
                          burn_in = 20)
)

test_that("nobs", 
  expect_equal(nobs(xs),length(FFR_benvo@subject_data$BMI))
)
