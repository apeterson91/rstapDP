test_that("parameter checks", {
  expect_error(fdp_staplm(BMI ~ sex + sap(FFR),benvo = FFR_benvo,iter_max = 0,burn_in = 5))
  expect_error(fdp_staplm(BMI ~ sex + sap(FFR),benvo = FFR_benvo,iter_max = 3,burn_in = 5))
  expect_error(fdp_staplm(BMI ~ sex + sap(FFR),benvo = FFR_benvo,alpha_a = -3))
  expect_error(fdp_staplm(BMI ~ sex + sap(FFR),benvo = FFR_benvo,alpha_b = -3))
  expect_error(fdp_staplm(BMI ~ sex + sap(FFR),benvo = FFR_benvo,tau_a = -3))
  expect_error(fdp_staplm(BMI ~ sex + sap(FFR),benvo = FFR_benvo,tau_b = -3))
  expect_error(fdp_staplm(BMI ~ sex + sap(FFR),benvo = FFR_benvo,sigma_a = -3))
  expect_error(fdp_staplm(BMI ~ sex + sap(FFR),benvo = FFR_benvo,sigma_b = -3))
})

test_that("formula errors correctly",{
  expect_error(fdp_staplm(BMI ~ sex,benvo=FFR_benvo,iter_max = 5,burn_in = 1))
  expect_error(fdp_staplm(BMI ~ sex + sap(HFS),benvo=FFR_benvo,iter_max = 5, burn_in = 1))
  expect_error(fdp_staplm(BMI ~ sex + tap(FFR),benvo=FFR_benvo,iter_max = 5 ,burn_in = 1))
})
