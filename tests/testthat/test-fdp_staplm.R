test_that("formula checks", {
  expect_error(fdp_staplm(BMI ~ sex + sap(FFR)))
})
