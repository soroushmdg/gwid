piggyback::pb_download(repo = "soroushmdg/gwid",
                       tag = "v0.0.1",
                       dest = tempdir())

case_control_data_file <- paste0(tempdir(),"//case-cont-RA.withmap.Rda")

case_control <- gwid::case_control(case_control_rda = case_control_data_file)


test_that("build section works", {
  expect_length(case_control,n = 6)
  expect_s3_class(case_control,"caco")
})
