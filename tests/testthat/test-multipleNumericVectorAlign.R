test_that("multipleNumericVectorAlign returns correct output", {

  pos_sets = readRDS(system.file("testdata", "coordinate_sets.RDS", package = "replicationOrigins"))
  res = readRDS(system.file("testdata", "mnva_res.RDS", package = "replicationOrigins"))

  pos_sets_num = lapply(pos_sets, function(x) x@ranges@start)
  mnva_res = multipleNumericVectorAlign(pos_sets_num, 100000, '-')

  expect_identical(mnva_res, res)
})
