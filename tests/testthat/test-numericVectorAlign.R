test_that("numericVectorAlign returns correct output", {

  pos_sets = readRDS(system.file("testdata", "coordinate_sets.RDS", package = "replicationOrigins"))
  res = readRDS(system.file("testdata", "nva_res.RDS", package = "replicationOrigins"))

  pos_sets_num = lapply(pos_sets, function(x) x@ranges@start)
  nva_res = numericVectorAlign(pos_sets_num$set1,pos_sets_num$set2, 100000)

  expect_identical(nva_res, res)
})
