test_that("pmaORIdetection returns correct output", {

  mutGR = readRDS(system.file("testdata", "somatic_mutations.RDS", package = "replicationOrigins"))
  res = readRDS(system.file("testdata", "ori_pos.RDS", package = "replicationOrigins"))

  ori_pos = pmaORIdetection(mutGR,
                            useChrs='chr22')

  expect_identical(ori_pos, res)
})
