test_that("detectPeaksOKseq returns correct output", {

  rfd = readRDS(system.file("testdata", "rfd_profile.RDS", package = "replicationOrigins"))
  res = readRDS(system.file("testdata", "rfd_ori_pos.RDS", package = "replicationOrigins"))

  rfd_ori_pos = detectPeaksOKseq(rfd, Verbose = FALSE)

  expect_identical(rfd_ori_pos, res)
})
