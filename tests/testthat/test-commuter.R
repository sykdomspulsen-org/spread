context("test-commuter")

test_that("seed works to keep results same", {
  set.seed(4)
  d1 <- commuter(r0 = 2)
  set.seed(4)
  d2 <- commuter(r0 = 2)

  testthat::expect_equal(d1, d2)
})

test_that("seed works to keep results different", {
  set.seed(4)
  d1 <- commuter(r0 = 2)
  set.seed(5)
  d2 <- commuter(r0 = 2)

  testthat::expect_false(isTRUE(all.equal(d1, d2)))
})
