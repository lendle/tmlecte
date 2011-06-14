context("Just some function checks")

test_that("qlogis and plogis work as I expect", {
  x <- -10:10
  y <- seq(.1, .9, .1)
  expect_that(logistic(x), equals(1/(1+exp(-x))))
  expect_that(logit(y), equals(log(y/(1-y))))
})

