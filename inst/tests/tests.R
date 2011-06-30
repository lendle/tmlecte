context("Just some function checks")

test_that("qlogis and plogis work as I expect", {
  x <- -10:10
  y <- seq(.1, .9, .1)
  expect_that(logistic(x), equals(1/(1+exp(-x))))
  expect_that(logit(y), equals(log(y/(1-y))))
})

test_that("regress and predict.regress", {
  X <- rnorm(100)
  Y <- rnorm(100, 10+X, .1)
  Ybin <- rbinom(100, 1, plogis(X))
  expect_that(r1 <- regress(Y, X, family=gaussian), is_a("regress"))
  expect_that(r2 <- regress(Ybin, X), is_a("regress"))
  expect_that(predict(r1), is_a("numeric"))
  expect_that(pr2 <- predict(r2), is_a("numeric"))
  expect_that(all(pr2 >= 0 & pr2 <=1), is_true())
  if (!("SuperLearner" %in% installed.packages()[,1L])) {
    expect_that(rnsl <- regress(Y, X, family=gaussian, method="SL", SL.library="SL.glm"), gives_warning("SuperLearner is not installed"))
    expect_that(rnsl$method, equals("glm"))
  } else {
    SL.version <- packageVersion("SuperLearner")$major
    if (SL.version==1) {
      expect_that(rsl <- regress(Y, X, method="SL", family=gaussian(), SL.library="SL.glm"), gives_warning("SuperLearner is out of date"))
    }
    expect_that(predict(rsl), is_a("matrix"))
  }
})


dat <- gendata(100)

t <- tmle.cte(dat$A, dat[,c("W1","W2", "W3")], dat$Y, a=1, family=binomial, Q.method="SL")
