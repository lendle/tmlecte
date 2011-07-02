context("Just some function checks")

test_that("qlogis and plogis work as I expect", {
  x <- -10:10
  y <- seq(.1, .9, .1)
  expect_that(logistic(x), equals(1/(1+exp(-x))))
  expect_that(logit(y), equals(log(y/(1-y))))
})

test_that("regress and predict.regress", {
  X <- matrix(rnorm(200), 100, 2)
  Y <- rnorm(100, 10+X[,1], .1)
  Ybin <- rbinom(100, 1, plogis(X[,2]))
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
      expect_that(rsl <- regress(Ybin, X, method="SL", family=binomial(), SL.library=c("SL.glm", "SL.knn")), gives_warning("SuperLearner is out of date"))
    } else {
      expect_that(rsl <- regress(Ybin, X, method="SL", family=binomial(), SL.library=c("SL.glm", "SL.knn")), is_a("regress"))
    }
    expect_that(predict(rsl), is_a("matrix"))
    expect_that(predict(rsl, newdata=X, X=X, Y=Ybin), is_a("matrix"))    
  }
})

