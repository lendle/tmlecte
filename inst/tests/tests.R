context("Checking regress")

test_that("regress and predict.regress", {
  X <- matrix(rnorm(200), 100, 2)
  Y <- rnorm(100, 10+X[,1], .1)
  A <- Y
  Ybin <- rbinom(100, 1, plogis(X[,2]))
  expect_that(r1 <- regress(Y, X, formula= A ~ X, family=gaussian), is_a("regress"))
  expect_that(regress(A, X, formula= blah ~ X, family=gaussian), is_a("regress"))  
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

context("Checking tmle.cte")

test_that("tmle.cte handles specification of missing outcomes correctly", {
  set.seed(123)
  Y <- rbinom(100, 1, .5)
  A <- rbinom(100, 1, .5)
  B <- data.frame(W1=runif(100))
  Delta <- c(rep(1, 90), rep(0, 10))
  expect_that(tmle.cte(A, B, Y, Delta=Delta), is_a("cte"))
  Y[Delta==0] <- NA
  expect_that(t <- tmle.cte(A, B, Y, Delta=Delta), is_a("cte"))
  expect_that(t <- tmle.cte(A, B, Y), gives_warning("but Delta is NULL"))
  expect_that(t <- tmle.cte(A, B, Y, Delta=1), gives_warning("Delta is the wrong length"))
  expect_that(t <- tmle.cte(A, B, Y, Delta=rep(1,100)), gives_warning("Y are NA where Delta is 1"))  
})
