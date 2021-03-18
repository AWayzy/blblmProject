test_that("blblm() returns a list of estimates of length B for each subsample 1 to m", {
  b <- 100
  m <- 20
  formula <- y ~ x
  data <- data.frame(x = rnorm(1000), y = rnorm(1000), stringsAsFactors = FALSE)
  # testing output lengths without a seed
  a <- blblm(y ~ x, data = data, m = m, B = b)
  c <- blblm(y ~ x, data = data, m = m, B = b, parallel = TRUE)

  # testing that function returns same values for set seed
  d1 <- blblm(y ~ x, data = data, m = m, B = b, seed = 100)
  d2 <- blblm(y ~ x, data = data, m = m, B = b, seed = 100)


  expect_equal(d1, d2)
  expect_equal(length(a$estimates), m)
  expect_equal(length(c$estimates), m)

  for (i in names(a$estimates)) {
    expect_equal(length(a$estimates[[i]]), b)
  }
  for (i in names(c$estimates)) {
    expect_equal(length(c$estimates[[i]]), b)
  }

  expect_equal(formula, a$formula)
  expect_equal(formula, c$formula)
})
