library(Rrobustsp)
library(MASS)

test_that('ladlasso1 p>1, no intcpt', {

  data('images')
  load(path_test('Blas'))

  load(path_test('sol1'))
  y20n <- unlist(unname(images['y20n']))
  lambda1 <- 124
  intcpt <- F

  Rsol <- ladlasso(y20n, diag(1, 400, 400), lambda1, intcpt, Blas)

  R_sol <- Rsol[[1]]
  R_it  <- Rsol[[2]]


  expect_equal(R_it, 11)
  expect_equal(R_sol, sol)
})

test_that('ladlasso2 intcpt, p>1', {

  data('images')

  y20n <- unlist(unname(images['y20n']))
  lambda1 <- 124
  intcpt <- T

  load(path_test('sol2'))

  Rsol <- ladlasso(y20n, diag(1, 400, 400), lambda1, intcpt)

  R_it  <- Rsol[[2]]
  R_sol <- Rsol[[1]] # N x 1 matrix, whereas sol is N


  expect_equal(R_sol, sol)
  expect_equal(R_it, 20)
})

test_that('ladlasso 3, intcpt, p=1'{
# eye-check, is correct
})
