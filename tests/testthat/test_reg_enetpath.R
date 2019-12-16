library(Rrobustsp)

test_that('enetpath test 1', {
  # Arguments ----
  load('~/Rrobustsp/data/images.RData')
  load(path_test('enetpath_test_1_stats'))
  load(path_test('enetpath_test_1_B'))

  y <- unlist(unname(images['y20n']))
  alpha <- 1

  X <- diag(1, 400, 400)
  L <- 20
  alpha <- 1

  sol <- enetpath(y, X, alpha, L)
  B_sol <- sol[[1]]
  stats_sol <- sol[[2]]

  expect_equal(B_sol, matrix(enetpath_test_1_B, 401, 21))
  expect_equal(stats_sol[[1]], NULL)
  expect_equal(stats_sol[[2]], NULL)
  expect_equal(stats_sol[[3]], c(enetpath_test_1_stats[[3]]))
  })

test_that('enetpath test 2', {
  # Arguments ----
  load('~/Rrobustsp/data/images.RData')
  load(path_test('enetpath_test_2_stats'))
  load(path_test('enetpath_test_2_B'))

  y <- unlist(unname(images['y20n']))
  alpha <- 1

  X <- diag(1, 400, 400)
  L <- 20
  alpha <- 1

  sol <- enetpath(y, X, alpha, L, eps = 1e-3, intcpt = F)
  B_sol <- sol[[1]]
  stats_sol <- sol[[2]]

  load(path_test('enetpath_test_2_ero'))

  x <- matrix(c(-0.3034,0.8884,-0.8095,0.2939,-1.1471,-2.9443,-0.7873,-1.0689,1.4384), 3, 3, byrow = T)

  scaledata <- function(x){
    nom <- sweep(x, 2, apply(x, 2, min)) # subtract col_mins from each respective column
    denom <- apply(x, 2, max) - apply(x, 2, min)
    return(3 * sweep(nom, 2, denom, FUN = '/'))
  }

  B_sol2 <- B_sol[, 2:ncol(B_sol)]
  ero <- scaledata(B_sol2) - repmat(y, ncol(B_sol2))

  expect_equal(ero, enetpath_test_2_ero)
  expect_equal(B_sol, matrix(enetpath_test_2_B, 400, 21))
  expect_equal(stats_sol[[1]], NULL)
  expect_equal(stats_sol[[2]], NULL)
  expect_equal(stats_sol[[3]], enetpath_test_2_stats)
})
