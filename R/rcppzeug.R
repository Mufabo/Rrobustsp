cppFunction(depends='RcppArmadillo', code = '

  void convertSparse(S4 mat) {

    // obtain dim, i, p. x from S4 object
    IntegerVector dims = mat.slot("Dim");
    arma::urowvec i = Rcpp::as<arma::urowvec>(mat.slot("i"));
    arma::urowvec p = Rcpp::as<arma::urowvec>(mat.slot("p"));
    arma::vec x     = Rcpp::as<arma::vec>(mat.slot("x"));

    int nrow = dims[0], ncol = dims[1];

    // use Armadillo sparse matrix constructor
    arma::sp_mat res(i, p, x, nrow, ncol);
    Rcout << "SpMat res:\n" << res << std::endl;
  }')

X <- sparseMatrix(1:4000, 1:4000, x = rnorm(4000))

convertSparse(X)
