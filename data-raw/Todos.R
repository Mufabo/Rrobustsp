# debug ----
'
ladlasso(randn(5), matrix(randn(10), 5, 2), 0.5)

> ranklasso(c(-0.7411  , -0.5078  , -0.3206    ,0.0125  , -3.0292),
matrix(c(-0.4570  ,  1.2424  , -1.0667 ,   0.9337 ,   0.3503)), 0.5)
'

# mising examples ----
'

  prepare_Rd: m_param_est.Rd:32-34: Dropping empty section \examples
  prepare_Rd: rankflasso.Rd:33-35: Dropping empty section \examples
  prepare_Rd: rankflassopath.Rd:36-38: Dropping empty section \examples
  prepare_Rd: ranklasso.Rd:30-32: Dropping empty section \examples
'
# ----
'
> checking R code for possible problems ... NOTE
  ar_est_bip_s: no visible binding for global variable 'ind_max2'
  ar_est_bip_s: no visible binding for global variable 'temp2'
  ar_est_bip_tau: no visible binding for global variable 'ind_max2'
  ar_est_bip_tau: no visible binding for global variable 'temp2'
  arma_est_bip_s: no visible binding for global variable 'a_bip_sc'
  arma_est_bip_s: no visible binding for global variable 'x_filt'
  arma_est_bip_tau: no visible binding for global variable 'a_bip_sc'
  arma_est_bip_tau: no visible binding for global variable 'x_filt'
  bip_ar1_s: no visible global function definition for 'polyval'
  bip_ar1_s: no visible global function definition for 'polyfit'
  bip_ar1_tau: no visible global function definition for 'polyval'
  bip_ar1_tau: no visible global function definition for 'polyfit'
  bip_tau_arma_PSD_book: no visible global function definition for
    'linspace'
  bip_tau_arma_PSD_book: no visible global function definition for
    'polyval'
  biweight_filter: no visible binding for global variable 'xFRM'
  biweight_filter: no visible binding for global variable 'ARM'
  biweight_filter: no visible binding for global variable 'BRM'
  biweight_filter: no visible global function definition for 'nls'
  biweight_filter: no visible binding for global variable 'N_Prime'
  biweight_filter: no visible binding for global variable 'xFb2'
  create_environment_book: no visible global function definition for
    'rnorm'
  ekf_toa: no visible global function definition for 'to.tensor'
  ekf_toa: no visible global function definition for 'inv'
  ekf_toa_Masreliez: no visible global function definition for
    'to.tensor'
  ekf_toa_Masreliez: no visible global function definition for 'inv'
  ekf_toa_Masreliez: no visible binding for global variable 'v'
  ekf_toa_Masreliez: no visible binding for global variable 'vp'
  ekf_toa_robust: no visible global function definition for 'fzero'
  ekf_toa_robust: no visible global function definition for 'to.tensor'
  ekf_toa_robust: no visible global function definition for 'inv'
  eval_track: no visible global function definition for 'plot'
  eval_track: no visible global function definition for 'lines'
  eval_track: no visible global function definition for 'legend'
  eval_track: no visible global function definition for 'grid'
  eval_track: no visible global function definition for 'abline'
  eval_track: no visible global function definition for 'ecdf'
  hublasso: no visible global function definition for 'pchisq'
  hublasso: no visible global function definition for 'rhofun'
  hublasso: no visible global function definition for 'sprinf'
  hubreg: no visible global function definition for 'pchisq'
  loss_Huber: no visible global function definition for 'qchisq'
  loss_Huber: no visible global function definition for 'pchisq'
  markov_chain_book: no visible global function definition for 'runif'
  markov_chain_book: no visible global function definition for 'rnorm'
  mat_stem: no visible global function definition for 'plot'
  mat_stem: no visible global function definition for 'lines'
  order_wk: no visible global function definition for 'tail'
  order_wk: no visible binding for global variable 'PSD_sort'
  orth: no visible global function definition for 'rankMatrix'
  rankflasso: no visible global function definition for 'sparseMatrix'
  rankflasso: no visible global function definition for 'sparseVector'
  repeated_median_filter: no visible binding for global variable 'w'
  repeated_median_filter: no visible binding for global variable 'yt'
  repeated_median_filter: no visible binding for global variable
    'N_Prime'
  res_scale_approx: no visible global function definition for 'polyfit'
  res_scale_approx: no visible global function definition for 'polyval'
  robust_starting_point: no visible global function definition for
    'arima'
  save_mat_to_rdata: no visible global function definition for 'readMat'
  spec_arma_est_bip_mm: no visible global function definition for
    'linspace'
  spec_arma_est_bip_mm: no visible global function definition for
    'polyval'
  spec_arma_est_bip_s: no visible global function definition for
    'linspace'
  spec_arma_est_bip_s: no visible global function definition for
    'polyval'
  spec_arma_est_bip_tau: no visible global function definition for
    'linspace'
  spec_arma_est_bip_tau: no visible global function definition for
    'polyval'
  split_into_prime: no visible global function definition for 'data'
  split_into_prime: no visible binding for global variable
    'prime_numbers'
  split_into_prime: no visible global function definition for 'tail'
  tloss_consistency_factor: no visible global function definition for
    'integrate'
  Undefined global functions or variables:
    ARM BRM N_Prime PSD_sort a_bip_sc abline arima data ecdf fzero grid
    ind_max2 integrate inv legend lines linspace nls pchisq plot polyfit
    polyval prime_numbers qchisq rankMatrix readMat rhofun rnorm runif
    sparseMatrix sparseVector sprinf tail temp2 to.tensor v vp w xFRM
    xFb2 x_filt yt
  Consider adding
    importFrom("graphics", "abline", "grid", "legend", "lines", "plot")
    importFrom("stats", "arima", "ecdf", "integrate", "nls", "pchisq",
               "qchisq", "rnorm", "runif")
    importFrom("utils", "data", "tail")
  to your NAMESPACE file.

'
# ----

# ----


'

Undocumented arguments in documentation object arma_s_resid'
     'x' 'beta_hat' 'p' 'q'
   Documented arguments not in \usage in documentation object 'arma_s_resid':
     'x:' 'beta_hat:' 'p:' 'q:'

   Undocumented arguments in documentation object 'arma_tau_resid_sc'
     'x' 'beta_hat' 'p' 'q'
   Documented arguments not in \usage in documentation object 'arma_tau_resid_sc':
     'x:' 'beta_hat:' 'p:' 'q:'

   Undocumented arguments in documentation object 'bip_ar1_s'
     'x' 'N' 'phi_grid' 'fine_grid' 'kap2'

   Undocumented arguments in documentation object 'bip_ar1_tau'
     'x' 'N' 'phi_grid' 'fine_grid' 'kap2'

   Undocumented arguments in documentation object 'bip_resid'
     'x' 'beta_hat' 'p' 'q'

   Undocumented arguments in documentation object 'bip_s_resid'
     'x' 'beta_hat' 'p' 'q'

   Undocumented arguments in documentation object 'bip_s_resid_sc'
     'x' 'beta_hat' 'p' 'q'

   Undocumented arguments in documentation object 'bip_tau_arma_PSD_book'
     'x' 'p' 'q'

   Undocumented arguments in documentation object 'bip_tau_resid_sc'
     'x' 'beta_hat' 'p' 'q'

   Undocumented arguments in documentation object 'biweight_filter'
     'x'
   Documented arguments not in \usage in documentation object 'biweight_filter':
     'x:'

   Undocumented arguments in documentation object 'create_environment_book'
     'parameter' 'start' 'sigma_v'

   Undocumented arguments in documentation object 'ekf_toa_Masreliez'
     'r_ges' 'theta_init' 'BS' 'parameter'

   Undocumented arguments in documentation object 'enet'
     'printitn' 'itermax'
   Documented arguments not in \usage in documentation object 'enet':
     'printitn:'

   Undocumented arguments in documentation object 'enetpath'
     'y' 'X' 'alpha' 'L' 'eps' 'intcpt' 'printitn'

   Undocumented arguments in documentation object 'eta'
     'x' 'c'
   Documented arguments not in \usage in documentation object 'eta':
     'x:' 'c:'

   Undocumented arguments in documentation object 'hublasso'
     'y' 'X' 'c' 'lambda' 'b0' 'sig0' 'reltol' 'printitn' 'iter_max'

   Undocumented arguments in documentation object 'hublassopath'
     'y' 'X' 'c' 'eps' 'intcpt' 'printitn'
   Documented arguments not in \usage in documentation object 'hublassopath':
     'y:' 'X:' 'c:' 'intcpt:' 'eps:' 'printitn:'

   Undocumented arguments in documentation object 'hubreg'
     'y' 'X' 'c' 'sig0' 'b0' 'printitn' 'iter_max' 'errortol'
   Documented arguments not in \usage in documentation object 'hubreg':
     'y:' 'X:' 'c:' 'sig0:' 'b0:' 'printitn:' 'iter_max:' 'errortol:'

   Undocumented arguments in documentation object 'inf_norm'
     'mat'

   Undocumented arguments in documentation object 'ladlassopath'
     'eps' 'intcpt' 'printitn'
   Documented arguments not in \usage in documentation object 'ladlassopath':
     'intcpt:' 'eps:' 'printitn:'

   Undocumented arguments in documentation object 'ladreg'
     'y' 'X' 'intcpt' 'b0' 'printitn'
   Documented arguments not in \usage in documentation object 'ladreg':
     'y:' 'X:' 'intcpt:' 'b0:' 'printitn:'

   Undocumented arguments in documentation object 'm_param_est'
     'receive' 'C' 'Theta' 'param'
   Documented arguments not in \usage in documentation object 'm_param_est':
     'recieve:' 'Theta:' 'C:' 'param.maxiters:' 'param.break:'

   Undocumented arguments in documentation object 'madn'
     'y'
   Documented arguments not in \usage in documentation object 'madn':
     'y:'

   Undocumented arguments in documentation object 'mat_sign'
     'x'

   Undocumented arguments in documentation object 'mat_stem'
     'x' 'y' 'pch' 'linecol' 'clinecol' '...'

   Undocumented arguments in documentation object 'muler_rho2'
     'x'

   Undocumented arguments in documentation object 'order_wk'
     'y'

   Undocumented arguments in documentation object 'psihub'
     'c'
   Documented arguments not in \usage in documentation object 'psihub':
     'threshold'

   Undocumented arguments in documentation object 'ranklassopath'
     'y' 'X' 'L' 'eps' 'reltol' 'printitn'
   Documented arguments not in \usage in documentation object 'ranklassopath':
     'y:' 'X:' 'L:' 'eps:' 'reltol:' 'printitn:'

   Undocumented arguments in documentation object 'repeated_median_filter'
     'x'
   Documented arguments not in \usage in documentation object 'repeated_median_filter':
     'x:'

   Undocumented arguments in documentation object 'repmat'
     'v' 's'

   Undocumented arguments in documentation object 'rladreg'
     'y' 'X' 'b0' 'printitn'
   Documented arguments not in \usage in documentation object 'rladreg':
     'y:' 'X:' 'b0:' 'printitn:'

   Undocumented arguments in documentation object 'robust_starting_point'
     'x' 'p' 'q' 'recursion_num'
   Documented arguments not in \usage in documentation object 'robust_starting_point':
     'x:' 'p:' 'q:'

   Undocumented arguments in documentation object 'signcm'
     'x'
   Documented arguments not in \usage in documentation object 'signcm':
     'X'

   Undocumented arguments in documentation object 'spec_arma_est_bip_mm'
     'x' 'p' 'q' 'tolx'

   Undocumented arguments in documentation object 'spec_arma_est_bip_s'
     'x' 'p' 'q' 'tolx'

   Undocumented arguments in documentation object 'spec_arma_est_bip_tau'
     'x' 'p' 'q' 'tolx'

   Undocumented arguments in documentation object 'sqrtm'
     'mat'

     Mloc: no visible global function definition for 'median'
     MlocHUB: no visible global function definition for 'median'
     MlocTUK: no visible global function definition for 'median'
     Mlocscale: no visible global function definition for 'median'
     Mreg: no visible global function definition for 'median'
     MscaleHUB: no visible global function definition for 'median'
     MscaleTUK: no visible global function definition for 'median'

     ar_est_bip_s: no visible binding for global variable 'ind_max2'
     ar_est_bip_s: no visible binding for global variable 'temp2'

     ar_est_bip_tau: no visible binding for global variable 'ind_max2'
     ar_est_bip_tau: no visible binding for global variable 'temp2'
     arma_est_bip_m: no visible global function definition for 'lsqnonlin'

     arma_est_bip_s: no visible binding for global variable 'a_bip_sc'
     arma_est_bip_s: no visible binding for global variable 'x_filt'

     arma_est_bip_tau: no visible binding for global variable 'a_bip_sc'
     arma_est_bip_tau: no visible binding for global variable 'x_filt'
     arma_s_resid: no visible global function definition for 'or'
     arma_s_resid: no visible global function definition for 'roots'

     bip_ar1_s: no visible global function definition for 'polyval'
     bip_ar1_s: no visible global function definition for 'polyfit'
     bip_ar1_tau: no visible global function definition for 'polyval'
     bip_ar1_tau: no visible global function definition for 'polyfit'
     bip_resid: no visible global function definition for 'or'
     bip_resid: no visible global function definition for 'roots'
     bip_s_resid: no visible global function definition for 'roots'
     bip_s_resid_sc: no visible global function definition for 'roots'
     bip_tau_arma_PSD_book: no visible global function definition for
     'median'
     bip_tau_arma_PSD_book: no visible global function definition for
     'linspace'
     bip_tau_arma_PSD_book: no visible global function definition for
     'polyval'
     bip_tau_resid_sc: no visible global function definition for 'roots'
     biweight_filter: no visible global function definition for '%<-%'
     biweight_filter: no visible binding for global variable 'xFRM'
     biweight_filter: no visible binding for global variable 'ARM'
     biweight_filter: no visible binding for global variable 'BRM'
     biweight_filter: no visible global function definition for 'nls'
     biweight_filter: no visible binding for global variable 'N_Prime'
     biweight_filter: no visible binding for global variable 'xFb2'
     create_environment_book: no visible global function definition for
     'rnorm'
     ekf_toa: no visible global function definition for 'to.tensor'
     ekf_toa: no visible global function definition for 'inv'
     ekf_toa_Masreliez: no visible global function definition for
     'to.tensor'
     ekf_toa_Masreliez: no visible global function definition for 'inv'
     ekf_toa_Masreliez: no visible global function definition for '%<-%'
     ekf_toa_Masreliez: no visible binding for global variable 'v'
     ekf_toa_Masreliez: no visible binding for global variable 'vp'
     ekf_toa_robust: no visible global function definition for 'fzero'
     ekf_toa_robust: no visible global function definition for 'to.tensor'
     ekf_toa_robust: no visible global function definition for 'inv'
     ekf_toa_robust: no visible global function definition for 'ginv'
     eval_track: no visible global function definition for 'plot'
     eval_track: no visible global function definition for 'lines'
     eval_track: no visible global function definition for 'legend'
     eval_track: no visible global function definition for 'grid'
     eval_track: no visible global function definition for 'abline'
     eval_track: no visible global function definition for 'ecdf'
     hublasso: no visible global function definition for 'pchisq'
     hublasso: no visible global function definition for 'rhofun'
     hublasso: no visible global function definition for 'sprinf'
     hubreg: no visible global function definition for 'pchisq'
     hubreg: no visible global function definition for 'ginv'
     ladlassopath: no visible global function definition for 'median'
     loss_Huber: no visible global function definition for 'qchisq'
     loss_Huber: no visible global function definition for 'pchisq'

     markov_chain_book: no visible global function definition for 'runif'
     markov_chain_book: no visible global function definition for 'rnorm'
     mat_stem: no visible global function definition for 'plot'
     mat_stem: no visible global function definition for 'lines'
     order_wk: no visible global function definition for 'median'
     order_wk: no visible binding for global variable 'median'
     order_wk: no visible global function definition for 'head'
     order_wk: no visible global function definition for 'tail'
     order_wk: no visible global function definition for '%<-%'
     order_wk: no visible binding for global variable 'PSD_sort'
     orth: no visible global function definition for 'rankMatrix'
     rankflasso: no visible global function definition for 'sparseMatrix'
     rankflasso: no visible global function definition for 'sparseVector'
     ranklassopath: no visible global function definition for 'median'
     repeated_median_filter: no visible global function definition for
     'median'
     repeated_median_filter: no visible binding for global variable 'w'
     repeated_median_filter: no visible binding for global variable 'yt'
     repeated_median_filter: no visible binding for global variable 'median'
     repeated_median_filter: no visible binding for global variable
     'N_Prime'
     res_scale_approx: no visible global function definition for 'polyfit'
     res_scale_approx: no visible global function definition for 'polyval'
     rladreg: no visible global function definition for 'ginv'
     robust_starting_point: no visible global function definition for
     'arima'
     robust_starting_point: no visible global function definition for
     'roots'
     robust_starting_point: no visible global function definition for '%<-%'
     save_mat_to_rdata: no visible global function definition for 'readMat'
     spatmed: no visible global function definition for 'median'
     spec_arma_est_bip_mm: no visible global function definition for
     'median'
     spec_arma_est_bip_mm: no visible global function definition for
     'linspace'
     spec_arma_est_bip_mm: no visible global function definition for
     'polyval'
     spec_arma_est_bip_s: no visible global function definition for 'median'
     spec_arma_est_bip_s: no visible global function definition for
     'linspace'
     spec_arma_est_bip_s: no visible global function definition for
     'polyval'
     spec_arma_est_bip_tau: no visible global function definition for
     'median'
     spec_arma_est_bip_tau: no visible global function definition for
     'linspace'
     spec_arma_est_bip_tau: no visible global function definition for
     'polyval'
     split_into_prime: no visible global function definition for 'data'
     split_into_prime: no visible binding for global variable
     'prime_numbers'
     split_into_prime: no visible global function definition for 'tail'
     tloss_consistency_factor: no visible global function definition for
     'integrate'
     wmed: no visible global function definition for 'median'
     Undefined global functions or variables:
        %<-% ARM BRM N_Prime PSD_sort a_bip_sc abline arima data ecdf fzero
     ginv grid head ind_max2 integrate inv legend lines linspace lsqnonlin
     mad median nls or pchisq plot polyfit polyval prime_numbers qchisq
     rankMatrix readMat rhofun rnorm roots runif sparseMatrix sparseVector
     sprinf tail temp2 to.tensor v vp w xFRM xFb2 x_filt yt


     > checking Rd files ... NOTE
     prepare_Rd: Mscat.Rd:42-44: Dropping empty section \examples
     prepare_Rd: ar_est_bip_s.Rd:27-29: Dropping empty section \examples
     prepare_Rd: ar_est_bip_tau.Rd:27-29: Dropping empty section \examples
     prepare_Rd: biweight_filter.Rd:35-37: Dropping empty section \examples
     prepare_Rd: ekf_toa.Rd:28-30: Dropping empty section \examples
     prepare_Rd: ekf_toa_robust.Rd:29-31: Dropping empty section \examples
     prepare_Rd: ladlasso.Rd:45-47: Dropping empty section \examples
     prepare_Rd: m_param_est.Rd:32-34: Dropping empty section \examples
     prepare_Rd: rankflasso.Rd:33-35: Dropping empty section \examples
     prepare_Rd: rankflassopath.Rd:36-38: Dropping empty section \examples
     prepare_Rd: ranklasso.Rd:30-32: Dropping empty section \examples
     prepare_Rd: robust_starting_point.Rd:28-30: Dropping empty section \examples
'
