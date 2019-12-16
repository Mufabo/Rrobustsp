library(R.matlab)

save_mat_to_rdata <- function(s){

  sol <- readMat(paste('C:/users/computer/desktop/testdaten_Rrobustsp/',s, '.mat', sep = ''))

  sol <- unlist(unname(sol))

  save(sol, file = paste('C:/users/computer/desktop/testdaten_Rrobustsp/', s, '.rdata', sep=""))
}
