eval_track <- function(ekf_th, ekf_Hc, param){

  # EKF ----

  bx <- ekf_th[1,] - param$thx
  by <- ekf_th[2,] - param$thy
  bvx <- ekf_th[3,] - param$thvx
  bvy <- ekf_th[4,] - param$thvy

  eex <- (ekf_th[1, ] - param$thx)^2
  eey <- (ekf_th[2, ] - param$thy)^2
  med <- sqrt(eex + eey)
  eevx <- (ekf_th[3, ] - param$thvx)^2
  eevy <- (ekf_th[4, ] - param$thvy)^2

  eval_filter <- list()

  eval_filter$rmsex <- sqrt(mean(eex[, param$discardN:param$N]))
  eval_filter$rmsey <- sqrt(mean(eey[, param$discardN:param$N]))

  eval_filter$rmsevx <- sqrt(mean(eevx[, param$discardN:param$N]))
  eval_filter$rmsevy <- sqrt(mean(eevy[, param$discardN:param$N]))

  eexy <- eex + eey

  eval_filter$med <- mean(med[, param$discardN:param$N])
  eval_filter$rmse <- sqrt(mean(eexy[, param$discardN:param$N]))
  eval_filter$biasx_v <- mean(bx[, param$discardN:param$N])
  eval_filter$biasy_v <- mean(by[, param$discardN:param$N])
  eval_filter$biasvx_v <- mean(bvx[, param$discardN:param$N])
  eval_filter$biasvy_v <- mean(bvy[, param$discardN:param$N])
  eval_filter$biasx <- mean(bx[, param$discardN:param$N])
  eval_filter$biasy <- mean(by[, param$discardN:param$N])
  eval_filter$biasvx <- mean(bvx[, param$discardN:param$N])
  eval_filter$biasvy <- mean(bvy[, param$discardN:param$N])

  MED <- med[, param$discardN:param$N]
  MED1 <- med

  thx <- ekf_th[1,]
  thy <- ekf_th[2,]

  # REKF ----

  bxR <- ekf_Hc[1,] - param$thx
  byR <- ekf_Hc[2,] - param$thy
  bvxR <- ekf_Hc[3,] - param$thvx
  bvyR <- ekf_Hc[4,] - param$thvy

  eexR <- (ekf_Hc[1, ] - param$thx)^2
  eeyR <- (ekf_Hc[2, ] - param$thy)^2
  medR <- sqrt(eexR + eeyR)
  eevxR <- (ekf_Hc[3, ] - param$thvx)^2
  eevyR <- (ekf_Hc[4, ] - param$thvy)^2

  eval_filterR <- list()

  eval_filterR$rmsex <- sqrt(mean(eexR[, param$discardN:param$N]))
  eval_filterR$rmsey <- sqrt(mean(eeyR[, param$discardN:param$N]))

  eval_filterR$rmsevx <- sqrt(mean(eevxR[, param$discardN:param$N]))
  eval_filterR$rmsevy <- sqrt(mean(eevyR[, param$discardN:param$N]))

  eexyR <- eexR + eeyR

  eval_filterR$med <- mean(medR[, param$discardN:param$N])
  eval_filterR$rmse <- sqrt(mean(eexyR[, param$discardN:param$N]))
  eval_filterR$biasx_v <- mean(bxR[, param$discardN:param$N])
  eval_filterR$biasy_v <- mean(byR[, param$discardN:param$N])
  eval_filterR$biasvx_v <- mean(bvxR[, param$discardN:param$N])
  eval_filterR$biasvy_v <- mean(bvyR[, param$discardN:param$N])
  eval_filterR$biasx <- mean(bxR[, param$discardN:param$N])
  eval_filterR$biasy <- mean(byR[, param$discardN:param$N])
  eval_filterR$biasvx <- mean(bvxR[, param$discardN:param$N])
  eval_filterR$biasvy <- mean(bvyR[, param$discardN:param$N])

  MEDR <- medR[, param$discardN:param$N]
  MED1R <- medR

  thxR <- ekf_Hc[1,]
  thyR <- ekf_Hc[2,]

  # fig 1 ----
  # plot top-down view of trajectory and base stations
  plot(param$BS[,1], param$BS[,2]
       , lwd = 2
       , xlab = 'x-position in m'
       , ylab = 'y-position in m'
       , col = 1
       , xlim = c(0, 12000)
       , ylim = c(2000, 14000))
  lines(param$thx, param$thy, lwd = 2, col = 1)
  lines(ekf_th[1,], ekf_th[2,], lwd = 2, col = 'blue')
  lines(ekf_Hc[1,], ekf_Hc[2,], lwd = 2, col = 'red')
  legend('topright',
         legend = c('BS', 'true', 'EKF', 'Robust EKF'),
         lty = 1, col = c(1,2, 'blue', 'red'))
  grid()

  # fig 2 ----

  indices <- colSums(param$noiseindices)
  plot(NULL, xlim = c(0,3000), ylim = c(0, 800), xlab = 'time', ylab = 'RMSE')
  abline(v = which(indices > 0), col = 'goldenrod1')


  lines(x = 1:3000, y = sqrt(eexy), col = 'blue')
  lines(x = 1:3000, y = sqrt(eexyR),col = 'red')
  legend('topright',
         legend = c('NLOS', 'EKF', 'REKF'),
         lty = 1, col = c('goldenrod1', 'blue', 'red'))

  grid()

  # fig 3 ----
  plot(ecdf(MED)
       , xlab = 'Error Distance in m'
       , ylab = 'Empirical CDF'
       , col = 'blue')
  lines(ecdf(MEDR), col = 'red')
  legend('topright',
         legend = c('EKF', 'REKF'),
         lty = 1, col = c('blue', 'red'))
  }
