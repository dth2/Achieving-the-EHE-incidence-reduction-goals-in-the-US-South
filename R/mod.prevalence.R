
#' @title Prevalence Calculations within Time Steps
#'
#' @description This module calculates demographic, transmission, and clinical
#'              statistics at each time step within the simulation.
#'
#' @inheritParams aging_msm
#'
#' @details
#' Summary statistic calculations are of two broad forms: prevalence and
#' incidence. This function establishes the summary statistic vectors for both
#' prevalence and incidence at time 1, and then calculates the prevalence
#' statistics for times 2 onward. Incidence statistics (e.g., number of new
#' infections or deaths) are calculated within the modules as they depend on
#' vectors that are not stored external to the module.
#'
#' @return
#' This function returns the \code{dat} object with an updated summary of current
#' attributes stored in \code{dat$epi}.
#'
#' @keywords module msm
#'
#' @export
#'
prevalence_msm <- function(dat, at) {

  active <- dat$attr$active
  status <- dat$attr$status
  diag.status <- dat$attr$diag.status
  diag.stage <- dat$attr$diag.stage
  diag.time <- dat$attr$diag.time
  aids.time <- dat$attr$aids.time
  inf.time <- dat$attr$inf.time
  race <- dat$attr$race
  age <- dat$attr$age
  sex <- dat$attr$sex
  msm <- dat$attr$msm
  dem.cat <- dat$attr$dem.cat
  tx.init.time <- dat$attr$tx.init.time
  tx.status <- dat$attr$tx.status
  vl <- dat$attr$vl
  vl.last.usupp <- dat$attr$vl.last.usupp
  last.neg.test <- dat$attr$last.neg.test
  stage <- dat$attr$stage

  prepElig <- dat$attr$prepElig
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass

  rGC <- dat$attr$rGC
  uGC <- dat$attr$uGC
  rCT <- dat$attr$rCT
  uCT <- dat$attr$uCT

  # Pop Size / Demog
  dat$epi$num[at] <- sum(active == 1, na.rm = TRUE)
  dat$epi$num.B[at] <- sum(race == 1, na.rm = TRUE)
  dat$epi$num.H[at] <- sum(race == 2, na.rm = TRUE)
  dat$epi$num.W[at] <- sum(race == 3, na.rm = TRUE)
  dat$epi$age.mean[at] <- mean(age, na.rm = TRUE)
  dat$epi$s.num[at] <- sum(status == 0, na.rm = TRUE)
  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)
  dat$epi$i.num.B[at] <- sum(status == 1 & race == 1, na.rm = TRUE)
  dat$epi$i.num.H[at] <- sum(status == 1 & race == 2, na.rm = TRUE)
  dat$epi$i.num.W[at] <- sum(status == 1 & race == 3, na.rm = TRUE)

  dat$epi$i.num.B.msm[at] <- sum(status == 1 & dem.cat == 1, na.rm = TRUE)
  dat$epi$i.num.H.msm[at] <- sum(status == 1 & dem.cat == 2, na.rm = TRUE)
  dat$epi$i.num.W.msm[at] <- sum(status == 1 & dem.cat == 3, na.rm = TRUE)
  dat$epi$i.num.B.m.het[at] <- sum(status == 1 & dem.cat == 4, na.rm = TRUE)
  dat$epi$i.num.H.m.het[at] <- sum(status == 1 & dem.cat == 5, na.rm = TRUE)
  dat$epi$i.num.W.m.het[at] <- sum(status == 1 & dem.cat == 6, na.rm = TRUE)
  dat$epi$i.num.B.f.het[at] <- sum(status == 1 & dem.cat == 7, na.rm = TRUE)
  dat$epi$i.num.H.f.het[at] <- sum(status == 1 & dem.cat == 8, na.rm = TRUE)
  dat$epi$i.num.W.f.het[at] <- sum(status == 1 & dem.cat == 9, na.rm = TRUE)

  dat$epi$i.num.MSM[at] <- sum(status == 1 & msm == 1, na.rm = TRUE)
  dat$epi$i.num.HET[at] <- sum(status == 1 & msm == 0, na.rm = TRUE)

  dat$epi$i.num.dx[at] <- sum(diag.status == 1, na.rm = TRUE)

  # Prev / Incid
  dat$epi$i.prev[at] <- dat$epi$i.num[at] / dat$epi$num[at]
  dat$epi$i.prev.B[at] <- sum(race == 1 & status == 1, na.rm = TRUE) / sum(race == 1, na.rm = TRUE)
  dat$epi$i.prev.H[at] <- sum(race == 2 & status == 1, na.rm = TRUE) / sum(race == 2, na.rm = TRUE)
  dat$epi$i.prev.W[at] <- sum(race == 3 & status == 1, na.rm = TRUE) / sum(race == 3, na.rm = TRUE)

  dat$epi$i.prev.B.msm[at] <- sum(dem.cat == 1 & status == 1, na.rm = TRUE) / sum(dem.cat == 1, na.rm = TRUE)
  dat$epi$i.prev.H.msm[at] <- sum(dem.cat == 2 & status == 1, na.rm = TRUE) / sum(dem.cat == 2, na.rm = TRUE)
  dat$epi$i.prev.W.msm[at] <- sum(dem.cat == 3 & status == 1, na.rm = TRUE) / sum(dem.cat == 3, na.rm = TRUE)
  dat$epi$i.prev.B.m.het[at] <- sum(dem.cat == 4 & status == 1, na.rm = TRUE) / sum(dem.cat == 4, na.rm = TRUE)
  dat$epi$i.prev.H.m.het[at] <- sum(dem.cat == 5 & status == 1, na.rm = TRUE) / sum(dem.cat == 5, na.rm = TRUE)
  dat$epi$i.prev.W.m.het[at] <- sum(dem.cat == 6 & status == 1, na.rm = TRUE) / sum(dem.cat == 6, na.rm = TRUE)
  dat$epi$i.prev.B.f.het[at] <- sum(dem.cat == 7 & status == 1, na.rm = TRUE) / sum(dem.cat == 7, na.rm = TRUE)
  dat$epi$i.prev.H.f.het[at] <- sum(dem.cat == 8 & status == 1, na.rm = TRUE) / sum(dem.cat == 8, na.rm = TRUE)
  dat$epi$i.prev.W.f.het[at] <- sum(dem.cat == 9 & status == 1, na.rm = TRUE) / sum(dem.cat == 9, na.rm = TRUE)

  dat$epi$i.prev.adol[at] <- sum(age < 18 & status == 1, na.rm = TRUE) / sum(age < 18, na.rm = TRUE)

  dat$epi$i.prev.prop.B.msm[at] <- sum(dem.cat == 1 & status == 1, na.rm = TRUE) / sum(status == 1, na.rm = TRUE)
  dat$epi$i.prev.prop.H.msm[at] <- sum(dem.cat == 2 & status == 1, na.rm = TRUE) / sum(status == 1, na.rm = TRUE)
  dat$epi$i.prev.prop.W.msm[at] <- sum(dem.cat == 3 & status == 1, na.rm = TRUE) / sum(status == 1, na.rm = TRUE)
  dat$epi$i.prev.prop.B.m.het[at] <- sum(dem.cat == 4 & status == 1, na.rm = TRUE) / sum(status == 1, na.rm = TRUE)
  dat$epi$i.prev.prop.H.m.het[at] <- sum(dem.cat == 5 & status == 1, na.rm = TRUE) / sum(status == 1, na.rm = TRUE)
  dat$epi$i.prev.prop.W.m.het[at] <- sum(dem.cat == 6 & status == 1, na.rm = TRUE) / sum(status == 1, na.rm = TRUE)
  dat$epi$i.prev.prop.B.f.het[at] <- sum(dem.cat == 7 & status == 1, na.rm = TRUE) / sum(status == 1, na.rm = TRUE)
  dat$epi$i.prev.prop.H.f.het[at] <- sum(dem.cat == 8 & status == 1, na.rm = TRUE) / sum(status == 1, na.rm = TRUE)
  dat$epi$i.prev.prop.W.f.het[at] <- sum(dem.cat == 9 & status == 1, na.rm = TRUE) / sum(status == 1, na.rm = TRUE)



  dat$epi$i.prev.dx[at] <- sum(diag.status == 1, na.rm = TRUE) / dat$epi$num[at]
  dat$epi$i.prev.dx.B[at] <- sum(race == 1 & diag.status == 1, na.rm = TRUE) / sum(race == 1, na.rm = TRUE)
  dat$epi$i.prev.dx.H[at] <- sum(race == 2 & diag.status == 1, na.rm = TRUE) / sum(race == 2, na.rm = TRUE)
  dat$epi$i.prev.dx.W[at] <- sum(race == 3 & diag.status == 1, na.rm = TRUE) / sum(race == 3, na.rm = TRUE)

  dat$epi$i.prev.dx.B.msm[at] <- sum(dem.cat == 1 & diag.status == 1, na.rm = TRUE) / sum(dem.cat == 1, na.rm = TRUE)
  dat$epi$i.prev.dx.H.msm[at] <- sum(dem.cat == 2 & diag.status == 1, na.rm = TRUE) / sum(dem.cat == 2, na.rm = TRUE)
  dat$epi$i.prev.dx.W.msm[at] <- sum(dem.cat == 3 & diag.status == 1, na.rm = TRUE) / sum(dem.cat == 3, na.rm = TRUE)
  dat$epi$i.prev.dx.B.m.het[at] <- sum(dem.cat == 4 & diag.status == 1, na.rm = TRUE) / sum(dem.cat == 4, na.rm = TRUE)
  dat$epi$i.prev.dx.H.m.het[at] <- sum(dem.cat == 5 & diag.status == 1, na.rm = TRUE) / sum(dem.cat == 5, na.rm = TRUE)
  dat$epi$i.prev.dx.W.m.het[at] <- sum(dem.cat == 6 & diag.status == 1, na.rm = TRUE) / sum(dem.cat == 6, na.rm = TRUE)
  dat$epi$i.prev.dx.B.f.het[at] <- sum(dem.cat == 7 & diag.status == 1, na.rm = TRUE) / sum(dem.cat == 7, na.rm = TRUE)
  dat$epi$i.prev.dx.H.het[at] <- sum(dem.cat == 8 & diag.status == 1, na.rm = TRUE) / sum(dem.cat == 8, na.rm = TRUE)
  dat$epi$i.prev.dx.W.het[at] <- sum(dem.cat == 9 & diag.status == 1, na.rm = TRUE) / sum(dem.cat == 9, na.rm = TRUE)

  dat$epi$ir100[at] <- (dat$epi$incid[at] / sum(status == 0, dat$epi$incid[at], na.rm = TRUE)) * 5200

  dat$epi$ir100.M[at] <- (dat$epi$incid.M[at] / sum(status == 0 & sex == 1, dat$epi$incid.M[at], na.rm = TRUE)) * 5200
  dat$epi$ir100.F[at] <- (dat$epi$incid.F[at] / sum(status == 0 & sex == 2, dat$epi$incid.F[at], na.rm = TRUE)) * 5200

  dat$epi$ir100.MSM[at] <- (dat$epi$incid.MSM[at] / sum(status == 0 & msm == 1, dat$epi$incid.MSM[at], na.rm = TRUE)) * 5200
  dat$epi$ir100.HET[at] <- (dat$epi$incid.HET[at] / sum(status == 0 & msm == 0, dat$epi$incid.HET[at], na.rm = TRUE)) * 5200

  dat$epi$ir100.B[at] <- (dat$epi$incid.B[at] / sum(status == 0 & race == 1, dat$epi$incid.B[at], na.rm = TRUE)) * 5200
  dat$epi$ir100.H[at] <- (dat$epi$incid.H[at] / sum(status == 0 & race == 2, dat$epi$incid.H[at], na.rm = TRUE)) * 5200
  dat$epi$ir100.W[at] <- (dat$epi$incid.W[at] / sum(status == 0 & race == 3, dat$epi$incid.W[at], na.rm = TRUE)) * 5200

  dat$epi$ir100.B.msm[at] <- (dat$epi$incid.B.msm[at] / sum(status == 0 & dem.cat == 1, dat$epi$incid.B.msm[at], na.rm = TRUE)) * 5200
  dat$epi$ir100.H.msm[at] <- (dat$epi$incid.H.msm[at] / sum(status == 0 & dem.cat == 2, dat$epi$incid.H.msm[at], na.rm = TRUE)) * 5200
  dat$epi$ir100.W.msm[at] <- (dat$epi$incid.W.msm[at] / sum(status == 0 & dem.cat == 3, dat$epi$incid.W.msm[at], na.rm = TRUE)) * 5200
  dat$epi$ir100.B.m.het[at] <- (dat$epi$incid.B.m.het[at] / sum(status == 0 & dem.cat == 4, dat$epi$incid.B.m.het[at], na.rm = TRUE)) * 5200
  dat$epi$ir100.H.m.het[at] <- (dat$epi$incid.H.m.het[at] / sum(status == 0 & dem.cat == 5, dat$epi$incid.H.m.het[at], na.rm = TRUE)) * 5200
  dat$epi$ir100.W.m.het[at] <- (dat$epi$incid.W.m.het[at] / sum(status == 0 & dem.cat == 6, dat$epi$incid.W.m.het[at], na.rm = TRUE)) * 5200
  dat$epi$ir100.B.f.het[at] <- (dat$epi$incid.B.f.het[at] / sum(status == 0 & dem.cat == 7, dat$epi$incid.B.f.het[at], na.rm = TRUE)) * 5200
  dat$epi$ir100.H.f.het[at] <- (dat$epi$incid.H.f.het[at] / sum(status == 0 & dem.cat == 8, dat$epi$incid.H.f.het[at], na.rm = TRUE)) * 5200
  dat$epi$ir100.W.f.het[at] <- (dat$epi$incid.W.f.het[at] / sum(status == 0 & dem.cat == 9, dat$epi$incid.W.f.het[at], na.rm = TRUE)) * 5200

  dat$epi$ir100.adol[at] <- (dat$epi$incid.adol[at] / sum(status == 0 & age < 18, dat$epi$incid.adol[at], na.rm = TRUE)) * 5200


  # Care continuum stats (primary)
  dat$epi$cc.dx[at] <- sum(diag.status == 1 & inf.time >= 2, na.rm = TRUE) /
                       sum(status == 1 & inf.time >= 2, na.rm = TRUE)
  dat$epi$cc.dx.B[at] <- sum(diag.status == 1 & inf.time >= 2 & race == 1, na.rm = TRUE) /
                         sum(status == 1 & inf.time >= 2 & race == 1, na.rm = TRUE)
  dat$epi$cc.dx.H[at] <- sum(diag.status == 1 & inf.time >= 2 & race == 2, na.rm = TRUE) /
                         sum(status == 1 & inf.time >= 2 & race == 2, na.rm = TRUE)
  dat$epi$cc.dx.W[at] <- sum(diag.status == 1 & inf.time >= 2 & race == 3, na.rm = TRUE) /
                         sum(status == 1 & inf.time >= 2 & race == 3, na.rm = TRUE)

  dat$epi$cc.dx.B.msm[at] <- sum(diag.status == 1 & inf.time >= 2 & dem.cat == 1, na.rm = TRUE) /
                              sum(status == 1 & inf.time >= 2 & dem.cat == 1, na.rm = TRUE)
  dat$epi$cc.dx.H.msm[at] <- sum(diag.status == 1 & inf.time >= 2 & dem.cat == 2, na.rm = TRUE) /
                              sum(status == 1 & inf.time >= 2 & dem.cat == 2, na.rm = TRUE)
  dat$epi$cc.dx.W.msm[at] <- sum(diag.status == 1 & inf.time >= 2 & dem.cat == 3, na.rm = TRUE) /
                              sum(status == 1 & inf.time >= 2 & dem.cat == 3, na.rm = TRUE)
  dat$epi$cc.dx.B.m.het[at] <- sum(diag.status == 1 & inf.time >= 2 & dem.cat == 4, na.rm = TRUE) /
                              sum(status == 1 & inf.time >= 2 & dem.cat == 4, na.rm = TRUE)
  dat$epi$cc.dx.H.m.het[at] <- sum(diag.status == 1 & inf.time >= 2 & dem.cat == 5, na.rm = TRUE) /
                              sum(status == 1 & inf.time >= 2 & dem.cat == 5, na.rm = TRUE)
  dat$epi$cc.dx.W.m.het[at] <- sum(diag.status == 1 & inf.time >= 2 & dem.cat == 6, na.rm = TRUE) /
                              sum(status == 1 & inf.time >= 2 & dem.cat == 6, na.rm = TRUE)
  dat$epi$cc.dx.B.f.het[at] <- sum(diag.status == 1 & inf.time >= 2 & dem.cat == 7, na.rm = TRUE) /
                              sum(status == 1 & inf.time >= 2 & dem.cat == 7, na.rm = TRUE)
  dat$epi$cc.dx.H.f.het[at] <- sum(diag.status == 1 & inf.time >= 2 & dem.cat == 8, na.rm = TRUE) /
                              sum(status == 1 & inf.time >= 2 & dem.cat == 8, na.rm = TRUE)
  dat$epi$cc.dx.W.f.het[at] <- sum(diag.status == 1 & inf.time >= 2 & dem.cat == 9, na.rm = TRUE) /
                              sum(status == 1 & inf.time >= 2 & dem.cat == 9, na.rm = TRUE)

  ##DX.rate
  dat$epi$cc.dx.r.B.msm[at] <- (sum(diag.time <= at & diag.time >= at-52 & dem.cat == 1, na.rm = TRUE) /
                                  sum(dem.cat == 1, na.rm = TRUE)) * 100000
  dat$epi$cc.dx.r.H.msm[at] <- (sum(diag.time <= at & diag.time >= at-52 & dem.cat == 2, na.rm = TRUE) /
                                  sum(dem.cat == 2, na.rm = TRUE)) * 100000
  dat$epi$cc.dx.r.W.msm[at] <- (sum(diag.time <= at & diag.time >= at-52 & dem.cat == 3, na.rm = TRUE) /
                                  sum(dem.cat == 3, na.rm = TRUE)) * 100000
  dat$epi$cc.dx.r.B.m.het[at] <- (sum(diag.time <= at & diag.time >= at-52 & dem.cat == 4, na.rm = TRUE) /
                                    sum(dem.cat == 4, na.rm = TRUE)) * 100000
  dat$epi$cc.dx.r.H.m.het[at] <- (sum(diag.time <= at & diag.time >= at-52 & dem.cat == 5, na.rm = TRUE) /
                                    sum(dem.cat == 5, na.rm = TRUE)) * 100000
  dat$epi$cc.dx.r.W.m.het[at] <- (sum(diag.time <= at & diag.time >= at-52 & dem.cat == 6, na.rm = TRUE) /
                                    sum(dem.cat == 6, na.rm = TRUE)) * 100000
  dat$epi$cc.dx.r.B.f.het[at] <- (sum(diag.time <= at & diag.time >= at-52 & dem.cat == 7, na.rm = TRUE) /
                                    sum(dem.cat == 7, na.rm = TRUE)) * 100000
  dat$epi$cc.dx.r.H.f.het[at] <- (sum(diag.time <= at & diag.time >= at-52 & dem.cat == 8, na.rm = TRUE) /
                                    sum(dem.cat == 8, na.rm = TRUE)) * 100000
  dat$epi$cc.dx.r.W.f.het[at] <- (sum(diag.time <= at & diag.time >= at-52 & dem.cat == 9, na.rm = TRUE) /
                                    sum(dem.cat == 9, na.rm = TRUE)) * 100000


  dat$epi$cc.dx.aids[at] <- sum(diag.status == 1 & stage == 4 & inf.time >= 2 &
                                  aids.time - diag.time <= 52, na.rm = TRUE) /
                                  sum(diag.status == 1 & inf.time >= 2, na.rm = TRUE)
  dat$epi$cc.dx.aids.B[at] <- sum(diag.status == 1 & stage == 4 & inf.time >= 2 &
                                    aids.time - diag.time <= 52 & race == 1, na.rm = TRUE) /
                                    sum(diag.status == 1 & inf.time >= 2 & race == 1, na.rm = TRUE)
  dat$epi$cc.dx.aids.H[at] <- sum(diag.status == 1 & stage == 4 &  inf.time >= 2 &
                                    aids.time - diag.time <= 52 & race == 2, na.rm = TRUE) /
                                    sum(diag.status == 1 & inf.time >= 2 & race == 2, na.rm = TRUE)
  dat$epi$cc.dx.aids.W[at] <- sum(diag.status == 1 & stage == 4 & inf.time >= 2 &
                                    aids.time - diag.time <= 52 & race == 3, na.rm = TRUE) /
                                    sum(diag.status == 1 & inf.time >= 2 & race == 3, na.rm = TRUE)

  dat$epi$cc.dx.aids.B.msm[at] <- sum(diag.status == 1 & stage == 4 & inf.time >= 2 &
                                    aids.time - diag.time <= 52 & dem.cat == 1, na.rm = TRUE) /
                                    sum(diag.status == 1 & inf.time >= 2 & dem.cat == 1, na.rm = TRUE)
  dat$epi$cc.dx.aids.H.msm[at] <- sum(diag.status == 1 & stage == 4 &  inf.time >= 2 &
                                    aids.time - diag.time <= 52 & dem.cat == 2, na.rm = TRUE) /
                                     sum(diag.status == 1 & inf.time >= 2 & dem.cat == 2, na.rm = TRUE)
  dat$epi$cc.dx.aids.W.msm[at] <- sum(diag.status == 1 & stage == 4 & inf.time >= 2 &
                                    aids.time - diag.time <= 52 & dem.cat == 3, na.rm = TRUE) /
                                    sum(diag.status == 1 & inf.time >= 2 & dem.cat == 3, na.rm = TRUE)
  dat$epi$cc.dx.aids.B.m.het[at] <- sum(diag.status == 1 & stage == 4 & inf.time >= 2 &
                                    aids.time - diag.time <= 52 & dem.cat == 4, na.rm = TRUE) /
                                    sum(diag.status == 1 & inf.time >= 2 & dem.cat == 4, na.rm = TRUE)
  dat$epi$cc.dx.aids.H.m.het[at] <- sum(diag.status == 1 & stage == 4 &  inf.time >= 2 &
                                    aids.time - diag.time <= 52 & dem.cat == 5, na.rm = TRUE) /
                                   sum(diag.status == 1 & inf.time >= 2 & dem.cat == 5, na.rm = TRUE)
  dat$epi$cc.dx.aids.W.m.het[at] <- sum(diag.status == 1 & stage == 4 & inf.time >= 2 &
                                    aids.time - diag.time <= 52 & dem.cat == 6, na.rm = TRUE) /
                                   sum(diag.status == 1 & inf.time >= 2 & dem.cat == 6, na.rm = TRUE)
  dat$epi$cc.dx.aids.B.f.het[at] <- sum(diag.status == 1 & stage == 4 & inf.time >= 2 &
                                    aids.time - diag.time <= 52 & dem.cat == 7, na.rm = TRUE) /
                                    sum(diag.status == 1 & inf.time >= 2 & dem.cat == 7, na.rm = TRUE)
  dat$epi$cc.dx.aids.H.f.het[at] <- sum(diag.status == 1 & stage == 4 &  inf.time >= 2 &
                                    aids.time - diag.time <= 52 & dem.cat == 8, na.rm = TRUE) /
                                    sum(diag.status == 1 & inf.time >= 2 & dem.cat == 8, na.rm = TRUE)
  dat$epi$cc.dx.aids.W.f.het[at] <- sum(diag.status == 1 & stage == 4 & inf.time >= 2 &
                                    aids.time - diag.time <= 52 & dem.cat == 9, na.rm = TRUE) /
                                    sum(diag.status == 1 & inf.time >= 2 & dem.cat == 9, na.rm = TRUE)



  dat$epi$cc.linked1m[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 2, na.rm = TRUE) /
                             sum(dat$attr$diag.status == 1 & diag.time >= 2, na.rm = TRUE)
  dat$epi$cc.linked1m.B[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 2 & race == 1, na.rm = TRUE) /
                               sum(diag.status == 1 & diag.time >= 2 & race == 1, na.rm = TRUE)
  dat$epi$cc.linked1m.H[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 2 & race == 2, na.rm = TRUE) /
                               sum(diag.status == 1 & diag.time >= 2 & race == 2, na.rm = TRUE)
  dat$epi$cc.linked1m.W[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 2 & race == 3, na.rm = TRUE) /
                               sum(diag.status == 1 & diag.time >= 2 & race == 3, na.rm = TRUE)

  dat$epi$cc.linked1m.B.msm[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 2 & dem.cat == 1, na.rm = TRUE) /
                                    sum(diag.status == 1 & diag.time >= 2 & dem.cat == 1, na.rm = TRUE)
  dat$epi$cc.linked1m.H.msm[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 2 & dem.cat == 2, na.rm = TRUE) /
                                    sum(diag.status == 1 & diag.time >= 2 & dem.cat == 2, na.rm = TRUE)
  dat$epi$cc.linked1m.W.msm[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 2 & dem.cat == 3, na.rm = TRUE) /
                                    sum(diag.status == 1 & diag.time >= 2 & dem.cat == 3, na.rm = TRUE)
  dat$epi$cc.linked1m.B.m.het[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 2 & dem.cat == 4, na.rm = TRUE) /
                                    sum(diag.status == 1 & diag.time >= 2 & dem.cat == 4, na.rm = TRUE)
  dat$epi$cc.linked1m.H.m.het[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 2 & dem.cat == 5, na.rm = TRUE) /
                                    sum(diag.status == 1 & diag.time >= 2 & dem.cat == 5, na.rm = TRUE)
  dat$epi$cc.linked1m.W.m.het[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 2 & dem.cat == 6, na.rm = TRUE) /
                                    sum(diag.status == 1 & diag.time >= 2 & dem.cat == 6, na.rm = TRUE)
  dat$epi$cc.linked1m.B.f.het[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 2 & dem.cat == 7, na.rm = TRUE) /
                                    sum(diag.status == 1 & diag.time >= 2 & dem.cat == 7, na.rm = TRUE)
  dat$epi$cc.linked1m.H.f.het[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 2 & dem.cat == 7, na.rm = TRUE) /
                                    sum(diag.status == 1 & diag.time >= 2 & dem.cat == 7, na.rm = TRUE)
  dat$epi$cc.linked1m.W.f.het[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 2 & dem.cat == 7, na.rm = TRUE) /
                                    sum(diag.status == 1 & diag.time >= 2 & dem.cat == 7, na.rm = TRUE)


#  dat$epi$cc.linked1m.int[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 3380, na.rm = TRUE) /
#                                 sum(dat$attr$diag.status == 1 & diag.time >= 3380, na.rm = TRUE)
#  dat$epi$cc.linked1m.int.B[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 3380 & race == 1, na.rm = TRUE) /
#                                   sum(diag.status == 1 & diag.time >= 3380 & race == 1, na.rm = TRUE)
#  dat$epi$cc.linked1m.int.H[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 3380 & race == 2, na.rm = TRUE) /
#                                   sum(diag.status == 1 & diag.time >= 3380 & race == 2, na.rm = TRUE)
#  dat$epi$cc.linked1m.int.W[at] <- sum(tx.init.time - diag.time <= 4 & diag.time >= 3380 & race == 3, na.rm = TRUE) /
#                                   sum(diag.status == 1 & diag.time >= 3380 & race == 3, na.rm = TRUE)

  dat$epi$cc.vsupp[at] <- sum(vl <= log10(200) & diag.status == 1 & diag.time >= 2, na.rm = TRUE) /
                          sum(diag.status == 1 & diag.time >= 2, na.rm = TRUE)
  dat$epi$cc.vsupp.B[at] <- sum(vl <= log10(200) & diag.status == 1 &
                                diag.time >= 2 & race == 1, na.rm = TRUE) /
                            sum(diag.status == 1 & diag.time >= 2 & race == 1, na.rm = TRUE)
  dat$epi$cc.vsupp.H[at] <- sum(vl <= log10(200) & diag.status == 1 &
                                diag.time >= 2 & race == 2, na.rm = TRUE) /
                            sum(diag.status == 1 & diag.time >= 2 & race == 2, na.rm = TRUE)
  dat$epi$cc.vsupp.W[at] <- sum(vl <= log10(200) & diag.status == 1 &
                                diag.time >= 2 & race == 3, na.rm = TRUE) /
                            sum(diag.status == 1 & diag.time >= 2 & race == 3, na.rm = TRUE)


  dat$epi$cc.vsupp.B.msm[at] <- sum(vl <= log10(200) & diag.status == 1 &
                                  diag.time >= 2 & dem.cat == 1, na.rm = TRUE) /
                                  sum(diag.status == 1 & diag.time >= 2 & dem.cat == 1, na.rm = TRUE)
  dat$epi$cc.vsupp.H.msm[at] <- sum(vl <= log10(200) & diag.status == 1 &
                                  diag.time >= 2 & dem.cat == 2, na.rm = TRUE) /
                                  sum(diag.status == 1 & diag.time >= 2 & dem.cat == 2, na.rm = TRUE)
  dat$epi$cc.vsupp.W.msm[at] <- sum(vl <= log10(200) & diag.status == 1 &
                                  diag.time >= 2 & dem.cat == 3, na.rm = TRUE) /
                                  sum(diag.status == 1 & diag.time >= 2 & dem.cat == 3, na.rm = TRUE)
  dat$epi$cc.vsupp.B.m.het[at] <- sum(vl <= log10(200) & diag.status == 1 &
                                  diag.time >= 2 & dem.cat == 4, na.rm = TRUE) /
                                  sum(diag.status == 1 & diag.time >= 2 & dem.cat == 4, na.rm = TRUE)
  dat$epi$cc.vsupp.H.m.het[at] <- sum(vl <= log10(200) & diag.status == 1 &
                                  diag.time >= 2 & dem.cat == 5, na.rm = TRUE) /
                                  sum(diag.status == 1 & diag.time >= 2 & dem.cat == 5, na.rm = TRUE)
  dat$epi$cc.vsupp.W.m.het[at] <- sum(vl <= log10(200) & diag.status == 1 &
                                  diag.time >= 2 & dem.cat == 6, na.rm = TRUE) /
                                  sum(diag.status == 1 & diag.time >= 2 & dem.cat == 6, na.rm = TRUE)
  dat$epi$cc.vsupp.B.f.het[at] <- sum(vl <= log10(200) & diag.status == 1 &
                                  diag.time >= 2 & dem.cat == 7, na.rm = TRUE) /
                                  sum(diag.status == 1 & diag.time >= 2 & dem.cat == 7, na.rm = TRUE)
  dat$epi$cc.vsupp.H.f.het[at] <- sum(vl <= log10(200) & diag.status == 1 &
                                  diag.time >= 2 & dem.cat == 8, na.rm = TRUE) /
                                  sum(diag.status == 1 & diag.time >= 2 & dem.cat == 8, na.rm = TRUE)
  dat$epi$cc.vsupp.W.f.het[at] <- sum(vl <= log10(200) & diag.status == 1 &
                                  diag.time >= 2 & dem.cat == 9, na.rm = TRUE) /
                                  sum(diag.status == 1 & diag.time >= 2 & dem.cat == 9, na.rm = TRUE)


  dat$epi$cc.vsupp.all[at] <- sum(vl <= log10(200) & status == 1 & inf.time >= 2, na.rm = TRUE) /
                              sum(status == 1 & inf.time >= 2, na.rm = TRUE)
  dat$epi$cc.vsupp.all.B[at] <- sum(vl <= log10(200) & status == 1 & inf.time >= 2 & race == 1, na.rm = TRUE) /
                                sum(status == 1 & inf.time >= 2 & race == 1, na.rm = TRUE)
  dat$epi$cc.vsupp.all.H[at] <- sum(vl <= log10(200) & status == 1 & inf.time >= 2 & race == 2, na.rm = TRUE) /
                                sum(status == 1 & inf.time >= 2 & race == 2, na.rm = TRUE)
  dat$epi$cc.vsupp.all.W[at] <- sum(vl <= log10(200) & status == 1 & inf.time >= 2 & race == 3, na.rm = TRUE) /
                                sum(status == 1 & inf.time >= 2 & race == 3, na.rm = TRUE)


  dat$epi$cc.vsupp.all.B.msm[at] <- sum(vl <= log10(200) & status == 1 & inf.time >= 2 & dem.cat == 1, na.rm = TRUE) /
                                    sum(status == 1 & inf.time >= 2 & dem.cat == 1, na.rm = TRUE)
  dat$epi$cc.vsupp.all.H.msm[at] <- sum(vl <= log10(200) & status == 1 & inf.time >= 2 & dem.cat == 2, na.rm = TRUE) /
                                    sum(status == 1 & inf.time >= 2 & dem.cat == 2, na.rm = TRUE)
  dat$epi$cc.vsupp.all.W.msm[at] <- sum(vl <= log10(200) & status == 1 & inf.time >= 2 & dem.cat == 3, na.rm = TRUE) /
                                    sum(status == 1 & inf.time >= 2 & dem.cat == 3, na.rm = TRUE)

  dat$epi$cc.vsupp.all.B.m.het[at] <- sum(vl <= log10(200) & status == 1 & inf.time >= 2 & dem.cat == 4, na.rm = TRUE) /
                                      sum(status == 1 & inf.time >= 2 & dem.cat == 4, na.rm = TRUE)
  dat$epi$cc.vsupp.all.H.m.het[at] <- sum(vl <= log10(200) & status == 1 & inf.time >= 2 & dem.cat == 5, na.rm = TRUE) /
                                      sum(status == 1 & inf.time >= 2 & dem.cat == 5, na.rm = TRUE)
  dat$epi$cc.vsupp.all.W.m.het[at] <- sum(vl <= log10(200) & status == 1 & inf.time >= 2 & dem.cat == 6, na.rm = TRUE) /
                                      sum(status == 1 & inf.time >= 2 & dem.cat == 6, na.rm = TRUE)

  dat$epi$cc.vsupp.all.B.f.het[at] <- sum(vl <= log10(200) & status == 1 & inf.time >= 2 & dem.cat == 7, na.rm = TRUE) /
                                      sum(status == 1 & inf.time >= 2 & dem.cat == 7, na.rm = TRUE)
  dat$epi$cc.vsupp.all.H.f.het[at] <- sum(vl <= log10(200) & status == 1 & inf.time >= 2 & dem.cat == 8, na.rm = TRUE) /
                                      sum(status == 1 & inf.time >= 2 & dem.cat == 8, na.rm = TRUE)
  dat$epi$cc.vsupp.all.W.f.het[at] <- sum(vl <= log10(200) & status == 1 & inf.time >= 2 & dem.cat == 9, na.rm = TRUE) /
                                      sum(status == 1 & inf.time >= 2 & dem.cat == 9, na.rm = TRUE)

  dat$epi$cc.vsupp.dur1y[at] <- 1 - (sum((at - vl.last.usupp) <= 52 & diag.time >= 2 &
                                         diag.status == 1, na.rm = TRUE) /
                                        sum(diag.status == 1 & diag.time >= 2, na.rm = TRUE))
  dat$epi$cc.vsupp.dur1y.B[at] <- 1 - (sum((at - vl.last.usupp) <= 52 & diag.time >= 2 &
                                         diag.status == 1 & race == 1, na.rm = TRUE) /
                                        sum(diag.status == 1 & diag.time >= 2 & race == 1, na.rm = TRUE))
  dat$epi$cc.vsupp.dur1y.H[at] <- 1 - (sum((at - vl.last.usupp) <= 52 & diag.time >= 2 &
                                           diag.status == 1 & race == 2, na.rm = TRUE) /
                                       sum(diag.status == 1 & diag.time >= 2 & race == 2, na.rm = TRUE))
  dat$epi$cc.vsupp.dur1y.W[at] <- 1 - (sum((at - vl.last.usupp) <= 52 & diag.time >= 2 &
                                           diag.status == 1 & race == 3, na.rm = TRUE) /
                                       sum(diag.status == 1 & diag.time >= 2 & race == 3, na.rm = TRUE))

  dat$epi$cc.HIV.mr[at] <- (dat$epi$dep.HIV[at]/dat$epi$i.num[at])*52


  #COVERAGE ON ART given infected

  dat$epi$cc.ART.cov.infected[at] <- sum(tx.status = 1, na.rm = TRUE) /
    sum(status == 1, na.rm = TRUE)

  dat$epi$cc.ART.cov.infected.B.msm[at] <- sum(tx.status = 1 & dem.cat == 1, na.rm = TRUE) /
          sum(status == 1 & dem.cat == 1, na.rm = TRUE)
  dat$epi$cc.ART.cov.infected.H.msm[at] <- sum(tx.status = 1 & dem.cat == 2, na.rm = TRUE) /
          sum(status == 1 & dem.cat == 2, na.rm = TRUE)
  dat$epi$cc.ART.cov.infected.W.msm[at] <- sum(tx.status = 1 & dem.cat == 3, na.rm = TRUE) /
          sum(status == 1 & dem.cat == 3, na.rm = TRUE)
  dat$epi$cc.ART.cov.infected.B.m.het[at] <- sum(tx.status = 1 & dem.cat == 4, na.rm = TRUE) /
          sum(status == 1 & dem.cat == 4, na.rm = TRUE)
  dat$epi$cc.ART.cov.infected.H.m.het[at] <- sum(tx.status = 1 & dem.cat == 5, na.rm = TRUE) /
          sum(status == 1 & dem.cat == 5, na.rm = TRUE)
  dat$epi$cc.ART.cov.infected.W.m.het[at] <- sum(tx.status = 1 & dem.cat == 6, na.rm = TRUE) /
          sum(status == 1 & dem.cat == 6, na.rm = TRUE)
  dat$epi$cc.ART.cov.infected.B.f.het[at] <- sum(tx.status = 1 & dem.cat == 7, na.rm = TRUE) /
          sum(status == 1 & dem.cat == 7, na.rm = TRUE)
  dat$epi$cc.ART.cov.infected.H.f.het[at] <- sum(tx.status = 1 & dem.cat == 8, na.rm = TRUE) /
          sum(status == 1 & dem.cat == 8, na.rm = TRUE)
  dat$epi$cc.ART.cov.infected.W.f.het[at] <- sum(tx.status = 1 & dem.cat == 9, na.rm = TRUE) /
          sum(status == 1 & dem.cat == 9, na.rm = TRUE)




  # Care continuum stats (secondary)
  dat$epi$cc.test.int[at] <- mean(at - last.neg.test[diag.status == 0], na.rm = TRUE)
  dat$epi$cc.test.int.B[at] <- mean(at - last.neg.test[diag.status == 0 & race == 1], na.rm = TRUE)
  dat$epi$cc.test.int.H[at] <- mean(at - last.neg.test[diag.status == 0 & race == 2], na.rm = TRUE)
  dat$epi$cc.test.int.W[at] <- mean(at - last.neg.test[diag.status == 0 & race == 3], na.rm = TRUE)

  dat$epi$cc.dx.delay[at] <- mean(diag.time[diag.time >= 2] - inf.time[diag.time >= 2], na.rm = TRUE)
  dat$epi$cc.dx.delay.B[at] <- mean(diag.time[diag.time >= 2 & race == 1] -
                                      inf.time[diag.time >= 2 & race == 1], na.rm = TRUE)
  dat$epi$cc.dx.delay.H[at] <- mean(diag.time[diag.time >= 2 & race == 2] -
                                      inf.time[diag.time >= 2 & race == 2], na.rm = TRUE)
  dat$epi$cc.dx.delay.W[at] <- mean(diag.time[diag.time >= 2 & race == 3] -
                                      inf.time[diag.time >= 2 & race == 3], na.rm = TRUE)


  dat$epi$cc.dx.delay.B.msm[at] <- mean(diag.time[diag.time >= 2 & dem.cat == 1] -
                                      inf.time[diag.time >= 2 & dem.cat == 1], na.rm = TRUE)
  dat$epi$cc.dx.delay.H.msm[at] <- mean(diag.time[diag.time >= 2 & dem.cat == 2] -
                                      inf.time[diag.time >= 2 & dem.cat == 2], na.rm = TRUE)
  dat$epi$cc.dx.delay.W.msm[at] <- mean(diag.time[diag.time >= 2 & dem.cat == 3] -
                                      inf.time[diag.time >= 2 & dem.cat == 3], na.rm = TRUE)
  dat$epi$cc.dx.delay.B.m.het[at] <- mean(diag.time[diag.time >= 2 & dem.cat == 4] -
                                      inf.time[diag.time >= 2 & dem.cat == 4], na.rm = TRUE)
  dat$epi$cc.dx.delay.H.m.het[at] <- mean(diag.time[diag.time >= 2 & dem.cat == 5] -
                                      inf.time[diag.time >= 2 & dem.cat == 5], na.rm = TRUE)
  dat$epi$cc.dx.delay.W.m.het[at] <- mean(diag.time[diag.time >= 2 & dem.cat == 6] -
                                      inf.time[diag.time >= 2 & dem.cat == 6], na.rm = TRUE)
  dat$epi$cc.dx.delay.B.f.het[at] <- mean(diag.time[diag.time >= 2 & dem.cat == 7] -
                                      inf.time[diag.time >= 2 & dem.cat == 7], na.rm = TRUE)
  dat$epi$cc.dx.delay.H.f.het[at] <- mean(diag.time[diag.time >= 2 & dem.cat == 8] -
                                      inf.time[diag.time >= 2 & dem.cat == 8], na.rm = TRUE)
  dat$epi$cc.dx.delay.W.f.het[at] <- mean(diag.time[diag.time >= 2 & dem.cat == 9] -
                                      inf.time[diag.time >= 2 & dem.cat == 9], na.rm = TRUE)


#  dat$epi$cc.dx.delay.int[at] <- mean(diag.time[diag.time >= 3380] - inf.time[diag.time >= 3380], na.rm = TRUE)
#  dat$epi$cc.dx.delay.int.B[at] <- mean(diag.time[diag.time >= 3380 & race == 1] -
#                                      inf.time[diag.time >= 3380 & race == 1], na.rm = TRUE)
#  dat$epi$cc.dx.delay.int.H[at] <- mean(diag.time[diag.time >= 3380 & race == 2] -
#                                      inf.time[diag.time >= 3380 & race == 2], na.rm = TRUE)
#  dat$epi$cc.dx.delay.int.W[at] <- mean(diag.time[diag.time >= 3380 & race == 3] -
#                                      inf.time[diag.time >= 3380 & race == 3], na.rm = TRUE)

  # same as above, but with medians
#  dat$epi$cc.dx.delay.med[at] <- median(diag.time[diag.time >= 2] - inf.time[diag.time >= 2], na.rm = TRUE)
#  dat$epi$cc.dx.delay.B.med[at] <- median(diag.time[diag.time >= 2 & race == 1] -
#                                      inf.time[diag.time >= 2 & race == 1], na.rm = TRUE)
#  dat$epi$cc.dx.delay.H.med[at] <- median(diag.time[diag.time >= 2 & race == 2] -
#                                      inf.time[diag.time >= 2 & race == 2], na.rm = TRUE)
#  dat$epi$cc.dx.delay.W.med[at] <- median(diag.time[diag.time >= 2 & race == 3] -
#                                      inf.time[diag.time >= 2 & race == 3], na.rm = TRUE)
#
#  dat$epi$cc.dx.delay.int.med[at] <- median(diag.time[diag.time >= 3380] - inf.time[diag.time >= 3380], na.rm = TRUE)
#  dat$epi$cc.dx.delay.int.B.med[at] <- median(diag.time[diag.time >= 3380 & race == 1] -
#                                          inf.time[diag.time >= 3380 & race == 1], na.rm = TRUE)
#  dat$epi$cc.dx.delay.int.H.med[at] <- median(diag.time[diag.time >= 3380 & race == 2] -
#                                          inf.time[diag.time >= 3380 & race == 2], na.rm = TRUE)
#  dat$epi$cc.dx.delay.int.W.med[at] <- median(diag.time[diag.time >= 3380 & race == 3] -
#                                          inf.time[diag.time >= 3380 & race == 3], na.rm = TRUE)

  # dat$epi$cc.tx.any1y[at] <- sum((at - dat$attr$tx.period.last <= 52), na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1, na.rm = TRUE)
  # dat$epi$cc.tx.any1y.B[at] <- sum((at - dat$attr$tx.period.last <= 52) & race == 1, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1 & race == 1, na.rm = TRUE)
  # dat$epi$cc.tx.any1y.H[at] <- sum((at - dat$attr$tx.period.last <= 52) & race == 2, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1 & race == 2, na.rm = TRUE)
  # dat$epi$cc.tx.any1y.W[at] <- sum((at - dat$attr$tx.period.last <= 52) & race == 3, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1 & race == 3, na.rm = TRUE)

  # dat$epi$cc.dx.delay[at] <- mean(dat$attr$diag.time - dat$attr$inf.time, na.rm = TRUE)
  # dat$epi$cc.testpy[at] <- 1-sum((at - dat$attr$last.neg.test) > 52 & status == 0,
  #     is.na(dat$attr$last.neg.test) & status == 0, na.rm = TRUE) /
  #   sum(status == 0)
  # dat$epi$cc.linked[at] <- sum(dat$attr$cuml.time.on.tx > 0, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1, na.rm = TRUE)
  # dat$epi$cc.tx[at] <- sum(dat$attr$tx.status == 1, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1, na.rm = TRUE)
  # dat$epi$cc.tx.ret3m[at] <- sum((at - dat$attr$tx.period.last) <= 52 &
  #       (dat$attr$tx.period.last - dat$attr$tx.period.first) > 13, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1, na.rm = TRUE)
  # dat$epi$cc.vsupp.tt1[at] <- sum(dat$attr$vl <= log10(200) & dat$attr$tt.traj == 1, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1 & dat$attr$tt.traj == 1, na.rm = TRUE)
  # dat$epi$cc.vsupp.tt2[at] <- sum(dat$attr$vl <= log10(200) & dat$attr$tt.traj == 2, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1 & dat$attr$tt.traj == 2, na.rm = TRUE)
  # dat$epi$cc.vsupp.tt3[at] <- sum(dat$attr$vl <= log10(200) & dat$attr$tt.traj == 3, na.rm = TRUE) /
  #   sum(dat$attr$diag.status == 1 & dat$attr$tt.traj == 3, na.rm = TRUE)


  # HIV screening outcomes
  dat$epi$mean.neg.tests[at] <- mean(dat$attr$num.neg.tests[diag.status == 0], na.rm = TRUE)
  dat$epi$mean.neg.tests.B[at] <- mean(dat$attr$num.neg.tests[diag.status == 0 & race == 1], na.rm = TRUE)
  dat$epi$mean.neg.tests.H[at] <- mean(dat$attr$num.neg.tests[diag.status == 0 & race == 2], na.rm = TRUE)
  dat$epi$mean.neg.tests.W[at] <- mean(dat$attr$num.neg.tests[diag.status == 0 & race == 3], na.rm = TRUE)

  dat$epi$mean.neg.tests.B.msm[at] <- mean(dat$attr$num.neg.tests[diag.status == 0 & dem.cat == 1], na.rm = TRUE)
  dat$epi$mean.neg.tests.H.msm[at] <- mean(dat$attr$num.neg.tests[diag.status == 0 & dem.cat == 2], na.rm = TRUE)
  dat$epi$mean.neg.tests.W.msm[at] <- mean(dat$attr$num.neg.tests[diag.status == 0 & dem.cat == 3], na.rm = TRUE)
  dat$epi$mean.neg.tests.B.m.het[at] <- mean(dat$attr$num.neg.tests[diag.status == 0 & dem.cat == 4], na.rm = TRUE)
  dat$epi$mean.neg.tests.H.m.het[at] <- mean(dat$attr$num.neg.tests[diag.status == 0 & dem.cat == 5], na.rm = TRUE)
  dat$epi$mean.neg.tests.W.m.het[at] <- mean(dat$attr$num.neg.tests[diag.status == 0 & dem.cat == 6], na.rm = TRUE)
  dat$epi$mean.neg.tests.B.f.het[at] <- mean(dat$attr$num.neg.tests[diag.status == 0 & dem.cat == 7], na.rm = TRUE)
  dat$epi$mean.neg.tests.H.f.het[at] <- mean(dat$attr$num.neg.tests[diag.status == 0 & dem.cat == 8], na.rm = TRUE)
  dat$epi$mean.neg.tests.W.f.het[at] <- mean(dat$attr$num.neg.tests[diag.status == 0 & dem.cat == 9], na.rm = TRUE)

  dat$epi$test.past.year[at] <- sum(at - dat$attr$last.neg.test <= 52 & diag.status == 0, na.rm = TRUE) /
    sum(diag.status == 0, na.rm = TRUE)
  dat$epi$test.past.year.B[at] <- sum(at - dat$attr$last.neg.test <= 52 & diag.status == 0 & race == 1, na.rm = TRUE) /
    sum(diag.status == 0 & race == 1, na.rm = TRUE)
  dat$epi$test.past.year.H[at] <- sum(at - dat$attr$last.neg.test <= 52 & diag.status == 0 & race == 2, na.rm = TRUE) /
    sum(diag.status == 0 & race == 2, na.rm = TRUE)
  dat$epi$test.past.year.W[at] <- sum(at - dat$attr$last.neg.test <= 52 & diag.status == 0 & race == 3, na.rm = TRUE) /
    sum(diag.status == 0 & race == 3, na.rm = TRUE)

  dat$epi$test.past.year.B.msm[at] <- sum(at - dat$attr$last.neg.test <= 52 & diag.status == 0 & dem.cat == 1, na.rm = TRUE) /
                                  sum(diag.status == 0 & dem.cat == 1, na.rm = TRUE)
  dat$epi$test.past.year.H.msm[at] <- sum(at - dat$attr$last.neg.test <= 52 & diag.status == 0 & dem.cat == 2, na.rm = TRUE) /
                                  sum(diag.status == 0 & dem.cat == 2, na.rm = TRUE)
  dat$epi$test.past.year.W.msm[at] <- sum(at - dat$attr$last.neg.test <= 52 & diag.status == 0 & dem.cat == 3, na.rm = TRUE) /
                                  sum(diag.status == 0 & dem.cat == 3, na.rm = TRUE)
  dat$epi$test.past.year.B.m.het[at] <- sum(at - dat$attr$last.neg.test <= 52 & diag.status == 0 & dem.cat == 4, na.rm = TRUE) /
                                  sum(diag.status == 0 & dem.cat == 4, na.rm = TRUE)
  dat$epi$test.past.year.H.m.het[at] <- sum(at - dat$attr$last.neg.test <= 52 & diag.status == 0 & dem.cat == 5, na.rm = TRUE) /
                                  sum(diag.status == 0 & dem.cat == 5, na.rm = TRUE)
  dat$epi$test.past.year.W.m.het[at] <- sum(at - dat$attr$last.neg.test <= 52 & diag.status == 0 & dem.cat == 6, na.rm = TRUE) /
                                  sum(diag.status == 0 & dem.cat == 6, na.rm = TRUE)
  dat$epi$test.past.year.B.f.het[at] <- sum(at - dat$attr$last.neg.test <= 52 & diag.status == 0 & dem.cat == 7, na.rm = TRUE) /
                                  sum(diag.status == 0 & dem.cat == 7, na.rm = TRUE)
  dat$epi$test.past.year.H.f.het[at] <- sum(at - dat$attr$last.neg.test <= 52 & diag.status == 0 & dem.cat == 8, na.rm = TRUE) /
                                  sum(diag.status == 0 & dem.cat == 8, na.rm = TRUE)
  dat$epi$test.past.year.W.f.het[at] <- sum(at - dat$attr$last.neg.test <= 52 & diag.status == 0 & dem.cat == 9, na.rm = TRUE) /
                                  sum(diag.status == 0 & dem.cat == 9, na.rm = TRUE)


  # HIV stage
  dat$epi$hstage.acute[at] <- sum(stage %in% 1:2 & diag.time >= 2, na.rm = TRUE) /
                              sum(status == 1 & diag.time >= 2, na.rm = TRUE)
  dat$epi$hstage.chronic[at] <- sum(stage == 3 & diag.time >= 2, na.rm = TRUE) /
                                sum(status == 1 & diag.time >= 2, na.rm = TRUE)
  dat$epi$hstage.aids[at] <- sum(stage == 4 & diag.time >= 2, na.rm = TRUE) /
                             sum(status == 1 & diag.time >= 2, na.rm = TRUE)

  dat$epi$prepElig[at] <- sum(prepElig == 1, na.rm = TRUE)
  dat$epi$prepElig.B[at] <- sum(prepElig == 1 & race == 1, na.rm = TRUE)
  dat$epi$prepElig.H[at] <- sum(prepElig == 1 & race == 2, na.rm = TRUE)
  dat$epi$prepElig.W[at] <- sum(prepElig == 1 & race == 3, na.rm = TRUE)

  dat$epi$prepElig.B.msm[at] <- sum(prepElig == 1 & dem.cat == 1, na.rm = TRUE)
  dat$epi$prepElig.H.msm[at] <- sum(prepElig == 1 & dem.cat == 2, na.rm = TRUE)
  dat$epi$prepElig.W.msm[at] <- sum(prepElig == 1 & dem.cat == 3, na.rm = TRUE)
  dat$epi$prepElig.B.m.het[at] <- sum(prepElig == 1 & dem.cat == 4, na.rm = TRUE)
  dat$epi$prepElig.H.m.het[at] <- sum(prepElig == 1 & dem.cat == 5, na.rm = TRUE)
  dat$epi$prepElig.W.m.het[at] <- sum(prepElig == 1 & dem.cat == 6, na.rm = TRUE)
  dat$epi$prepElig.B.f.het[at] <- sum(prepElig == 1 & dem.cat == 7, na.rm = TRUE)
  dat$epi$prepElig.H.f.het[at] <- sum(prepElig == 1 & dem.cat == 8, na.rm = TRUE)
  dat$epi$prepElig.W.f.het[at] <- sum(prepElig == 1 & dem.cat == 9, na.rm = TRUE)

  dat$epi$prepCurr[at] <- sum(prepStat == 1, na.rm = TRUE)
  dat$epi$prepCurr.B[at] <- sum(prepStat == 1 & race == 1, na.rm = TRUE)
  dat$epi$prepCurr.H[at] <- sum(prepStat == 1 & race == 2, na.rm = TRUE)
  dat$epi$prepCurr.W[at] <- sum(prepStat == 1 & race == 3, na.rm = TRUE)

  dat$epi$prepCurr.B.msm[at] <- sum(prepStat == 1 & dem.cat == 1, na.rm = TRUE)
  dat$epi$prepCurr.H.msm[at] <- sum(prepStat == 1 & dem.cat == 2, na.rm = TRUE)
  dat$epi$prepCurr.W.msm[at] <- sum(prepStat == 1 & dem.cat == 3, na.rm = TRUE)
  dat$epi$prepCurr.B.m.het[at] <- sum(prepStat == 1 & dem.cat == 4, na.rm = TRUE)
  dat$epi$prepCurr.H.m.het[at] <- sum(prepStat == 1 & dem.cat == 5, na.rm = TRUE)
  dat$epi$prepCurr.W.m.het[at] <- sum(prepStat == 1 & dem.cat == 6, na.rm = TRUE)
  dat$epi$prepCurr.B.f.het[at] <- sum(prepStat == 1 & dem.cat == 7, na.rm = TRUE)
  dat$epi$prepCurr.H.f.het[at] <- sum(prepStat == 1 & dem.cat == 8, na.rm = TRUE)
  dat$epi$prepCurr.W.f.het[at] <- sum(prepStat == 1 & dem.cat == 9, na.rm = TRUE)

  dat$epi$prepCurr.B.msm[at] <- max(0,(dat$epi$prepCurr.B.msm[at] / dat$epi$prepElig.B.msm[at]) )
  dat$epi$prepCurr.H.msm[at] <- max(0,(dat$epi$prepCurr.H.msm[at] / dat$epi$prepElig.H.msm[at]) )
  dat$epi$prepCurr.W.msm[at] <- max(0,(dat$epi$prepCurr.W.msm[at] / dat$epi$prepElig.W.msm[at]) )
  dat$epi$prepCurr.B.m.het[at] <- max(0,(dat$epi$prepCurr.B.m.het[at] / dat$epi$prepElig.B.m.het[at]) )
  dat$epi$prepCurr.H.m.het[at] <- max(0,(dat$epi$prepCurr.H.m.het[at] / dat$epi$prepElig.H.m.het[at]) )
  dat$epi$prepCurr.W.m.het[at] <- max(0,(dat$epi$prepCurr.W.m.het[at] / dat$epi$prepElig.W.m.het[at]) )
  dat$epi$prepCurr.B.f.het[at] <- max(0,(dat$epi$prepCurr.B.f.het[at] / dat$epi$prepElig.B.f.het[at]) )
  dat$epi$prepCurr.H.f.het[at] <- max(0,(dat$epi$prepCurr.H.f.het[at] / dat$epi$prepElig.H.f.het[at]) )
  dat$epi$prepCurr.W.f.het[at] <- max(0,(dat$epi$prepCurr.W.f.het[at] / dat$epi$prepElig.W.f.het[at]) )




  dat$epi$prepCurr.hadr[at] <- sum(prepStat == 1 & prepClass == 3, na.rm = TRUE)

  # STIs
#  dat$epi$prev.gc[at] <- sum((rGC == 1 | uGC == 1), na.rm = TRUE) / dat$epi$num[at]
#  dat$epi$prev.ct[at] <- sum((rCT == 1 | uCT == 1), na.rm = TRUE) / dat$epi$num[at]
#  ir100.rgc <- (dat$epi$incid.rgc[at]/sum(rGC == 0, dat$epi$incid.rgc[at], na.rm = TRUE))*5200
#  ir100.ugc <- (dat$epi$incid.ugc[at]/sum(uGC == 0, dat$epi$incid.ugc[at], na.rm = TRUE))*5200
#  dat$epi$ir100.gc[at] <- ir100.rgc + ir100.ugc
#  ir100.rct <- (dat$epi$incid.rct[at]/sum(rCT == 0, dat$epi$incid.rct[at], na.rm = TRUE))*5200
#  ir100.uct <- (dat$epi$incid.uct[at]/sum(uCT == 0, dat$epi$incid.uct[at], na.rm = TRUE))*5200
#  dat$epi$ir100.ct[at] <- ir100.rct + ir100.uct
#  dat$epi$ir100.sti[at] <- dat$epi$ir100.gc[at] + dat$epi$ir100.ct[at]

  return(dat)
}

#' @export
#' @rdname prevalence_msm
prevalence_het <- function(dat, at) {

  status <- dat$attr$status
  male <- dat$attr$male
  age <- dat$attr$age

  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)

  # Initialize vectors
  if (at == 1) {
    dat$epi$i.num <- rNA
    dat$epi$num <- rNA

    dat$epi$i.num.male <- rNA
    dat$epi$i.num.feml <- rNA
    dat$epi$i.prev.male <- rNA
    dat$epi$i.prev.feml <- rNA

    dat$epi$num.male <- rNA
    dat$epi$num.feml <- rNA
    dat$epi$meanAge <- rNA
    dat$epi$propMale <- rNA

    dat$epi$si.flow <- rNA
    dat$epi$si.flow.male <- rNA
    dat$epi$si.flow.feml <- rNA

    dat$epi$b.flow <- rNA
    dat$epi$ds.flow <- dat$epi$di.flow <- rNA
  }

  dat$epi$i.num[at] <- sum(status == 1, na.rm = TRUE)
  dat$epi$num[at] <- length(status)

  dat$epi$i.num.male[at] <- sum(status == 1 & male == 1, na.rm = TRUE)
  dat$epi$i.num.feml[at] <- sum(status == 1 & male == 0, na.rm = TRUE)
  dat$epi$i.prev.male[at] <- sum(status == 1 & male == 1, na.rm = TRUE) /
    sum(male == 1, na.rm = TRUE)
  dat$epi$i.prev.feml[at] <- sum(status == 1 & male == 0, na.rm = TRUE) /
    sum(male == 0, na.rm = TRUE)

  dat$epi$num.male[at] <- sum(male == 1, na.rm = TRUE)
  dat$epi$num.feml[at] <- sum(male == 0, na.rm = TRUE)
  dat$epi$meanAge[at] <- mean(age, na.rm = TRUE)
  dat$epi$propMale[at] <- mean(male, na.rm = TRUE)

  return(dat)
}


whichVlSupp <- function(attr, param) {
  which(attr$status == 1 &
        attr$vlLevel <= log10(50) &
        (attr$age - attr$ageInf) * (365 / param$time.unit) >
        (param$vl.acute.topeak + param$vl.acute.toset))
}
