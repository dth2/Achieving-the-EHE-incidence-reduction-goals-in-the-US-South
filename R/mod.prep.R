
#' @title PrEP Module
#'
#' @description Module function for implementation and uptake of pre-exposure
#'              prophylaxis (PrEP) to prevent HIV infection.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
prep_msm <- function(dat, at) {


  # Function Selection ------------------------------------------------------

  if (at >= dat$param$riskh.start) {
    dat <- riskhist_msm(dat, at)
  } else {
    return(dat)
  }

  if (at < dat$param$prep.start) {
    return(dat)
  }

  # Attributes --------------------------------------------------------------

  # Core Attributes
  active <- dat$attr$active
  status <- dat$attr$status
  diag.status <- dat$attr$diag.status
  lnt <- dat$attr$last.neg.test
  dem.cat <- dat$attr$dem.cat
  msm <- dat$attr$msm
  age <- dat$attr$age

  # PrEP Attributes
  prepElig <- dat$attr$prepElig
  prepStat <- dat$attr$prepStat
  prepClass <- dat$attr$prepClass
  prepLastRisk <- dat$attr$prepLastRisk
  prepStartTime <- dat$attr$prepStartTime
  prepLastStiScreen <- dat$attr$prepLastStiScreen


  # Parameters --------------------------------------------------------------

  prep.start.prob <- dat$param$prep.start.prob
  prep.adhr.dist.msm <- dat$param$prep.adhr.dist.msm
  prep.adhr.dist.het <- dat$param$prep.adhr.dist.het
  prep.require.lnt <- dat$param$prep.require.lnt
  prep.risk.reassess.method <- dat$param$prep.risk.reassess.method
  prep.discont.rate <- dat$param$prep.discont.rate

  prep.base.cov <- dat$param$prep.base.cov
  dem.cat.prep.fixed <- c(dat$param$prep.int.1,dat$param$prep.int.2,dat$param$prep.int.3,
                          dat$param$prep.int.4,dat$param$prep.int.5,dat$param$prep.int.6,
                          dat$param$prep.int.7,dat$param$prep.int.8,dat$param$prep.int.9)
  dem.cat.prep.fixed.prop <- dem.cat.prep.fixed * dat$param$prep.int.cov

#  prep.start.prob = 0.2,
#  prep.adhr.dist = c(0.089, 0.127, 0.784),
#  prep.adhr.hr = c(0.69, 0.19, 0.01),
#  prep.discont.rate = 1 - (2^(-1/(224.4237/7))),
#  prep.tst.int = 90/7,
#  prep.risk.int = 182/7,
#  prep.sti.screen.int = 182/7,
#  prep.sti.prob.tx = 1,
#  prep.risk.reassess.method = "year",
#  prep.fixed = TRUE,
  #after prep.time 1 selection is from ART naive to fix the coverage
  #after prep.time 2 selection is from ART experienced to fix the coverage
#  fixed.prep.time = c(1,4),
#  dem.cat.prep.fixed = c(0,0,0,0,0,0,1,1,0),
#  dem.cat.prep.fixed.prop  = c(.95,.95,.95,.95,.95,.95,.95,.95,.95),


  # Indications -------------------------------------------------------------

  ind1 <- dat$attr$prep.ind.ui.mono
  # ind2 <- dat$attr$prep.ind.uai.nmain
  ind2 <- dat$attr$prep.ind.ui.conc
  ind3 <- dat$attr$prep.ind.sti



  twind <- at - dat$param$prep.risk.int

  # No indications in window
  idsNoIndic <- which((ind1 < twind | is.na(ind1)) &
                      (ind2 < twind | is.na(ind2)) &
                      (ind3 < twind | is.na(ind3)))
  base.cond.no <- which(active == 0 | diag.status == 1)
  idsNoIndic <- union(idsNoIndic, base.cond.no)

  # Indications in window
  idsIndic <- which(ind1 >= twind | ind2 >= twind | ind3 >= twind)
  base.cond.yes <- which(active == 1 & diag.status == 0)
  idsIndic <- intersect(idsIndic, base.cond.yes)

  # Set eligibility to 1 if indications
  prepElig[idsIndic] <- 1

  # Set eligibility to 0 if no indications
  prepElig[idsNoIndic] <- 0


  ## Stoppage ------------------------------------------------------------------

  # Indication lapse
  # Rules = None, instant, yearly (CDC guidelines)
  if (prep.risk.reassess.method == "none") {
    idsStpInd <- NULL
  } else if (prep.risk.reassess.method == "inst") {
    idsRiskAssess <- which(active == 1 & prepStat == 1)
    prepLastRisk[idsRiskAssess] <- at
    idsStpInd <- intersect(idsNoIndic, idsRiskAssess)
  } else if (prep.risk.reassess.method == "year") {
    idsRiskAssess <- which(active == 1 &
                           prepStat == 1 &
                           lnt == at &
                           (at - prepLastRisk) >= 52)
    prepLastRisk[idsRiskAssess] <- at
    idsStpInd <- intersect(idsNoIndic, idsRiskAssess)
  }

  # Random discontinuation
  #Make this dem.group specific
  idsStpRand <- NULL
  for(i in 1:9){
  idsEligStpRand <- which(active == 1 & prepStat == 1 & dem.cat == i)
  vecStpRand <- rbinom(length(idsEligStpRand), 1, prep.discont.rate[i])
  idsStpRand.temp <- idsEligStpRand[which(vecStpRand == 1)]
  idsStpRand <- c(idsStpRand, idsStpRand.temp)
  }

  # Diagnosis
  idsStpDx <- which(active == 1 & prepStat == 1 & diag.status == 1)

  # Death
  idsStpDth <- which(active == 0 & prepStat == 1)

  # Reset PrEP status
  idsStp <- c(idsStpInd, idsStpRand, idsStpDx, idsStpDth)

  # Update attributes for stoppers
  prepStat[idsStp] <- 0
  prepLastRisk[idsStp] <- NA
  prepStartTime[idsStp] <- NA
  prepLastStiScreen[idsStp] <- NA


  ## Initiation ----------------------------------------------------------------

  ## Eligibility ##

  # Indications to start
  if (prep.require.lnt == TRUE) {
    idsEligStart <- which(prepStat == 0 & lnt == at)
  } else {
    idsEligStart <- which(prepStat == 0 & status == 0)
  }

  idsEligStart <- intersect(idsIndic, idsEligStart)
  prepElig[idsEligStart] <- 1


  #Baseline
  idsStart <-NULL

    for(i in 1:9){
      #prep.base.cov
      #get ids.ELIG.start with demcat i
      dem.list<- dem.cat[idsEligStart]
      cov <- prep.base.cov[i]

      if(i < 4){
      count <- (sum(dat$attr$msm==1)/100000) * cov
      on.prep <- which(dat$attr$dem.cat==i & dat$attr$prepStat == 1)
      count <- count - length(on.prep)
      }

      if(i > 3){
        count <- (sum(dat$attr$msm==0)/100000) * cov
        on.prep <- which(dat$attr$dem.cat==i & dat$attr$prepStat == 1)
        count <- count - length(on.prep)
      }

      if (count > 0){
      idsStart.temp <- sample(idsEligStart[dem.list==i],count,FALSE)
      idsStart <- c(idsStart,idsStart.temp)
      }

}
  ## IF PREP COVERAGE IS GOING TO BE FIXED FOR EFFECT SIZE ANALYSIS

  if(dat$param$prep.fixed==TRUE){

for(i in 1:9){
  if(dem.cat.prep.fixed[i] == 1){
    dem.list<- dem.cat[idsEligStart]
    cov <- dem.cat.prep.fixed.prop[i]

    count <- sum(dem.cat==i &  prepElig == 1) * cov
    #subtract of those on PrEP
    on.prep <- which(dat$attr$dem.cat==i & dat$attr$prepStat == 1)
    count <- count - length(on.prep)

    idsStart.temp <- NULL
    if(count > 0){
    idsStart.temp <- sample(idsEligStart[dem.list==i],count,FALSE)
    }

    if(count > 0){
    idsStart <- c(idsStart,idsStart.temp)
    idsStart <- unique(idsStart)
    }

}
    }
  }


  # Set attributes for starters
  if (length(idsStart) > 0) {
    prepStat[idsStart] <- 1
    prepStartTime[idsStart] <- at
    prepLastRisk[idsStart] <- at

    # PrEP adherence class
    # Split out by MSM/HET

    # HET
    needPC <- which(is.na(prepClass[idsStart]) & msm[idsStart]==0)
    prepClass[idsStart[needPC]] <- sample(x = 1:3, size = length(needPC),
                                          replace = TRUE, prob = prep.adhr.dist.het)

    # MSM
    needPC <- which(is.na(prepClass[idsStart]) & msm[idsStart]==1)
    prepClass[idsStart[needPC]] <- sample(x = 1:3, size = length(needPC),
                                          replace = TRUE, prob = prep.adhr.dist.msm)



    }


  ## Output --------------------------------------------------------------------

  # Random discontinuation
  dat$epi$prep.rand.stop[at] <- length(idsStpRand)

  # Attributes
  dat$attr$prepElig <- prepElig
  dat$attr$prepStat <- prepStat
  dat$attr$prepClass <- prepClass

  dat$attr$prepStartTime <- prepStartTime
  dat$attr$prepLastRisk <- prepLastRisk
  dat$attr$prepLastStiScreen <- prepLastStiScreen

  return(dat)
}


#' @title Risk History Sub-Module
#'
#' @description Sub-Module function to track the risk history of uninfected persons
#'              for purpose of PrEP targeting.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
riskhist_msm <- function(dat, at) {

  ## Attributes
  n <- length(dat$attr$active)
  dx <- dat$attr$diag.status
  since.test <- at - dat$attr$last.neg.test
  rGC.tx <- dat$attr$rGC.tx
  uGC.tx <- dat$attr$uGC.tx
  rCT.tx <- dat$attr$rCT.tx
  uCT.tx <- dat$attr$uCT.tx

  ## Parameters

  ## Edgelist, adds ui summation per partnership from act list
  #try sorting al and el by p1*1e7 + p2
  pid <- NULL # For R CMD Check
  al <- as.data.frame(dat$temp$al)
  by_pid <- group_by(al, pid)
  ui <- summarise(by_pid, ui = sum(ui))[, 2]
  el <- as.data.frame(cbind(dat$temp$el, ui))

  if (max(el[, 1:2]) > n) stop("riskhist max(el) > n")

  # Remove concordant positive edges
  el2 <- el[el$st2 == 0, ]

  # Initialize attributes
  if (is.null(dat$attr$prep.ind.ui.mono)) {
    dat$attr$prep.ind.ui.mono <- rep(NA, n)
    dat$attr$prep.ind.ui.nmain <- rep(NA, n)
    dat$attr$prep.ind.sti <- rep(NA, n)
  }
  if (is.null(dat$attr$prep.ind.ui.conc)) {
    dat$attr$prep.ind.ui.conc <- rep(NA, n)
  }

  ## Degree ##
  main.deg.het <- get_degree(dat$el[[1]])
  casl.deg.het <- get_degree(dat$el[[2]])
  inst.deg.het <- get_degree(dat$el[[3]])

  main.deg.msm <- get_degree(dat$el[[4]])
  casl.deg.msm <- get_degree(dat$el[[5]])
  inst.deg.msm <- get_degree(dat$el[[6]])

  ## Preconditions ##

  # Any UI
  ui.any <- unique(c(el2$p1[el2$ui > 0],
                      el2$p2[el2$ui > 0]))

  # Monogamous partnerships: 1-sided
  tot.deg <- main.deg.het + casl.deg.het + inst.deg.het + main.deg.msm + casl.deg.msm + inst.deg.msm
  ui.mono1 <- intersect(which(tot.deg == 1), ui.any)

  # "Negative" partnerships
  tneg <- unique(c(el2$p1[el2$st1 == 0], el2$p2[el2$st1 == 0]))
  fneg <- unique(c(el2$p1[which(dx[el2$p1] == 0)], el2$p2[which(dx[el2$p1] == 0)]))
  all.neg <- c(tneg, fneg)

  ## Condition 1b: UI in 1-sided "monogamous" "negative" partnership,
  ##               partner not tested in past 6 months
  ui.mono1.neg <- intersect(ui.mono1, all.neg)
  part.id1 <- c(el2[el2$p1 %in% ui.mono1.neg, 2], el2[el2$p2 %in% ui.mono1.neg, 1])
  not.tested.6mo <- since.test[part.id1] > (180/7)
  part.not.tested.6mo <- ui.mono1.neg[which(not.tested.6mo == TRUE)]
  dat$attr$prep.ind.ui.mono[part.not.tested.6mo] <- at

  ## Condition 2a: UI + concurrency
  el2.ui <- el2[el2$ui > 0, ]
  vec <- c(el2.ui[, 1], el2.ui[, 2])
  ui.conc <- unique(vec[duplicated(vec)])
  dat$attr$prep.ind.ui.conc[ui.conc] <- at

  ## Condition 2b: UAI in non-main partnerships
  ui.nmain.het <- unique(c(el2$p1[el2$st1 == 0 & el2$ui > 0 & el2$ptype %in% 2:3],
                        el2$p2[el2$ui > 0 & el2$ptype %in% 2:3]))

  ui.nmain.msm <- unique(c(el2$p1[el2$st1 == 0 & el2$ui > 0 & el2$ptype %in% 5:6],
                           el2$p2[el2$ui > 0 & el2$ptype %in% 5:6]))

  ui.nmain <- c(ui.nmain.het, ui.nmain.msm)

  dat$attr$prep.ind.ui.nmain[ui.nmain] <- at

  ## Condition 4, any STI diagnosis
  idsDx <- which(rGC.tx == 1 | uGC.tx == 1 | rCT.tx == 1 | uCT.tx == 1)
  dat$attr$prep.ind.sti[idsDx] <- at

  return(dat)
}


#' @title Proportionally Reallocate PrEP Adherence Class Probability
#'
#' @description Shifts probabilities from the high-adherence category to the lower
#'              three adherence categories while maintaining the proportional
#'              distribution of those lower categories.
#'
#' @param in.pcp Input vector of length four for the \code{prep.adhr.dist}
#'        parameters.
#' @param reall The pure percentage points to shift from the high adherence
#'        group to the lower three groups.
#'
#' @export
#'
reallocate_pcp <- function(in.pcp = c(0.089, 0.127, 0.784), reall = 0) {

  dist <- in.pcp[1]/sum(in.pcp[1:2])
  dist[2] <- in.pcp[2]/sum(in.pcp[1:2])

  out.pcp <- rep(NA, 3)
  out.pcp[1:2] <- in.pcp[1:2] - (dist * reall)
  out.pcp[3] <- 1 - sum(out.pcp[1:2])

  return(out.pcp)
}
