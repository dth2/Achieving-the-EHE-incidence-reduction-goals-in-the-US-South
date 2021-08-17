

# MSM -----------------------------------------------------------------

#' @title Network Resimulation Module
#'
#' @description Module function for resimulating the sexual networks for one
#'              time step.
#'
#' @inheritParams aging_msm
#'
#' @keywords module msm
#'
#' @export
#'
simnet_msm <- function(dat, at) {



  ## Edges correction
  dat <- edges_correct_msm(dat, at)

  ## Main heterosexual network
  nwparam.m <- EpiModel::get_nwparam(dat, network = 1)

  dat$attr$deg.casl.het <- get_degree(dat$el[[2]])
  dat$attr$deg.casl.c.het <- ifelse(dat$attr$deg.casl.het > 1, 1,dat$attr$deg.casl.het)

   dat <- tergmLite::updateModelTermInputs(dat, network = 1)

  rv_1 <- tergmLite::simulate_network(state = dat$p[[1]]$state,
                                      coef = c(nwparam.m$coef.form, nwparam.m$coef.diss$coef.adj),
                                      control = dat$control$mcmc.control[[1]],
                                      save.changes = TRUE)

  dat$el[[1]] <- rv_1$el

  if(dat$control$tergmLite.track.duration) {
    dat$p[[1]]$state$nw0 %n% "time" <- rv_1$state$nw0 %n% "time"
    dat$p[[1]]$state$nw0 %n% "lasttoggle" <- rv_1$state$nw0 %n% "lasttoggle"
  }

  plist1 <- update_plist(dat, at, ptype = 1)


  ## Main MSM network
  nwparam.m <- EpiModel::get_nwparam(dat, network = 4)

  dat$attr$deg.casl.msm <- get_degree(dat$el[[5]])
  dat$attr$deg.casl.c.msm <- ifelse(dat$attr$deg.casl.msm > 1, 1, dat$attr$deg.casl.msm)

  dat <- tergmLite::updateModelTermInputs(dat, network = 4)

  rv_4 <- tergmLite::simulate_network(state = dat$p[[4]]$state,
                                      coef = c(nwparam.m$coef.form, nwparam.m$coef.diss$coef.adj),
                                      control = dat$control$mcmc.control[[4]],
                                      save.changes = TRUE)
  dat$el[[4]] <- rv_4$el

  if(dat$control$tergmLite.track.duration) {
    dat$p[[4]]$state$nw0 %n% "time" <- rv_4$state$nw0 %n% "time"
    dat$p[[4]]$state$nw0 %n% "lasttoggle" <- rv_4$state$nw0 %n% "lasttoggle"
  }

  plist4 <- update_plist(dat, at, ptype = 4)


  ## Casual heterosexual network

  nwparam.p <- EpiModel::get_nwparam(dat, network = 2)

  dat$attr$deg.main.het <- get_degree(dat$el[[1]])
  dat$attr$deg.main.c.het <- ifelse(dat$attr$deg.main.het > 1, 1, dat$attr$deg.main.het)

  dat <- tergmLite::updateModelTermInputs(dat, network = 2)

  rv_2 <- tergmLite::simulate_network(state = dat$p[[2]]$state,
                                      coef = c(nwparam.p$coef.form, nwparam.p$coef.diss$coef.adj),
                                      control = dat$control$mcmc.control[[2]],
                                      save.changes = TRUE)
  dat$el[[2]] <- rv_2$el

  dat$attr$deg.casl.het <- get_degree(dat$el[[2]])
  dat$attr$deg.casl.c.het <- ifelse(dat$attr$deg.casl.het > 1, 1, dat$attr$deg.casl.het)

  if(dat$control$tergmLite.track.duration) {
    dat$p[[2]]$state$nw0 %n% "time" <- rv_2$state$nw0 %n% "time"
    dat$p[[2]]$state$nw0 %n% "lasttoggle" <- rv_2$state$nw0 %n% "lasttoggle"
  }

  plist2 <- update_plist(dat, at, ptype = 2)



  ## Casual msm network

  nwparam.p <- EpiModel::get_nwparam(dat, network = 5)

  dat$attr$deg.main.msm <- get_degree(dat$el[[4]])
  dat$attr$deg.main.c.msm <- ifelse(dat$attr$deg.main.msm > 1, 1, dat$attr$deg.main.msm)

  dat <- tergmLite::updateModelTermInputs(dat, network = 5)

  rv_5 <- tergmLite::simulate_network(state = dat$p[[5]]$state,
                                      coef = c(nwparam.p$coef.form, nwparam.p$coef.diss$coef.adj),
                                      control = dat$control$mcmc.control[[5]],
                                      save.changes = TRUE)
  dat$el[[5]] <- rv_5$el

  dat$attr$deg.casl.msm <- get_degree(dat$el[[5]])
  dat$attr$deg.casl.c.msm <- ifelse(dat$attr$deg.casl.msm > 1, 1, dat$attr$deg.casl.msm)

  if(dat$control$tergmLite.track.duration) {
    dat$p[[5]]$state$nw0 %n% "time" <- rv_5$state$nw0 %n% "time"
    dat$p[[5]]$state$nw0 %n% "lasttoggle" <- rv_5$state$nw0 %n% "lasttoggle"
  }

  plist5 <- update_plist(dat, at, ptype = 5)

  dat$temp$plist <- rbind(plist1, plist2, plist4, plist5)
  if (dat$control$truncate.plist == TRUE) {
    to.keep <- which(is.na(dat$temp$plist[, "stop"]))
    dat$temp$plist <- dat$temp$plist[to.keep, ]
  }


  ## One-off heterosexual network
  nwparam.i <- EpiModel::get_nwparam(dat, network = 3)

  dat$attr$deg.tot.het <- pmin(dat$attr$deg.main.het + get_degree(dat$el[[2]]), 3)
  dat$attr$deg.tot.c.het <- ifelse(dat$attr$deg.tot.het > 2, 2, dat$attr$deg.tot.het)

  dat <- tergmLite::updateModelTermInputs(dat, network = 3)

  rv_3 <- tergmLite::simulate_ergm(state = dat$p[[3]]$state,
                                   coef = nwparam.i$coef.form,
                                   control = dat$control$mcmc.control[[3]])
  dat$el[[3]] <- rv_3$el


  ## One-off MSM network
  nwparam.i <- EpiModel::get_nwparam(dat, network = 6)

  dat$attr$deg.tot.msm <- pmin(dat$attr$deg.main.msm + get_degree(dat$el[[5]]), 3)
  dat$attr$deg.tot.c.msm <- ifelse(dat$attr$deg.tot.msm > 2, 2,  dat$attr$deg.tot.msm)

  dat <- tergmLite::updateModelTermInputs(dat, network = 6)

  rv_6 <- tergmLite::simulate_ergm(state = dat$p[[6]]$state,
                                   coef = nwparam.i$coef.form,
                                   control = dat$control$mcmc.control[[6]])
  dat$el[[6]] <- rv_6$el

  if (dat$control$save.nwstats == TRUE) {
    for (i in 1:6) {
      nwL <- networkLite(dat$el[[i]], dat$attr)
      if (dat$control$tergmLite.track.duration) {
        nwL %n% "time" <- dat$p[[i]]$state$nw0 %n% "time"
        nwL %n% "lasttoggle" <- dat$p[[i]]$state$nw0 %n% "lasttoggle"
      }
      dat$stats$nwstats[[i]] <- rbind(dat$stats$nwstats[[i]],
                                      summary(dat$control$nwstats.formulas[[i]],
                                              basis = nwL,
                                              term.options = dat$control$mcmc.control[[i]]$term.options))
    }
  }


  return(dat)
}

# updates the partnership list
update_plist <- function(dat, at, ptype) {
  # pull existing partner type specific list
  plist1 <- dat$temp$plist[dat$temp$plist[, "ptype"] == ptype, ]

  # look up dissolutions, update stop time
  uid <- dat$attr$uid
  news <- attr(dat$el[[ptype]], "changes")
  news_uid <- cbind(matrix(uid[news[, 1:2]], ncol = 2), news[, 3])
  news_uid_stop <- news_uid[news_uid[, 3] == 0, , drop = FALSE]
  pid_plist1 <- plist1[, 1]*1e7 + plist1[, 2]
  pid_stop <- news_uid_stop[, 1]*1e7 + news_uid_stop[, 2]
  matches_stop <- match(pid_stop, pid_plist1)
  plist1[matches_stop, "stop"] <- at

  # look up new formations, row bind them
  news_uid_start <- news_uid[news_uid[, 3] == 1, , drop = FALSE]
  plist1 <- rbind(plist1, cbind(news_uid_start[, 1:2, drop = FALSE], ptype, at, NA))

  return(plist1)
}


calc_nwstats <- function(dat, at) {

  for (nw in 1:6) {
    n <- attr(dat$el[[nw]], "n")
    edges <- nrow(dat$el[[nw]])
    meandeg <- round(edges * (2/n), 3)
    concurrent <- round(mean(get_degree(dat$el[[nw]]) > 1), 3)
    mat <- matrix(c(edges, meandeg, concurrent), ncol = 3, nrow = 1)
    if (at == 2) {
      dat$stats$nwstats[[nw]] <- mat
      colnames(dat$stats$nwstats[[nw]]) <- c("edges", "mdeg", "conc")
    }
    if (at > 2) {
      dat$stats$nwstats[[nw]] <- rbind(dat$stats$nwstats[[nw]], mat)
    }
  }

  return(dat)
}



#' @title Adjustment for the Edges Coefficient with Changing Network Size
#'
#' @description Adjusts the edges coefficients in a dynamic network model
#'              to preserve the mean degree.
#'
#' @inheritParams aging_msm
#'
#' @details
#' In HIV/STI modeling, there is typically an assumption that changes in
#' population size do not affect one's number of partners, specified as the
#' mean degree for network models. A person would not have 10 times the number
#' of partners should he move from a city 10 times as large. This module uses
#' the adjustment of Krivitsky et al. to adjust the edges coefficients on the
#' three network models to account for varying population size in order to
#' preserve that mean degree.
#'
#' @return
#' The network model parameters stored in \code{dat$nwparam} are updated for
#' each of the three network models.
#'
#' @references
#' Krivitsky PN, Handcock MS, and Morris M. "Adjusting for network size and
#' composition effects in exponential-family random graph models." Statistical
#' Methodology. 2011; 8.4: 319-339.
#'
#' @keywords module msm
#'
#' @export
#'
edges_correct_msm <- function(dat, at) {

  old.num <- dat$epi$num[at - 1]
  new.num <- sum(dat$attr$active == 1, na.rm = TRUE)
  adjust <- log(old.num) - log(new.num)

  coef.form.m.het <- get_nwparam(dat, network = 1)$coef.form
  coef.form.m.het[1] <- coef.form.m.het[1] + adjust
  dat$nwparam[[1]]$coef.form <- coef.form.m.het

  coef.form.m.msm <- get_nwparam(dat, network = 4)$coef.form
  coef.form.m.msm[1] <- coef.form.m.msm[1] + adjust
  dat$nwparam[[4]]$coef.form <- coef.form.m.msm

  coef.form.p.het <- get_nwparam(dat, network = 2)$coef.form
  coef.form.p.het[1] <- coef.form.p.het[1] + adjust
  dat$nwparam[[2]]$coef.form <- coef.form.p.het

  coef.form.p.msm <- get_nwparam(dat, network = 5)$coef.form
  coef.form.p.msm[1] <- coef.form.p.msm[1] + adjust
  dat$nwparam[[5]]$coef.form <- coef.form.p.msm

  coef.form.i.het <- get_nwparam(dat, network = 3)$coef.form
  coef.form.i.het[1] <- coef.form.i.het[1] + adjust
  dat$nwparam[[3]]$coef.form <- coef.form.i.het

  coef.form.i.msm <- get_nwparam(dat, network = 6)$coef.form
  coef.form.i.msm[1] <- coef.form.i.msm[1] + adjust
  dat$nwparam[[6]]$coef.form <- coef.form.i.msm

  return(dat)
}




# HET -----------------------------------------------------------------


#' @export
#' @rdname simnet_msm
simnet_het <- function(dat, at) {

  # Update edges coefficients
  dat <- edges_correct_het(dat, at)

  # Update internal ergm data
  dat <- tergmLite::updateModelTermInputs(dat, network = 1)

  # Pull network parameters
  nwparam <- get_nwparam(dat, network = 1)

  # Simulate edgelist
  rv_1 <- tergmLite::simulate_network(state = dat$p[[1]]$state,
                                      coef = c(nwparam$coef.form, nwparam$coef.diss$coef.adj),
                                      control = dat$control$mcmc.control[[1]])
  dat$el[[1]] <- rv_1$el
  dat$p[[1]]$state$el <- rv_1$state$el

  return(dat)
}


#' @export
#' @rdname edges_correct_msm
edges_correct_het <- function(dat, at) {

  # Popsize
  old.num <- dat$epi$num[at - 1]
  new.num <- sum(dat$attr$active == 1, na.rm = TRUE)

  # New Coefs
  coef.form <- get_nwparam(dat)$coef.form
  coef.form[1] <- coef.form[1] + log(old.num) - log(new.num)
  dat$nwparam[[1]]$coef.form <- coef.form

  return(dat)
}

