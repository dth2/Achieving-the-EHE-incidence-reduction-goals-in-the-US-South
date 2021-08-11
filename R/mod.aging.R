
#' @title Aging Module
#'
#' @description Module for aging over time for active nodes in the population.
#'
#' @param dat Master data list object of class \code{dat} containing networks,
#'        individual-level attributes, and summary statistics.
#' @param at Current time step.
#'
#' @return
#' This function returns \code{dat} after updating the nodal attribute
#' \code{age} and \code{sqrt.age}. The \code{sqrt.age} vertex attribute is also
#' updated on the three networks.
#'
#' @keywords module msm
#' @export
#'
aging_msm <- function(dat, at) {

  age <- dat$attr$age
  age.grp <- dat$attr$age.grp

  age <- age + 7 / 365

  age.grp[age < 19] <- 1
  age.grp[age >= 19 & age < 25] <- 2
  age.grp[age >= 25 & age < 35] <- 3
  age.grp[age >= 35 & age < 45] <- 4
  age.grp[age >= 45 & age < 55] <- 5
  age.grp[age >= 55] <- 6

  sqrt.age <- sqrt(age)

  dat$attr$age.grp <- age.grp
  dat$attr$age <- age
  dat$attr$sqrt.age <- sqrt.age
  dat$attr$age15 <- ifelse(dat$attr$age < 16, 1, 0)

  return(dat)
}


#' @export
#' @rdname aging_msm
aging_het <- function(dat, at) {

  ## Parameters
  time.unit <- dat$param$time.unit

  ## Attributes
  age <- dat$attr$age
  active <- dat$attr$active

  ## Updates
  age[active == 1] <- age[active == 1] + time.unit/365

  ## Save out
  dat$attr$age <- age

  return(dat)
}
