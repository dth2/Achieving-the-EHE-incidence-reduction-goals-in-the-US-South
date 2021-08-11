
#' @title Sexual Acts Module
#'
#' @description Module function for setting the number of sexual acts on the
#'              discordant edgelist.
#'
#' @inheritParams aging_msm
#'
#' @details
#' The number of acts at each time step is specified as a function of the race of
#' both members in a pair and the expected values within black-black, black-white,
#' and white-white combinations. For one-off partnerships, this is deterministically
#' set at 1, whereas for main and causal partnerships it is a stochastic draw
#' from a Poisson distribution. The number of total acts may further be modified
#' by the level of HIV viral suppression in an infected person.
#'
#' @keywords module msm
#' @export
#'
acts_msm <- function(dat, at) {


  # Attributes
  status <- dat$attr$status
  diag.status <- dat$attr$diag.status
  race <- dat$attr$race
  age <- dat$attr$age
  stage <- dat$attr$stage
  vl <- dat$attr$vl
  uid <- dat$attr$uid
  msm <- dat$attr$msm
  het <- dat$attr$het
  sex <- dat$attr$sex
  age <- dat$attr$age

  plist <- dat$temp$plist

  # Parameters
  acts.mod.het <- dat$param$epistats$acts.mod.het
  acts.mod.msm <- dat$param$epistats$acts.mod.msm
  acts.aids.vl <- dat$param$acts.aids.vl
  acts.scale.het <- dat$param$acts.scale.het
  acts.scale.msm <- dat$param$acts.scale.msm

  # Construct edgelist
  el <- rbind(dat$el[[1]], dat$el[[2]], dat$el[[3]], dat$el[[4]], dat$el[[5]], dat$el[[6]])


  ptype  <- rep(1:6, times = c(nrow(dat$el[[1]]),
                                   nrow(dat$el[[2]]),
                                   nrow(dat$el[[3]]),
                                   nrow(dat$el[[4]]),
                                   nrow(dat$el[[5]]),
                                   nrow(dat$el[[6]])))

  ptype.het  <- rep(c(1:3, 1:3), times = c(nrow(dat$el[[1]]),
                               nrow(dat$el[[2]]),
                               nrow(dat$el[[3]]),
                               nrow(dat$el[[4]]),
                               nrow(dat$el[[5]]),
                               nrow(dat$el[[6]])))

  ptype.msm  <- rep(c(4:6, 4:6), times = c(nrow(dat$el[[1]]),
                                   nrow(dat$el[[2]]),
                                   nrow(dat$el[[3]]),
                                   nrow(dat$el[[4]]),
                                   nrow(dat$el[[5]]),
                                   nrow(dat$el[[6]])))

  st1 <- status[el[, 1]]
  st2 <- status[el[, 2]]
  sex1 <- sex[el[, 1]]
  sex2 <- sex[el[, 2]]
  age1 <- age[el[, 1]]
  age2 <- age[el[, 2]]
  adol <- rep(0,length(age1))
  adol <- ifelse(age1 < 19 | age2 < 19, 1, adol)

  el <- cbind(el, st1, st2, ptype, ptype.het, ptype.msm, sex1, sex2, adol)
  colnames(el) <- c("p1", "p2", "st1", "st2", "ptype", "ptype.het", "ptype.msm", "sex1", "sex2", "adol")
  el <- apply(el, 2 ,as.numeric)

  # Subset to main/casual HET/MSM
  #el.mc.het <- el[el[, "ptype"] != 3, ]
  el.mc <- el[el[, "ptype"] != 3 & el[, "ptype"] != 6, ]

  # Base Act rates based on Poisson model for main/casual
  race.combo <- rep(NA, nrow(el.mc))
  race.combo[race[el.mc[, 1]] == 1 & race[el.mc[, 2]] == 1] <- 1
  race.combo[race[el.mc[, 1]] == 1 & race[el.mc[, 2]] %in% 2:3] <- 2
  race.combo[race[el.mc[, 1]] == 2 & race[el.mc[, 2]] %in% c(1, 3)] <- 3
  race.combo[race[el.mc[, 1]] == 2 & race[el.mc[, 2]] == 2] <- 4
  race.combo[race[el.mc[, 1]] == 3 & race[el.mc[, 2]] %in% 1:2] <- 5
  race.combo[race[el.mc[, 1]] == 3 & race[el.mc[, 2]] == 3] <- 6


  comb.age <- age[el.mc[, 1]] + age[el.mc[, 2]]


  # Current partnership duration
  pid_plist <- plist[, 1]*1e7 + plist[, 2]
  pid_el <- uid[el.mc[, 1]]*1e7 + uid[el.mc[, 2]]
  matches <- match(pid_el, pid_plist)
  durations <- (at - plist[, "start"])[matches]


  # HIV-positive concordant
  hiv.concord.pos <- rep(0, nrow(el.mc))
  cp <- which(diag.status[el.mc[, 1]] == 1 & diag.status[el.mc[, 2]] == 1)
  hiv.concord.pos[cp] <- 1

  # Model predictions VI
  x <- data.frame(ptype = el.mc[, "ptype.het"],
                  duration = durations,
                  race.combo = race.combo,
                  comb.age = comb.age,
                  adol = el.mc[, "adol"])

  rates <- unname(predict(acts.mod.het, newdata = x, type = "response"))/52
  rates <- rates * acts.scale.het
  vi <- rpois(length(rates), rates)
  el.mc <- cbind(el.mc, durations, vi)

  # Model predictions AI
  x <- data.frame(ptype = el.mc[, "ptype.msm"],
                  duration = durations,
                  race.combo = race.combo,
                  comb.age = comb.age,
                  adol = el.mc[, "adol"])

  rates <- unname(predict(acts.mod.msm, newdata = x, type = "response"))/52
  rates <- rates * acts.scale.msm
  ai <- rpois(length(rates), rates)
  el.mc <- cbind(el.mc, ai)

  ##Zero out AI for het and VI for MSM
  el.mc[,"ai"] <- ifelse(el.mc[, "ptype"] == 1,0,el.mc[,"ai"])
  el.mc[,"ai"] <- ifelse(el.mc[, "ptype"] == 2,0,el.mc[,"ai"])

  el.mc[,"vi"] <- ifelse(el.mc[, "ptype"] == 4,0,el.mc[,"vi"])
  el.mc[,"vi"] <- ifelse(el.mc[, "ptype"] == 5,0,el.mc[,"vi"])

  # Add one-time het partnerships
  el.oo.het <- el[el[, "ptype"] == 3, ]
  vi <- durations <- rep(1, nrow(el.oo.het))
  ai <- durations <- rep(0, nrow(el.oo.het))
  el.oo.het <- cbind(el.oo.het, durations, vi, ai)

  # Add one-time msm partnerships
  el.oo.msm <- el[el[, "ptype"] == 6, ]
  vi <- durations <- rep(0, nrow(el.oo.msm))
  ai <- durations <- rep(1, nrow(el.oo.msm))
  el.oo.msm <- cbind(el.oo.msm, durations, vi, ai)


  # Bind el back together
  el <- rbind(el.mc, el.oo.het, el.oo.msm)

  # For AIDS cases with VL above acts.aids.vl, reduce their their acts to 0
  p1HIV <- which(el[, "st1"] == 1)
  p1AIDS <- stage[el[p1HIV, "p1"]] == 4 & vl[el[p1HIV, "p1"]] >= acts.aids.vl
  el[p1HIV[p1AIDS == TRUE], "vi"] <- 0
  el[p1HIV[p1AIDS == TRUE], "ai"] <- 0

  p2HIV <- which(el[, "st2"] == 1)
  p2AIDS <- stage[el[p2HIV, "p2"]] == 4 & vl[el[p2HIV, "p2"]] >= acts.aids.vl
  el[p2HIV[p2AIDS == TRUE], "vi"] <- 0
  el[p2HIV[p2AIDS == TRUE], "ai"] <- 0

  # Flip order of discordant edges
  disc <- abs(el[, "st1"] - el[, "st2"]) == 1
  disc.st2pos <- which(disc == TRUE & el[, "st2"] == 1)
  el[disc.st2pos, 1:4] <- el[disc.st2pos, c(2, 1, 4, 3)]

  # Remove inactive edges from el
  el <- el[-which(el[, "ai"] == 0 & el[, "vi"] == 0), ]


  # Save out
  dat$temp$el <- el

  return(dat)
}
