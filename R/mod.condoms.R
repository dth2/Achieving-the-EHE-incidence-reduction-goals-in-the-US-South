
#' @title Condom Use Module
#'
#' @description Module function stochastically simulates potential condom use
#'              for each act on the discordant edgelist.
#'
#' @inheritParams aging_msm
#'
#' @details
#' For each act on the discordant edgelist, condom use is stochastically simulated
#' based on the partnership type and racial combination of the dyad. Other
#' modifiers for the probability of condom use in that pair are diagnosis of
#' disease, and full or partial HIV viral suppression
#' given HIV anti-retroviral therapy.
#'
#' @return
#' Updates the discordant edgelist with a \code{uai} variable indicating whether
#' condoms were used in that act.
#'
#' @keywords module msm
#' @export
#'
condoms_msm <- function(dat, at) {


  # Attributes
  race <- dat$attr$race
  age <- dat$attr$age
  diag.status <- dat$attr$diag.status
  prepStat <- dat$attr$prepStat

  # Condom Use Models
  cond.mc.mod.het <- dat$param$epistats$cond.mc.mod.het
  cond.oo.mod.het <- dat$param$epistats$cond.oo.mod.het

  cond.mc.mod.msm <- dat$param$epistats$cond.mc.mod.msm
  cond.oo.mod.msm <- dat$param$epistats$cond.oo.mod.msm

  cond.scale.het <- dat$param$cond.scale.het
  cond.scale.msm <- dat$param$cond.scale.msm

  # Temp edgelist
  el <- dat$temp$el

  race.combo <- rep(NA, nrow(el))
  race.combo[race[el[, 1]] == 1 & race[el[, 2]] == 1] <- 1
  race.combo[race[el[, 1]] == 1 & race[el[, 2]] %in% 2:3] <- 2
  race.combo[race[el[, 1]] == 2 & race[el[, 2]] %in% c(1, 3)] <- 3
  race.combo[race[el[, 1]] == 2 & race[el[, 2]] == 2] <- 4
  race.combo[race[el[, 1]] == 3 & race[el[, 2]] %in% 1:2] <- 5
  race.combo[race[el[, 1]] == 3 & race[el[, 2]] == 3] <- 6

  comb.age <- age[el[, 1]] + age[el[, 2]]

  hiv.concord.pos <- rep(0, nrow(el))
  cp <- which(diag.status[el[, 1]] == 1 & diag.status[el[, 2]] == 1)
  hiv.concord.pos[cp] <- 1

  any.prep <- as.numeric((prepStat[el[, 1]] + prepStat[el[, 2]]) > 0)

  ## Main/casual partnerships het
  mc.parts <- which(el[, "ptype"] != 3)
  mc.parts <- which(el[, "ptype"] != 3 & el[, "ptype"] != 6)
  el.mc <- el[mc.parts, ]

  x <- data.frame(ptype = el.mc[, "ptype.het"],
                  duration = el.mc[, "durations"],
                  race.combo = race.combo[mc.parts],
                  comb.age = comb.age[mc.parts],
                  adol = el.mc[, "adol"])
  cond.prob.het <- unname(predict(cond.mc.mod.het, newdata = x, type = "response"))
  el.mc <- cbind(el.mc, cond.prob.het)


  ## Main/casual partnerships MSM

  x <- data.frame(ptype = el.mc[, "ptype.msm"],
                  duration = el.mc[, "durations"],
                  race.combo = race.combo[mc.parts],
                  comb.age = comb.age[mc.parts],
                  adol = el.mc[, "adol"])
  cond.prob.msm <- unname(predict(cond.mc.mod.msm, newdata = x, type = "response"))
  el.mc <- cbind(el.mc, cond.prob.msm)



  ## One-off partnerships het
  oo.parts <- which(el[, "ptype"] == 3 | el[, "ptype"] == 6)
  el.oo <- el[oo.parts, ]

  x <- data.frame(race.combo = race.combo[oo.parts],
                  comb.age = comb.age[oo.parts],
                  adol = el.oo[, "adol"])
  cond.prob.het <- unname(predict(cond.oo.mod.het, newdata = x, type = "response"))
  el.oo <- cbind(el.oo, cond.prob.het)

  ## One-off partnerships msm

  x <- data.frame(race.combo = race.combo[oo.parts],
                  comb.age = comb.age[oo.parts],
                  adol = el.oo[, "adol"])
  cond.prob.msm <- unname(predict(cond.oo.mod.msm, newdata = x, type = "response"))
  el.oo <- cbind(el.oo, cond.prob.msm)


  ## Bind el together
  el <- rbind(el.mc, el.oo)

  el.het <- el[el[, "ptype"] < 4, ]
  el.msm <- el[el[, "ptype"] > 3, ]


  # Acts VI
  vi.vec <- el.het[, "vi"]
  pid.het <- rep(1:length(vi.vec), vi.vec)
  p1.het <- rep(el.het[, "p1"], vi.vec)
  p2.het <- rep(el.het[, "p2"], vi.vec)
  ptype.het <- rep(el.het[, "ptype"], vi.vec)
  cond.prob.het <- rep(el.het[, "cond.prob.het"], vi.vec)

  cond.prob.het <- cond.prob.het * cond.scale.het

  # Acts AI
  ai.vec <- el.msm[, "ai"]
  pid.msm <- rep((length(vi.vec)+1):((length(vi.vec)+length(ai.vec))), ai.vec)
  p1.msm <- rep(el.msm[, "p1"], ai.vec)
  p2.msm <- rep(el.msm[, "p2"], ai.vec)
  ptype.msm <- rep(el.msm[, "ptype"], ai.vec)
  cond.prob.msm <- rep(el.msm[, "cond.prob.msm"], ai.vec)

  cond.prob.msm <- cond.prob.msm * cond.scale.msm

  # UI draw per act
  uvi <- rbinom(length(cond.prob.het), 1, 1 - cond.prob.het)
  uai <- rbinom(length(cond.prob.msm), 1, 1 - cond.prob.msm)

  # Act list construction
  p1 <- c(p1.het,p1.msm)
  p2 <- c(p2.het,p2.msm)
  ptype <- c(ptype.het, ptype.msm)
  ui <- c(uvi, uai)
  pid <- c(pid.het, pid.msm)

  al <- cbind(p1, p2, ptype, ui, pid)
  dat$temp$al <- al

  return(dat)
}
