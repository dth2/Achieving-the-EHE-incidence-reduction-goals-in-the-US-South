
#' @title Curing Module
#'
#' @description Module for curing at set timepoints to hit prevalence.
#'
#' @param dat Master data list object of class \code{dat} containing networks,
#'        individual-level attributes, and summary statistics.
#' @param at Current time step.
#'
#' @return
#' This function returns \code{dat} after updating the nodal attribute
#' \code{diag.status} and \code{status}, \code{inf.time}, \code{stage} \code{stage.time} \code{aids.time} \code{diag.stage}
#' \code{vl} \code{vl.last.usupp} \code{vl.last.supp} \code{diag.time}  \code{last.neg.test} \code{tx.status} \code{cuml.time.on.tx}
#' \code{cuml.time.off.tx}  \code{tx.period.first}   \code{tx.period.last}    \code{tx.init.time}      \code{count.trans} .
#'
#' @keywords module
#' @export
#'

cure_prev <- function(dat, at) {


#parameters
cure.time <- dat$param$cure.time
cure.time.w <- dat$param$cure.time.w
cure.time.msm <- dat$param$cure.time.msm
prev.targ <- dat$param$prev.targ


if(at %in% cure.time){

  #Attributes
  dem.cat <- dat$attr$dem.cat
  status <- dat$attr$status

el <- dat$el

cure.list <- NULL
cure.ids <- NULL

for(i in 1:9){

  #Select indexes to cure
  #select half the difference between curent case counts and targets to because partners will be cured increasing counts
  #and some partner will cross catagories and we don't want to over cure

  group <- which(dem.cat==i & status ==1)
  pos.count <- length(group)

  if(i > 3){
  cure.count <- round((max(0,pos.count - prev.targ[i]))/4)
  cure.ids <- sample(group,cure.count,replace=FALSE)
  }

  if(i <= 3){
    cure.count <- round((max(0,pos.count - prev.targ[i]))/3)
    cure.ids <- sample(group,cure.count,replace=FALSE)
  }


##el is index not uid

  #main partners
  #Heterosexual

  if(i > 3){
    p1.id.list <- c(dat$el[[1]][,1])
    p2.id.list <- c(dat$el[[1]][,2])
    p1 <- intersect( p1.id.list,cure.ids)
    p2 <- intersect( p2.id.list,cure.ids)

    alters1 <- which(dat$el[[1]][,1] %in% p1)
    alters1 <- dat$el[[1]][,2][alters1]
    alters1 <- as.vector(alters1)

    alters2 <- which(dat$el[[1]][,2] %in% p2)
    alters2 <- dat$el[[1]][,1][alters2]
    alters2 <-as.vector(alters2)

    cure.list <- c(cure.list, cure.ids, alters1, alters2)
    cure.list <- unique(cure.list)
  }

  #MSM
  if(i < 4){
  p1.id.list <- c(dat$el[[4]][,1])
  p2.id.list <- c(dat$el[[4]][,2])
  p1 <- intersect( p1.id.list,cure.ids)
  p2 <- intersect( p2.id.list,cure.ids)

  alters1 <- which(dat$el[[4]][,1] %in% p1)
  alters1 <- dat$el[[4]][,2][alters1]
  alters1 <- as.vector(alters1)

  alters2 <- which(dat$el[[4]][,2] %in% p2)
  alters2 <- dat$el[[4]][,1][alters2]
  alters2 <-as.vector(alters2)

  cure.list <- c(cure.list, cure.ids, alters1, alters2)
  cure.list <- unique(cure.list)
  }

  #casual partners
  #Heterosexual

  if(i > 3){
    p1.id.list <- c(dat$el[[2]][,1])
    p2.id.list <- c(dat$el[[2]][,2])
    p1 <- intersect( p1.id.list,cure.ids)
    p2 <- intersect( p2.id.list,cure.ids)

    alters1 <- which(dat$el[[2]][,1] %in% p1)
    alters1 <- dat$el[[2]][,2][alters1]
    alters1 <- as.vector(alters1)

    alters2 <- which(dat$el[[2]][,2] %in% p2)
    alters2 <- dat$el[[2]][,1][alters2]
    alters2 <- as.vector(alters2)

    cure.list <- c(cure.list, cure.ids, alters1, alters2)
    cure.list <- unique(cure.list)
  }

  #MSM
  if(i < 4){
    p1.id.list <- c(dat$el[[5]][,1])
    p2.id.list <- c(dat$el[[5]][,2])
    p1 <- intersect( p1.id.list,cure.ids)
    p2 <- intersect( p2.id.list,cure.ids)

    alters1 <- which(dat$el[[5]][,1] %in% p1)
    alters1 <- dat$el[[5]][,2][alters1]
    alters1 <- as.vector(alters1)

    alters2 <- which(dat$el[[5]][,2] %in% p2)
    alters2 <- dat$el[[5]][,1][alters2]
    alters2 <- as.vector(alters2)

    cure.list <- c(cure.list, cure.ids, alters1, alters2)
    cure.list <- unique(cure.list)
  }

}




dat$attr$status[cure.list] <- 0
dat$attr$inf.time[cure.list] <- NA
dat$attr$vl[cure.list] <- NA
dat$attr$stage[cure.list] <- NA
dat$attr$stage.time[cure.list] <- NA
dat$attr$diag.status[cure.list] <- 0
dat$attr$tx.status[cure.list] <- 0
dat$attr$cuml.time.on.tx[cure.list] <- 0
dat$attr$cuml.time.off.tx[cure.list] <- 0

dat$attr$tx.period.first[cure.list] <- NA
dat$attr$tx.period.last[cure.list] <- NA
dat$attr$tx.init.time[cure.list] <- NA

dat$attr$vl.last.usupp[cure.list] <- NA
dat$attr$vl.last.supp[cure.list] <- NA




}


if(at %in% cure.time.w){


    #Attributes
  dem.cat <- dat$attr$dem.cat
  status <- dat$attr$status

  el <- dat$el
  cure.list <- NULL
  cure.ids <- NULL
  cure.ids.het <- NULL
  cure.ids.msm <- NULL


  #Select indexes to cure
  #select half the difference between curent case counts and targets to because partners will be cured increasing counts
  #and some partner will cross catagories and we don't want to over cure

  group.msm <- which(dem.cat==3 & status ==1)
  pos.count.msm <- length(group.msm)

  group.het <- which((dem.cat==6 |dem.cat==9)  & status ==1)
  pos.count.het <- length(group.het)


  cure.count.msm <- round((max(0,pos.count.msm - prev.targ[3]))/3)
  cure.ids.msm <- sample(group.msm,cure.count.msm,replace=FALSE)


  cure.count.het <- round((max(0,pos.count.het - (prev.targ[6] + prev.targ[9])))/3)
  cure.ids.het <- sample(group.het,cure.count.het,replace=FALSE)



  ##el is index not uid

  #main partners
  #Heterosexual

  p1.id.list <- c(dat$el[[1]][,1])
  p2.id.list <- c(dat$el[[1]][,2])
  p1 <- intersect( p1.id.list,cure.ids.het)
  p2 <- intersect( p2.id.list,cure.ids.het)

  alters1 <- which(dat$el[[1]][,1] %in% p1)
  alters1 <- dat$el[[1]][,2][alters1]
  alters1 <- as.vector(alters1)

  alters2 <- which(dat$el[[1]][,2] %in% p2)
  alters2 <- dat$el[[1]][,1][alters2]
  alters2 <-as.vector(alters2)

  cure.list <- c(cure.list, cure.ids.het, alters1, alters2)
  cure.list <- unique(cure.list)


  #MSM
  p1.id.list <- c(dat$el[[4]][,1])
  p2.id.list <- c(dat$el[[4]][,2])
  p1 <- intersect( p1.id.list,cure.ids.msm)
  p2 <- intersect( p2.id.list,cure.ids.msm)

  alters1 <- which(dat$el[[4]][,1] %in% p1)
  alters1 <- dat$el[[4]][,2][alters1]
  alters1 <- as.vector(alters1)

  alters2 <- which(dat$el[[4]][,2] %in% p2)
  alters2 <- dat$el[[4]][,1][alters2]
  alters2 <-as.vector(alters2)

  cure.list <- c(cure.list, cure.ids.msm, alters1, alters2)
  cure.list <- unique(cure.list)


  #casual partners
  #Heterosexual

  p1.id.list <- c(dat$el[[2]][,1])
  p2.id.list <- c(dat$el[[2]][,2])
  p1 <- intersect( p1.id.list,cure.ids.het)
  p2 <- intersect( p2.id.list,cure.ids.het)

  alters1 <- which(dat$el[[2]][,1] %in% p1)
  alters1 <- dat$el[[2]][,2][alters1]
  alters1 <- as.vector(alters1)

  alters2 <- which(dat$el[[2]][,2] %in% p2)
  alters2 <- dat$el[[2]][,1][alters2]
  alters2 <- as.vector(alters2)

  cure.list <- c(cure.list, cure.ids.het, alters1, alters2)
  cure.list <- unique(cure.list)


  #MSM

  p1.id.list <- c(dat$el[[5]][,1])
  p2.id.list <- c(dat$el[[5]][,2])
  p1 <- intersect( p1.id.list,cure.ids.msm)
  p2 <- intersect( p2.id.list,cure.ids.msm)

  alters1 <- which(dat$el[[5]][,1] %in% p1)
  alters1 <- dat$el[[5]][,2][alters1]
  alters1 <- as.vector(alters1)

  alters2 <- which(dat$el[[5]][,2] %in% p2)
  alters2 <- dat$el[[5]][,1][alters2]
  alters2 <- as.vector(alters2)

  cure.list <- c(cure.list, cure.ids.msm, alters1, alters2)
  cure.list <- unique(cure.list)



  dat$attr$status[cure.list] <- 0
  dat$attr$inf.time[cure.list] <- NA
  dat$attr$vl[cure.list] <- NA
  dat$attr$stage[cure.list] <- NA
  dat$attr$stage.time[cure.list] <- NA
  dat$attr$diag.status[cure.list] <- 0
  dat$attr$tx.status[cure.list] <- 0
  dat$attr$cuml.time.on.tx[cure.list] <- 0
  dat$attr$cuml.time.off.tx[cure.list] <- 0

  dat$attr$tx.period.first[cure.list] <- NA
  dat$attr$tx.period.last[cure.list] <- NA
  dat$attr$tx.init.time[cure.list] <- NA

  dat$attr$vl.last.usupp[cure.list] <- NA
  dat$attr$vl.last.supp[cure.list] <- NA

}


if(at %in% cure.time.msm){

  #Attributes
  dem.cat <- dat$attr$dem.cat
  status <- dat$attr$status

  el <- dat$el

  cure.list <- NULL
  cure.ids <- NULL

  for(i in 1:3){

    #Select indexes to cure
    #select half the difference between curent case counts and targets to because partners will be cured increasing counts
    #and some partner will cross catagories and we don't want to over cure

    group <- which(dem.cat==i & status ==1)
    pos.count <- length(group)


    if(i <= 3){
      cure.count <- round((max(0,pos.count - prev.targ[i]))/2)
      cure.ids <- sample(group,cure.count,replace=FALSE)
    }


    ##el is index not uid


    #MSM
    if(i < 4){
      p1.id.list <- c(dat$el[[4]][,1])
      p2.id.list <- c(dat$el[[4]][,2])
      p1 <- intersect( p1.id.list,cure.ids)
      p2 <- intersect( p2.id.list,cure.ids)

      alters1 <- which(dat$el[[4]][,1] %in% p1)
      alters1 <- dat$el[[4]][,2][alters1]
      alters1 <- as.vector(alters1)

      alters2 <- which(dat$el[[4]][,2] %in% p2)
      alters2 <- dat$el[[4]][,1][alters2]
      alters2 <-as.vector(alters2)

      cure.list <- c(cure.list, cure.ids, alters1, alters2)
      cure.list <- unique(cure.list)
    }

    #casual partners
    #MSM
    if(i < 4){
      p1.id.list <- c(dat$el[[5]][,1])
      p2.id.list <- c(dat$el[[5]][,2])
      p1 <- intersect( p1.id.list,cure.ids)
      p2 <- intersect( p2.id.list,cure.ids)

      alters1 <- which(dat$el[[5]][,1] %in% p1)
      alters1 <- dat$el[[5]][,2][alters1]
      alters1 <- as.vector(alters1)

      alters2 <- which(dat$el[[5]][,2] %in% p2)
      alters2 <- dat$el[[5]][,1][alters2]
      alters2 <- as.vector(alters2)

      cure.list <- c(cure.list, cure.ids, alters1, alters2)
      cure.list <- unique(cure.list)
    }

  }




  dat$attr$status[cure.list] <- 0
  dat$attr$inf.time[cure.list] <- NA
  dat$attr$vl[cure.list] <- NA
  dat$attr$stage[cure.list] <- NA
  dat$attr$stage.time[cure.list] <- NA
  dat$attr$diag.status[cure.list] <- 0
  dat$attr$tx.status[cure.list] <- 0
  dat$attr$cuml.time.on.tx[cure.list] <- 0
  dat$attr$cuml.time.off.tx[cure.list] <- 0

  dat$attr$tx.period.first[cure.list] <- NA
  dat$attr$tx.period.last[cure.list] <- NA
  dat$attr$tx.init.time[cure.list] <- NA

  dat$attr$vl.last.usupp[cure.list] <- NA
  dat$attr$vl.last.supp[cure.list] <- NA




}


  return(dat)
}

