# TO DO:
## double-sided (done 2/18)
## loss of follow-up
## negative concordance

library(survival)
#library(PWEALL)
source("https://raw.githubusercontent.com/chengs94/SurvPET/master/survROC.R")

surv_enrichment <- function (formula, data, hr = 0.8, end.of.trial=NULL, a=NULL, f=NULL,
                             method = "KM", lambda = 0.05,
                             cost.screening = NULL, cost.keeping = NULL, cost.unit.keeping = NULL,
                             power = 0.9, alpha = 0.05, one.sided = F,
                             selected.biomarker.quantiles = seq(from = 0, to = 0.95, by = 0.05),
                             do.bootstrap = FALSE, n.bootstrap = 1000, seed = 2333,
                             print.summary.tables = FALSE){

  ##### check arguments #####
  if (class(formula) != "formula"){
    stop("invalid formula")
  }
  if (!is.data.frame(data)) {
    stop("dataset must be a data frame")
  }
  #if (length(cost.keeping)!=length(end.of.trial)){
  #  stop("length of trial costs does not match the length of trial duration")
  #}
  if (is.null(end.of.trial) & (is.null(a) | is.null(f))){
    stop("must specify either accrual + follow-up time or length of trial")
  }
  acc.fu <- FALSE
  if (!is.null(a) & !is.null(f))
    acc.fu <- TRUE
  if (method=="NNE") acc.fu <- FALSE
  comp.cost <- TRUE # indicator for whether to compute costs
  if (is.null(cost.screening) | (is.null(cost.keeping) & is.null(cost.unit.keeping)))
    comp.cost <- FALSE
  if (method=="NNE" & (is.null(cost.keeping) | is.null(cost.screening)))
    comp.cost <- FALSE

  formula.vars <- all.vars(formula)
  count.vars <- length(formula.vars)
  data <- data[,which(colnames(data) %in% formula.vars)]
  #ind.missing <- apply(data[, formula.vars], 1, function(x) sum(is.na(x)) > 0)
  ind.missing <- apply(as.matrix(data), 1, function(x) sum(is.na(x)) > 0)
  count.any.missing <- sum(ind.missing)
  if (count.any.missing > 0) {
    data <- data[!ind.missing, ]
    warning(paste(count.any.missing, "observation(s) with missing values were removed"))
  }
  if (!all(formula.vars %in% names(data))) {
    stop(paste("Variable(s) named", formula.vars[which(formula.vars %in%
                                                         names(data) == FALSE)], "were not found in the data"))
  }
  response.name <- formula.vars[1]
  response <- data[, response.name]
  #if (class(response)!="Surv"){
  #  stop("response must be a survival object as returned by the Surv function")
  #}
  features.names <- formula.vars[2:count.vars]
  features <- as.matrix(data[, features.names])
  if (!(power > 0 & power < 1)) {
    stop("power should be between 0 and 1")
  }
  if (!(alpha > 0 & alpha < 1)) {
    stop("alpha should be between 0 and 1")
  }
  if (!(hr > 0 & hr < 1)) {
    stop("hazard ratio should be between 0 and 1")
  }
  if (!all(selected.biomarker.quantiles >= 0 & selected.biomarker.quantiles < 1)) {
    stop("quantiles of the biomarker measured in controls must be at least 0 and less than 1")
  }

  selected.biomarker.quantiles <- c(0, selected.biomarker.quantiles) # no enrichment for reference
  if (one.sided == T) alpha <- alpha*2

  if (method!="KM" & method!="NNE")
    stop("method can only be 'KM' or 'NNE'")
  if (method=="NNE" & (is.null(end.of.trial)))
    stop("trial length missing; NNE can only deal with fixed length trial")
  if (method=="NNE" & !is.null(cost.unit.keeping))
    message("NNE does not deal with unit cost; 'cost.unit.keeping' will be ignored")

  # "pseudo-biomarker" if multiple biomarkers are used
  if(count.vars == 2){
    biomarker.name <- formula.vars[2]
    biomarker <- data[, biomarker.name]
  }
  if (count.vars > 2) {
    biomarker.name <- "combined biomarker"
    coxfit <- do.call(what = coxph, args = list(formula = formula, data = data))
    biomarker <- as.numeric(features%*%coxfit$coefficients)
  }

  ##### separate functions #####
  simps <- function(vec) (vec[1]+4*vec[2]+vec[3])/6 # Simpson's rule
  nevent.calc <- function(hr) 4*(qnorm(1-alpha/2)+qnorm(power))^2/(log(hr))^2 # sample size formula (# events)

  event.rate <- function(enr.level, eot){
    q <- quantile(biomarker,prob=enr.level)
    sobj <- response[biomarker>q]
    km <- survfit(sobj~1,error="greenwood")
    survest <- stepfun(km$time, c(1, km$surv))
    eprob <- 1-survest(eot)
    reftm <- tryCatch(
      {
        max(summary(km)$time[summary(km)$time<=eot]) # find a point where the sd of event rate is equal to that at t0
      },
      error=function(cond) {
        message("length of trial is too short; cannot compute sd of event rate for all enrichment levels (NA returned)")
      },
      warning=function(cond) {
        message("length of trial is too short; cannot compute sd of event rate for all enrichment levels (NA returned)")
      },
      finally={
      }
    )
    esd <- NULL
    esd <- summary(km)$std.err[which(summary(km)$time==reftm)]
    #rmst.obj <- rmsth(y=sobj[,1],d=sobj[,2],tcut=eot,eps=1.0e-06)
    #tmean <- rmst.obj$rmst
    tmean <- survival:::survmean(km, rmean=eot) [[1]]["*rmean"]
    #tmean.sd <- sqrt(rmst.obj$var)
    return(list(eprob = eprob, esd = esd, tmean = tmean))
  }

  event.rate.acc.fu <- function(enr.level,biomarker,response){
    q <- quantile(biomarker,prob=enr.level)
    sobj <- response[biomarker>q]
    km <- survfit(sobj~1,error="greenwood")
    survest <- stepfun(km$time, c(1, km$surv))
    p <- c(survest(a),survest(a+f/2),survest(a+f))
    d.ctrl <- 1-simps(p)
    d.trt <- 1-simps(p^hr)
    d <- (d.ctrl+d.trt)/2
    eprob <- d.ctrl
    tmean <- rep(NA, 3)
    tmean[1] <- survival:::survmean(km, rmean=f) [[1]]["*rmean"]
    tmean[2] <- survival:::survmean(km, rmean=a/2+f) [[1]]["*rmean"]
    tmean[3] <- survival:::survmean(km, rmean=a+f) [[1]]["*rmean"]
    avgtm <- simps(tmean)
    return(list(eprob = eprob, tmean = avgtm, d=d))
  }

  boot.acc.fu <- function(idx){
    rslt.boot <- sapply(selected.biomarker.quantiles, event.rate.acc.fu,
                        biomarker=biomarker.orig[idx], response=response.orig[idx,])
    eprob.boot <- as.numeric(rslt.boot[1,])
    tmean.boot <- as.numeric(rslt.boot[2,])
    d.boot <- as.numeric(rslt.boot[3,])
    return(list(eprob.boot = eprob.boot, tmean.boot = tmean.boot, d.boot=d.boot))
  }

  event.rate.NNE <- function(enr.level, eot){
    q <- quantile(biomarker,prob=enr.level)
    nne <- survROC(Stime=response[,1], status=response[,2], marker=biomarker,
                   predict.time=eot, cut.values=q,
                   method = method, lambda = lambda)
    return(eprob = 1-nne)
  }

  ##### the main function #####
  # 1. K-M left for "plot_surv_enrichment_summaries"
  # 2. event rate and sd
  if (!acc.fu){
    arg <- expand.grid(selected.biomarker.quantiles, end.of.trial)
    if (method!="NNE"){
      rslt <- mapply(event.rate, arg[,1], arg[,2])
      eprob <- matrix(as.numeric(rslt[1,]), ncol=length(end.of.trial))
      esd <- matrix(as.numeric(rslt[2,]), ncol=length(end.of.trial))
      tmean <- matrix(as.numeric(rslt[3,]), ncol=length(end.of.trial))
      #tmeansd <- matrix(as.numeric(rslt[4,]), ncol=length(end.of.trial))
    }
    if (method=="NNE"){
      rslt <- mapply(event.rate.NNE, arg[,1], arg[,2])
      eprob <- matrix(as.numeric(rslt), ncol=length(end.of.trial))
      esd <- NULL
      tmean <- NULL
    }
  }
  if (acc.fu){
    rslt <- sapply(selected.biomarker.quantiles, event.rate.acc.fu,
                   biomarker = biomarker, response = response)
    eprob <- as.numeric(rslt[1,])
    tmean <- as.numeric(rslt[2,])
    esd <- NULL
    d <- as.numeric(rslt[3,])
    biomarker.orig <- biomarker
    response.orig <- response
    if (do.bootstrap){
      set.seed(seed)
      idx.boot <- replicate(n.bootstrap,
                            sample(seq(1,length(biomarker)), length(biomarker),replace=T))
      idx.boot <- split(idx.boot, rep(1:n.bootstrap, each = length(biomarker)))
      stat.boot <- lapply(idx.boot, boot.acc.fu)
      eprob.boot <- do.call(rbind, lapply(stat.boot, '[[', 1))
      tmean.boot <- do.call(rbind, lapply(stat.boot, '[[', 2))
      d.boot <- do.call(rbind, lapply(stat.boot, '[[', 3))
      esd <- apply(eprob.boot, 2, sd)
    }
  }
  # 3. total sample size
  nevent <- nevent.calc(hr)
  if (!acc.fu){
    prob.all <- 1-(1-eprob)^hr+eprob
    npat <- ceiling(nevent*2/prob.all)
    sd.npat <- NULL
    if (method!="NNE"){
      sd.all <- (1+hr*(1-eprob)^(hr-1))*esd
      sd.npat <- npat/prob.all*sd.all
    }
  }
  if (acc.fu){
    nevent <- nevent.calc(hr)
    npat<-ceiling(nevent/d)
    sd.npat <- NULL
    if (do.bootstrap){
      npat.boot <- ceiling(nevent/d.boot)
      sd.npat <- apply(npat.boot, 2, sd)
    }
  }
  # 4. total patients screened = trial ss/(1-p)
  sd.nscr <-NULL
  nscr <- npat/(1-selected.biomarker.quantiles)
  if (!is.null(sd.npat)){
    sd.nscr <- sd.npat/(1-selected.biomarker.quantiles)
  }
  nscr.orig <- nscr
  sd.nscr.orig <- sd.nscr

  ## 5 & 6: costs
  if (!comp.cost){
    cost <- sd.cost <- reduc <- sd.reduc <- NULL
  }
  if (comp.cost){
    # 5. total costs for screening and trial
    nscr[selected.biomarker.quantiles == 0] <- 0
    if (!is.null(sd.npat)){
      sd.nscr[selected.biomarker.quantiles == 0] <- 0
    }
    if (!acc.fu){
      cost <- sd.cost <- matrix(NA, ncol = length(end.of.trial), nrow = length(selected.biomarker.quantiles))
      if (is.null(cost.keeping)){#cost.keeping <- end.of.trial*cost.unit.keeping
        # follow-up to expected occurrence of event
        sd.cost <- NULL
        for (k in 1:length(end.of.trial)){
          cost[,k] <- cost.unit.keeping[k]*npat[,k]*tmean[,k]+cost.screening*nscr[,k]
          #sd.cost[,k] <- (cost.keeping[k]+cost.screening/(1-selected.biomarker.quantiles))*sd.npat[,k]
        }
        idx <- which(selected.biomarker.quantiles==0) # no screening needed
        if (length(idx)!=0){
          for (k in 1:length(end.of.trial)){
            cost[idx,k] <- cost.unit.keeping[k]*npat[idx,k]*tmean[idx,k]
            #sd.cost[idx,k] <- sd.npat[idx,k]*cost.keeping[k]
          }
        }
      }
      if (!is.null(cost.keeping)){
        if (is.null(sd.npat)) sd.cost <- NULL
        # follow-up to the end of trial
        for (k in 1:length(end.of.trial)){
          cost[,k] <- cost.keeping[k]*npat[,k]+cost.screening*nscr[,k]
          if (!is.null(sd.npat))
            sd.cost[,k] <- (cost.keeping[k]+cost.screening/(1-selected.biomarker.quantiles))*sd.npat[,k]
        }
        idx <- which(selected.biomarker.quantiles==0) # no screening needed
        if (length(idx)!=0){
          for (k in 1:length(end.of.trial)){
            cost[idx,k] <- cost.keeping[k]*npat[idx,k]
            if (!is.null(sd.npat))
              sd.cost[idx,k] <- sd.npat[idx,k]*cost.keeping[k]
          }
        }
      }
    }
    if (acc.fu){
      sd.cost <- NULL
      if (!is.null(cost.unit.keeping)){
        cost <- cost.screening*nscr + cost.unit.keeping*tmean*npat
        if (do.bootstrap){
          nscr.boot <- t(apply(npat.boot, 1, function(vec) vec/(1-selected.biomarker.quantiles)))
          cscr.boot <- cost.screening*nscr.boot
          cscr.boot[,selected.biomarker.quantiles == 0] <- 0
          ckeep.boot <- cost.unit.keeping*tmean.boot*npat.boot
          cost.boot <- cscr.boot+ckeep.boot
          sd.cost <- apply(cost.boot, 2, sd)
        }
      }
      if (is.null(cost.unit.keeping)){
        cost <- cost.screening*nscr + cost.keeping*npat
        if (do.bootstrap){
          nscr.boot <- t(apply(npat.boot, 1, function(vec) vec/(1-selected.biomarker.quantiles)))
          cscr.boot <- cost.screening*nscr.boot
          cscr.boot[,selected.biomarker.quantiles == 0] <- 0
          ckeep.boot <- cost.keeping*npat.boot
          cost.boot <- cscr.boot+ckeep.boot
          sd.cost <- apply(cost.boot, 2, sd)
        }
      }
    }
    # 6. % reduction in total cost
    sd.reduc <- NULL
    if (!acc.fu){
      reduc <- matrix(NA, ncol = length(end.of.trial), nrow = length(selected.biomarker.quantiles))
      for (k in 1:length(end.of.trial)){
        if (!is.null(cost.keeping))
          reduc[,k] <- (cost.keeping[k]*npat[1,k]-cost[,k])/cost.keeping[k]/npat[1,k]
        if (is.null(cost.keeping))
          reduc[,k] <- (cost.unit.keeping[k]*npat[1,k]*tmean[1,k]-cost[,k])/cost.unit.keeping[k]/npat[1,k]/tmean[1,k]
      }
    }
    if (acc.fu){
      reduc <- (cost[1]-cost)/cost[1]
      if (do.bootstrap){
        reduc.boot <- apply(cost.boot, 2, function(vec) (cost.boot[,1]-vec)/cost.boot[,1])
        sd.reduc <- apply(reduc.boot, 2,sd)
      }
    }
  }

  # remove reference values
  ind.matrix <- class(eprob)=="matrix"
  if (ind.matrix){
    eprob <- eprob[-1,]
    esd <- esd[-1,]
    npat <- npat[-1,]
    sd.npat <- sd.npat[-1,]
    nscr.orig <- nscr.orig[-1,]
    sd.nscr.orig <- sd.nscr.orig[-1,]
    cost <- cost[-1,]
    sd.cost <- sd.cost[-1,]
    reduc <- reduc[-1,]
  }
  if (!ind.matrix){
    eprob <- eprob[-1]
    esd <- esd[-1]
    npat <- npat[-1]
    sd.npat <- sd.npat[-1]
    nscr.orig <- nscr.orig[-1]
    sd.nscr.orig <- sd.nscr.orig[-1]
    cost <- cost[-1]
    sd.cost <- sd.cost[-1]
    reduc <- reduc[-1]
    sd.reduc <- sd.reduc[-1]
  }
  selected.biomarker.quantiles <- selected.biomarker.quantiles[-1]

  # print tables
  cnames <- c("level.enrichment","prob","prob.sd",
              "n.patients","n.patients.sd",
              "n.screened","n.screened.sd",
              "cost","cost.sd",
              "reduction","reduction.sd")
  if (acc.fu){
    tab <- cbind(selected.biomarker.quantiles, eprob, esd, npat, sd.npat, nscr.orig, sd.nscr.orig,
                 cost, sd.cost, reduc, sd.reduc)
    if (comp.cost & do.bootstrap)
      colnames(tab) <- cnames
    if (comp.cost & !do.bootstrap)
      colnames(tab) <- cnames[c(1,2,4,6,8,10)]
    if (!comp.cost & do.bootstrap)
      colnames(tab) <- cnames[1:7]
    if (!comp.cost & !do.bootstrap)
      colnames(tab) <- cnames[c(1,2,4,6)]
  }
  if (!acc.fu){
    tab <- rep(list(NULL),length(end.of.trial))
    if (length(end.of.trial) > 1){
      for (j in 1:length(end.of.trial)){
        tab[[j]] <- cbind(selected.biomarker.quantiles, eprob[,j], esd[,j],
                          npat[,j], sd.npat[,j],
                          nscr.orig[,j], sd.nscr.orig[,j],
                          cost[,j], sd.cost[,j], reduc[,j])
        if (comp.cost & !is.null(cost.keeping) & !is.null(esd)) colnames(tab[[j]]) <- cnames[1:10]
        if (comp.cost & is.null(cost.keeping) & !is.null(esd)) colnames(tab[[j]]) <- cnames[c(1:8,10)]
        if (!comp.cost & !is.null(esd)) colnames(tab[[j]]) <- cnames[1:7]
        if (comp.cost & is.null(esd)) colnames(tab[[j]]) <- cnames[c(1,2,4,6,8,10)]
        if (!comp.cost & is.null(esd)) colnames(tab[[j]]) <- cnames[c(1,2,4,6)]
      }
      names(tab) <- paste("trial length=", end.of.trial)
    }
    if (length(end.of.trial) == 1){
      tab <- cbind(selected.biomarker.quantiles, eprob, esd, npat, sd.npat, nscr.orig, sd.nscr.orig,
                   cost, sd.cost, reduc)
      if (comp.cost & !is.null(esd)) colnames(tab) <- cnames[1:10]
      if (!comp.cost & !is.null(esd)) colnames(tab) <- cnames[1:7]
      if (comp.cost & is.null(esd)) colnames(tab) <- cnames[c(1,2,4,6,8,10)]
      if (!comp.cost & is.null(esd)) colnames(tab) <- cnames[c(1,2,4,6)]
    }
  }
  if(print.summary.tables) print(tab)

  return(list(summary.table=tab,
              event.prob=eprob,event.prob.sd=esd,
              n.patients=npat, n.patients.sd=sd.npat,
              num.screened=nscr.orig, num.screened.sd=sd.nscr.orig,
              total.cost=cost, total.cost.sd=sd.cost,
              cost.reduction=reduc, cost.reduction.sd=sd.reduc,
              response=response, biomarker=biomarker, biomarker.name=biomarker.name,
              selected.biomarker.quantiles=selected.biomarker.quantiles,
              end.of.trial=end.of.trial,
              a=a, f=f, acc.fu=acc.fu))
}
