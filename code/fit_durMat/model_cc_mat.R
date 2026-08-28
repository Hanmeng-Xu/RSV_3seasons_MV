
# in this model, maternal vaccine was administered to mothers and infants born to vaccinated mothers (aged 20-39 years old pregnant women) are protected (mother receiving the vaccine are also protected)
# this model separates the immunized and unimmunized cohort to track incidence separately

# ==============================================================================
# model_cc_mat -- maternal RSV vaccine, immunised / unimmunised cohorts tracked
# separately (no SES stratification).
#
# CHANGE FROM THE PREVIOUS VERSION -------------------------------------------
# Mothers entering the tracked cohort are now drawn from S3 and W3 in PROPORTION
# to occupancy, and the 0.5 factor means only what it says: half of the drawn
# mothers are vaccinated (-> Vmat) and half are not (-> SU3 / WU3).
#
#   total mothers drawn   T      = v_mat * total births
#   share from S3                = S3 / (S3 + W3)
#   share from W3                = W3 / (S3 + W3)
#   S3 -> Vmat  = 0.5 * share_S3 * T      S3 -> SU3 = 0.5 * share_S3 * T
#   W3 -> Vmat  = 0.5 * share_W3 * T      W3 -> WU3 = 0.5 * share_W3 * T
#
# Two consequences, both intended:
#   (a) TOTAL mothers moved is now T, not 2*T. The previous version removed
#       0.5*T from each of S3 and W3 for EACH destination, i.e. 2*T in total,
#       while only T newborns entered the cohort (0.5*T to MVmat, 0.5*T to MU).
#       Results will differ from earlier runs.
#   (b) Because the draw is proportional, neither S3 nor W3 can be pushed
#       negative by the transfer, which the absolute version could do.
#
# Mass balance: out = share_S3*T + share_W3*T = T
#               in  = Vmat 0.5*T + SU3 0.5*share_S3*T + WU3 0.5*share_W3*T = T
#
# although i tested, this doesn't really matter as it only influences the infections in adults (pregnant mothers aged 20-39 yo) but not infants under 6m
# ==============================================================================

model_cc_mat <- function(t, y, parms, time.step = "month") {
  
  States <- array(y, dim = dim(parms$yinit.matrix))
  dimnames(States) <- dimnames(parms$yinit.matrix)
  
  if (parms$time.step == 'month') {
    period <- 12;       length.step <- 30.44
  } else if (parms$time.step == 'week') {
    period <- 52.1775;  length.step <- 7
  } else {
    stop("time.step must be 'month' or 'week'")
  }
  
  N.ages <- nrow(States)
  
  # waning rates
  omega     <- 1/(parms$DurationMatImmunityDays/length.step)  # maternal protection
  omega_mat <- 1/(parms$DurationAbrImmunityDays/length.step)  # vaccine: SVmat->SV0, Vmat->SV3
  
  # ---- time index -----------------------------------------------------------
  # lsoda evaluates at continuous t, sometimes slightly beyond the current
  # output point. Clamping prevents RR[t] / VacPro[t] silently returning NA
  # (a vector indexed past its end returns NA without warning).
  ti <- min(max(floor(t), 1),
            length(parms$RR), length(parms$rrM), length(parms$VacPro),
            nrow(as.matrix(parms$PerCapitaBirthsYear)))
  
  RR  <- parms$RR[ti]      # relative risk of infection post-vax, infants
  rrM <- parms$rrM[ti]     # relative risk of infection post-vax, mothers
  
  V     <- parms$VacPro[ti]
  c_mat <- parms$cover_mat
  v_mat <- V * c_mat
  
  # birth, death, aging
  period.birth.rate <- log(parms$PerCapitaBirthsYear[ti, ] + 1)/period
  um <- parms$um
  mu <- 1/parms$WidthAgeClassMonth
  if (parms$time.step == 'week') mu <- 1/(parms$WidthAgeClassMonth * 4.345)
  Aging.Prop <- c(0, mu[1:(N.ages - 1)])
  
  gamma1 <- parms$gamma1
  gamma2 <- parms$gamma2
  gamma3 <- parms$gamma3
  gamma4 <- parms$gamma4
  
  sigma1 <- parms$sigma1
  sigma2 <- parms$sigma2
  sigma3 <- parms$sigma3
  
  theta1 <- if (parms$AllowWaning == 'Yes') {
    1/(parms$dur.immunity1/length.step)
  } else 0
  
  # ---- which age class(es) are the childbearing mothers? --------------------
  # Distribution over ages, summing to 1. Defaults to the 11th class (20-39Y)
  # for a 13-class structure, reproducing the old c(rep(0,10), x, rep(0,2)).
  if (!is.null(parms$mother_ages)) {
    mother_age <- parms$mother_ages
    if (is.character(mother_age)) {
      bad <- setdiff(mother_age, rownames(States))
      if (length(bad) > 0)
        stop("mother_ages not found among age labels: ",
             paste(bad, collapse = ", "))
      w <- as.numeric(rownames(States) %in% mother_age)
      mother_age <- w/sum(w)
    }
    if (length(mother_age) != N.ages) stop("mother_ages must have length ", N.ages)
  } else {
    if (N.ages != 13)
      stop("no parms$mother_ages supplied and N.ages != 13, so the default ",
           "(11th age class = 20-39Y) cannot be assumed")
    mother_age <- c(rep(0, 10), 1, rep(0, 2))
  }
  
  # ---- states ---------------------------------------------------------------
  # status quo
  M  <- States[, 'M']
  S0 <- States[, 'S0']; I1 <- States[, 'I1']
  S1 <- States[, 'S1']; I2 <- States[, 'I2']
  S2 <- States[, 'S2']; I3 <- States[, 'I3']
  S3 <- States[, 'S3']; I4 <- States[, 'I4']
  W1 <- States[, 'W1']; W2 <- States[, 'W2']; W3 <- States[, 'W3']
  
  # unimmunised arm
  MU  <- States[, 'MU']
  SU0 <- States[, 'SU0']; IU1 <- States[, 'IU1']
  SU1 <- States[, 'SU1']; IU2 <- States[, 'IU2']
  SU2 <- States[, 'SU2']; IU3 <- States[, 'IU3']
  SU3 <- States[, 'SU3']; IU4 <- States[, 'IU4']
  WU1 <- States[, 'WU1']; WU2 <- States[, 'WU2']; WU3 <- States[, 'WU3']
  
  # immunised arm
  MVmat <- States[, 'MVmat']   # newborns of vaccinated mothers
  SVmat <- States[, 'SVmat']   # after maternal immunity wanes, still protected
  Vmat  <- States[, 'Vmat']    # vaccinated mothers
  SV0 <- States[, 'SV0']; IV1 <- States[, 'IV1']
  SV1 <- States[, 'SV1']; IV2 <- States[, 'IV2']
  SV2 <- States[, 'SV2']; IV3 <- States[, 'IV3']
  SV3 <- States[, 'SV3']; IV4 <- States[, 'IV4']
  WV1 <- States[, 'WV1']; WV2 <- States[, 'WV2']; WV3 <- States[, 'WV3']
  
  # ---- force of infection ---------------------------------------------------
  seasonal.txn <- 1 + parms$b1 * cos(2*pi*(t - parms$phi*period)/period)
  transmission_unittime <- parms$baseline.txn.rate/(parms$dur.days1/length.step)
  beta     <- transmission_unittime * parms$c2
  beta_a_i <- seasonal.txn * beta/(sum(States)^parms$q)
  
  infectiousN <- I1  + parms$rho1*I2  + parms$rho2*I3  + parms$rho2*I4 +
    IV1 + parms$rho1*IV2 + parms$rho2*IV3 + parms$rho2*IV4 +
    IU1 + parms$rho1*IU2 + parms$rho2*IU3 + parms$rho2*IU4
  
  lambda <- as.vector(infectiousN %*% beta_a_i)
  
  # ==== MOTHERS ENTERING THE TRACKED COHORT =================================
  # Total births this step. period.birth.rate[1] in the old code is equivalent
  # to sum() whenever a single age column of PerCapitaBirthsYear is nonzero,
  # which dy[,'M'] already assumes.
  total_birth <- sum(period.birth.rate) * sum(States)
  
  # one mother per tracked pregnancy
  T_mothers <- mother_age * (v_mat * total_birth)
  
  # proportional draw between S3 and W3 (both zero if the pool is empty)
  den     <- S3 + W3
  share_S <- ifelse(den > 0, S3/den, 0)
  share_W <- ifelse(den > 0, W3/den, 0)
  
  from_S3 <- share_S * T_mothers          # leaves S3 in total
  from_W3 <- share_W * T_mothers          # leaves W3 in total
  
  # half of each stream is vaccinated, half is not
  S3_to_Vmat <- 0.5 * from_S3
  S3_to_SU3  <- 0.5 * from_S3
  W3_to_Vmat <- 0.5 * from_W3
  W3_to_WU3  <- 0.5 * from_W3
  # ===========================================================================
  
  dy <- matrix(NA_real_, nrow = N.ages, ncol = ncol(States))
  colnames(dy) <- colnames(States)
  
  # ------------ status quo ---------------------------------------------------
  dy[, 'M'] <- (1 - v_mat) * period.birth.rate * sum(States) -
    (omega + (mu + um))*M +
    Aging.Prop*c(0, M[1:(N.ages - 1)])
  
  dy[, 'S0'] <- theta1*W1 +
    omega*M -
    lambda*S0 -
    (mu + um)*S0 +
    Aging.Prop*c(0, S0[1:(N.ages - 1)])
  
  dy[, 'I1'] <- lambda*S0 -
    (gamma1 + mu + um)*I1 +
    Aging.Prop*c(0, I1[1:(N.ages - 1)])
  
  dy[, 'S1'] <- gamma1*I1 -
    sigma1*lambda*S1 -
    (mu + um)*S1 +
    Aging.Prop*c(0, S1[1:(N.ages - 1)]) +
    theta1*W2 -
    theta1*S1
  
  dy[, 'W1'] <- theta1*S1 -
    sigma1*lambda*W1 -
    (mu + um)*W1 +
    Aging.Prop*c(0, W1[1:(N.ages - 1)]) -
    theta1*W1
  
  dy[, 'I2'] <- sigma1*lambda*S1 -
    gamma2*I2 -
    (mu + um)*I2 +
    Aging.Prop*c(0, I2[1:(N.ages - 1)]) +
    sigma1*lambda*W1
  
  dy[, 'S2'] <- gamma2*I2 -
    sigma2*lambda*S2 -
    (mu + um)*S2 +
    Aging.Prop*c(0, S2[1:(N.ages - 1)]) +
    theta1*W3 -
    theta1*S2
  
  dy[, 'W2'] <- theta1*S2 -
    sigma2*lambda*W2 -
    (mu + um)*W2 +
    Aging.Prop*c(0, W2[1:(N.ages - 1)]) -
    theta1*W2
  
  dy[, 'I3'] <- sigma2*lambda*S2 -
    (gamma3 + mu + um)*I3 +
    Aging.Prop*c(0, I3[1:(N.ages - 1)]) +
    sigma2*lambda*W2
  
  dy[, 'S3'] <- gamma3*I3 +
    gamma4*I4 -
    sigma3*lambda*S3 -
    (mu + um)*S3 +
    Aging.Prop*c(0, S3[1:(N.ages - 1)]) -
    theta1*S3 -
    S3_to_Vmat -                      # S3 -> Vmat   (vaccinated mothers)
    S3_to_SU3                         # S3 -> SU3    (unvaccinated arm)
  
  dy[, 'W3'] <- theta1*S3 -
    sigma3*lambda*W3 -
    (mu + um)*W3 +
    Aging.Prop*c(0, W3[1:(N.ages - 1)]) -
    theta1*W3 -
    W3_to_Vmat -                      # W3 -> Vmat
    W3_to_WU3                         # W3 -> WU3
  
  dy[, 'I4'] <- sigma3*lambda*S3 +
    sigma3*lambda*W3 -
    gamma4*I4 -
    (mu + um)*I4 +
    Aging.Prop*c(0, I4[1:(N.ages - 1)])
  
  # ------------ immunised arm ------------------------------------------------
  dy[, 'MVmat'] <- -(mu + um)*MVmat +
    Aging.Prop*c(0, MVmat[1:(N.ages - 1)]) +
    v_mat * period.birth.rate * sum(States) * 0.5 -   # half born to vax mothers
    omega * MVmat
  
  dy[, 'SVmat'] <- MVmat*omega -
    (mu + um)*SVmat +
    Aging.Prop*c(0, SVmat[1:(N.ages - 1)]) -
    omega_mat*SVmat -
    RR * lambda * SVmat                              # SVmat -> IV1
  
  dy[, 'SV0'] <- theta1*WV1 -
    lambda*SV0 -
    (mu + um)*SV0 +
    Aging.Prop*c(0, SV0[1:(N.ages - 1)]) +
    SVmat*omega_mat
  
  dy[, 'IV1'] <- lambda*SV0 -
    (gamma1 + mu + um)*IV1 +
    Aging.Prop*c(0, IV1[1:(N.ages - 1)]) +
    RR * lambda * SVmat
  
  dy[, 'SV1'] <- gamma1*IV1 -
    sigma1*lambda*SV1 -
    (mu + um)*SV1 +
    Aging.Prop*c(0, SV1[1:(N.ages - 1)]) +
    theta1*WV2 -
    theta1*SV1
  
  dy[, 'WV1'] <- theta1*SV1 -
    sigma1*lambda*WV1 -
    (mu + um)*WV1 +
    Aging.Prop*c(0, WV1[1:(N.ages - 1)]) -
    theta1*WV1
  
  dy[, 'IV2'] <- sigma1*lambda*SV1 -
    gamma2*IV2 -
    (mu + um)*IV2 +
    Aging.Prop*c(0, IV2[1:(N.ages - 1)]) +
    sigma1*lambda*WV1
  
  dy[, 'SV2'] <- gamma2*IV2 -
    sigma2*lambda*SV2 -
    (mu + um)*SV2 +
    Aging.Prop*c(0, SV2[1:(N.ages - 1)]) +
    theta1*WV3 -
    theta1*SV2
  
  dy[, 'WV2'] <- theta1*SV2 -
    sigma2*lambda*WV2 -
    (mu + um)*WV2 +
    Aging.Prop*c(0, WV2[1:(N.ages - 1)]) -
    theta1*WV2
  
  dy[, 'IV3'] <- sigma2*lambda*SV2 -
    (gamma3 + mu + um)*IV3 +
    Aging.Prop*c(0, IV3[1:(N.ages - 1)]) +
    sigma2*lambda*WV2
  
  dy[, 'SV3'] <- gamma3*IV3 +
    gamma4*IV4 -
    sigma3*lambda*SV3 -
    (mu + um)*SV3 +
    Aging.Prop*c(0, SV3[1:(N.ages - 1)]) -
    theta1*SV3 +
    Vmat*omega_mat
  
  dy[, 'WV3'] <- theta1*SV3 -
    sigma3*lambda*WV3 -
    (mu + um)*WV3 +
    Aging.Prop*c(0, WV3[1:(N.ages - 1)]) -
    theta1*WV3
  
  dy[, 'IV4'] <- sigma3*lambda*SV3 +
    sigma3*lambda*WV3 -
    gamma4*IV4 -
    (mu + um)*IV4 +
    Aging.Prop*c(0, IV4[1:(N.ages - 1)]) +
    rrM * sigma3 * lambda * Vmat
  
  dy[, 'Vmat'] <- Aging.Prop*c(0, Vmat[1:(N.ages - 1)]) -
    (mu + um)*Vmat -
    omega_mat*Vmat -
    rrM * sigma3 * lambda * Vmat +
    S3_to_Vmat +                      # S3 -> Vmat
    W3_to_Vmat                        # W3 -> Vmat
  
  # ------------ unimmunised arm ---------------------------------------------
  dy[, 'MU'] <- v_mat * period.birth.rate * sum(States) * 0.5 -
    (omega + (mu + um))*MU +
    Aging.Prop*c(0, MU[1:(N.ages - 1)])
  
  dy[, 'SU0'] <- theta1*WU1 +
    omega*MU -
    lambda*SU0 -
    (mu + um)*SU0 +
    Aging.Prop*c(0, SU0[1:(N.ages - 1)])
  
  dy[, 'IU1'] <- lambda*SU0 -
    (gamma1 + mu + um)*IU1 +
    Aging.Prop*c(0, IU1[1:(N.ages - 1)])
  
  dy[, 'SU1'] <- gamma1*IU1 -
    sigma1*lambda*SU1 -
    (mu + um)*SU1 +
    Aging.Prop*c(0, SU1[1:(N.ages - 1)]) +
    theta1*WU2 -
    theta1*SU1
  
  dy[, 'WU1'] <- theta1*SU1 -
    sigma1*lambda*WU1 -
    (mu + um)*WU1 +
    Aging.Prop*c(0, WU1[1:(N.ages - 1)]) -
    theta1*WU1
  
  dy[, 'IU2'] <- sigma1*lambda*SU1 -
    gamma2*IU2 -
    (mu + um)*IU2 +
    Aging.Prop*c(0, IU2[1:(N.ages - 1)]) +
    sigma1*lambda*WU1
  
  dy[, 'SU2'] <- gamma2*IU2 -
    sigma2*lambda*SU2 -
    (mu + um)*SU2 +
    Aging.Prop*c(0, SU2[1:(N.ages - 1)]) +
    theta1*WU3 -
    theta1*SU2
  
  dy[, 'WU2'] <- theta1*SU2 -
    sigma2*lambda*WU2 -
    (mu + um)*WU2 +
    Aging.Prop*c(0, WU2[1:(N.ages - 1)]) -
    theta1*WU2
  
  dy[, 'IU3'] <- sigma2*lambda*SU2 -
    (gamma3 + mu + um)*IU3 +
    Aging.Prop*c(0, IU3[1:(N.ages - 1)]) +
    sigma2*lambda*WU2
  
  dy[, 'SU3'] <- gamma3*IU3 +
    gamma4*IU4 -
    sigma3*lambda*SU3 -
    (mu + um)*SU3 +
    Aging.Prop*c(0, SU3[1:(N.ages - 1)]) -
    theta1*SU3 +
    S3_to_SU3                         # S3 -> SU3
  
  dy[, 'WU3'] <- theta1*SU3 -
    sigma3*lambda*WU3 -
    (mu + um)*WU3 +
    Aging.Prop*c(0, WU3[1:(N.ages - 1)]) -
    theta1*WU3 +
    W3_to_WU3                         # W3 -> WU3
  
  dy[, 'IU4'] <- sigma3*lambda*SU3 +
    sigma3*lambda*WU3 -
    gamma4*IU4 -
    (mu + um)*IU4 +
    Aging.Prop*c(0, IU4[1:(N.ages - 1)])
  
  list(as.vector(dy))
}


# ------------------------------------------------------------------------------
# Checks worth running once
# ------------------------------------------------------------------------------
if (FALSE) {
  
  # (1) every column of yinit.matrix must be assigned a derivative, or the
  #     NA initialisation propagates through the whole solution
  assigned <- c("M","S0","I1","S1","I2","S2","I3","S3","I4","W1","W2","W3",
                "MU","SU0","IU1","SU1","IU2","SU2","IU3","SU3","IU4",
                "WU1","WU2","WU3",
                "MVmat","SVmat","Vmat","SV0","IV1","SV1","IV2","SV2","IV3",
                "SV3","IV4","WV1","WV2","WV3")
  setdiff(colnames(parms.test$yinit.matrix), assigned)   # both should be
  setdiff(assigned, colnames(parms.test$yinit.matrix))   # character(0)
  
  d1 <- model_cc_mat(1, yinit.vector, parms.test)[[1]]
  sum(is.na(d1))                                          # must be 0
  
  # (2) mothers moved == tracked pregnancies (was 2x before this change)
  #     Compare the Vmat inflow against MVmat inflow in a campaign month.
  #     Vmat should gain 0.5*T and MVmat should gain 0.5*T -- equal.
  
  # (3) no negative compartments
  results <- deSolve::ode(yinit.vector, run_times, model_cc_mat, parms.test)
  min(results[, -1])                                      # >= 0
  
  # (4) the 20-39Y class should no longer drain at 2x
  plot(results[, 1], results[, grep("20-39Y S3", colnames(results))],
       type = "l", xlab = "month", ylab = "S3, 20-39Y")
}





##########################################################################################
# old model using 0.5 share between S3 and W3 outflow of vaccinated pregnant mothers 
##########################################################################################

# model_cc_mat <- function(t, y, parms, time.step = "month"){
#   
#   
#   States <- array(y, dim=dim(parms$yinit.matrix))
#   dimnames(States) <- dimnames(parms$yinit.matrix)
#   
#   
#   if(parms$time.step =='month'){
#     period=12
#     length.step=30.44 
#   }else if(parms$time.step=='week'){
#     period=52.1775
#     length.step=7 
#   }
#   
#   
#   omega = 1/(parms$DurationMatImmunityDays/length.step) # waning rate of maternal protection after birth
#   omega_mat = 1/(parms$DurationAbrImmunityDays/length.step) # waning rate of MV in both infants (SVmat-->SV0) and vaccinated mothers (Vmat-->SV3)
#   RR = parms$RR[t] # relative risk of infection after vaccination (for infants); equal to 1 before vaccine introduction
#   rrM = parms$rrM[t] # relative risk of infection after vaccination (for vaxed moms)
#   
#   # vaccine indicator and coverage
#   V = parms$VacPro[t] # VacPro is an indicator for mat doses (binary, =1 during RSV season)
#   # all infants born between October-March will receive a dose at birth (routine birth dose)
#   c_mat <- parms$cover_mat # coverage of birth dose (set to 1 for model fitting)
#   v_mat <- V * c_mat
#   
#   N.ages <- nrow(States)
#   
#   # birth, death, aging in and aging out 
#   period.birth.rate <- log(parms$PerCapitaBirthsYear[t,]+1)/period # birth rate 
#   um = parms$um # death rate
#   mu = 1/parms$WidthAgeClassMonth # aging from current age group to the next age group 
#   if(parms$time.step =='week'){mu = 1/(WidthAgeClassMonth*4.345)}
#   Aging.Prop <- c(0, mu[1:(N.ages-1)]) # aging from last age group to the current age group
#   
#   # recovery rate
#   gamma1 = parms$gamma1  
#   gamma2 = parms$gamma2 
#   gamma3 = parms$gamma3
#   gamma4 = parms$gamma4  
#   
#   # relative risk of infection
#   sigma1 <- parms$sigma1
#   sigma2 <- parms$sigma2
#   sigma3 <- parms$sigma3
#   
#   
#   # rate of waning back to higher susceptible state (default: allows waning)
#   if(parms$AllowWaning=='Yes'){ 
#     theta1 = 1/(parms$dur.immunity1/length.step) 
#   }else{
#     theta1= 0 
#   }
#   
#   
#   # Pull out the initial states for the model as vectors
#   
#   # status quo compartments 
#   M  <-  States[,'M']   
#   
#   S0 <-  States[,'S0']
#   I1 <-  States[,'I1']
#   
#   S1 <-  States[,'S1']
#   I2 <-  States[,'I2']
#   
#   S2 <-  States[,'S2']
#   I3 <-  States[,'I3']
#   
#   S3 <-  States[,'S3']
#   I4 <-  States[,'I4']
#   
#   W1 <-  States[,'W1'] 
#   W2 <-  States[,'W2']
#   W3 <-  States[,'W3']
# 
#   # unimmunized arm
#   MU  <-  States[,'MU']   
#   
#   SU0 <-  States[,'SU0']
#   IU1 <-  States[,'IU1']
#   
#   SU1 <-  States[,'SU1']
#   IU2 <-  States[,'IU2']
#   
#   SU2 <-  States[,'SU2']
#   IU3 <-  States[,'IU3']
#   
#   SU3 <-  States[,'SU3']
#   IU4 <-  States[,'IU4']
#   
#   WU1 <-  States[,'WU1'] 
#   WU2 <-  States[,'WU2']
#   WU3 <-  States[,'WU3']
#   
#   # immunized arm
#   MVmat <- States[,'MVmat'] # newborns born to vaccinated mothers
#   SVmat <- States[,"SVmat"] # newborns born to vaccinated mothers after waning of maternal immunity and cocooning effects (omega)
#   Vmat <- States[,"Vmat"] # vaccinated mothers during pregnancy 
#   
#   SV0 <-  States[,'SV0']
#   IV1 <-  States[,'IV1']
#   
#   SV1 <-  States[,'SV1']
#   IV2 <-  States[,'IV2']
#   
#   SV2 <-  States[,'SV2']
#   IV3 <-  States[,'IV3']
#   
#   SV3 <-  States[,'SV3']
#   IV4 <-  States[,'IV4']
#   
#   WV1 <-  States[,'WV1'] 
#   WV2 <-  States[,'WV2']
#   WV3 <-  States[,'WV3']
#   
#   
#   # force of infection 
#   seasonal.txn <- (1 + parms$b1*cos(2*pi*(t-parms$phi*period)/period)) # seasonality waves
#   transmission_unittime <-  parms$baseline.txn.rate/(parms$dur.days1/length.step) # baseline.txn.rate is the probability of transmission given contact per capita
#   beta = transmission_unittime*parms$c2
#   beta_a_i <- seasonal.txn * beta/(sum(States)^parms$q) # if q = 1: frequency-dependent; q = 0, density-dependent
#   
#   # force of infection 
#   infectiousN <- I1 + parms$rho1*I2 + parms$rho2*I3 + parms$rho2*I4 +
#                  IV1 + parms$rho1*IV2 + parms$rho2*IV3 + parms$rho2*IV4  +
#                  IU1 + parms$rho1*IU2 + parms$rho2*IU3 + parms$rho2*IU4 
#   
#   
#   lambda <- infectiousN %*% beta_a_i 
#   lambda <- as.vector(lambda)
#   
#   
#   dy <- matrix(NA, nrow=N.ages, ncol=ncol(States))
#   colnames(dy) <- colnames(States)
#   
#   
#   # ------------ status quo (before / after one time immunization) ------------------
#   dy[,'M'] <- (1-v_mat) * period.birth.rate*sum(States) - 
#     (omega+(mu+um))*M + 
#     Aging.Prop*c(0, M[1:(N.ages-1)])
#   
#   dy[,'S0'] <- theta1*W1 +
#     omega*M -
#     lambda*S0 -
#     (mu + um)*S0 + 
#     Aging.Prop*c(0,S0[1:(N.ages-1)]) 
#   
#   dy[,'I1'] <-   lambda*S0 - 
#     (gamma1 + mu + um)*I1 + 
#     Aging.Prop*c(0,I1[1:(N.ages-1)]) 
#   
#   dy[,'S1'] <- gamma1*I1 - 
#     sigma1*lambda*S1 - 
#     (mu+um)*S1 + 
#     Aging.Prop*c(0,S1[1:(N.ages-1)]) +
#     theta1*W2-
#     theta1*S1 
#   
#   dy[,'W1'] <-  theta1*S1-
#     sigma1*lambda*W1 -
#     (mu + um)*W1 + 
#     Aging.Prop*c(0,W1[1:(N.ages-1)])-
#     theta1*W1 
#   
#   dy[,'I2'] <- sigma1*lambda*S1 - 
#     gamma2*I2-
#     (mu + um)*I2 + 
#     Aging.Prop*c(0,I2[1:(N.ages-1)]) +
#     sigma1*lambda*W1
#   
#   dy[,'S2'] <- gamma2*I2 - 
#     sigma2*lambda*S2 -
#     (mu+um)*S2 + 
#     Aging.Prop*c(0,S2[1:(N.ages-1)])+
#     theta1*W3 -
#     theta1*S2 
#   
#   dy[,'W2'] <-  theta1*S2-
#     sigma2*lambda*W2 -
#     (mu + um)*W2 + 
#     Aging.Prop*c(0,W2[1:(N.ages-1)])-
#     theta1*W2 
#   
#   dy[,'I3'] <- sigma2*lambda*S2 -
#     (gamma3 + mu+um)*I3 +  
#     Aging.Prop*c(0,I3[1:(N.ages-1)]) +
#     sigma2*lambda*W2
#   
#   dy[,'S3'] <- gamma3*I3 +  
#     gamma4*I4 -
#     sigma3*lambda*S3 -
#     (mu + um)*S3 + 
#     Aging.Prop*c(0,S3[1:(N.ages-1)]) -
#     theta1*S3 - 
#     c(rep(0,10), (v_mat * 0.5) * period.birth.rate[1]*sum(States),rep(0,2)) - # (S3-->Vmat) half of pregnant mothers (aged 20-39) are vaccinated and enter the vax arm
#     c(rep(0,10), (v_mat * 0.5) * period.birth.rate[1]*sum(States),rep(0,2))   # (S3-->SU3) half of pregnant mothers (aged 20-39) are not vaccinated and enter the unvax arm
#     
#   
#   dy[,'W3'] <-  theta1*S3 -
#     sigma3*lambda*W3 -
#     (mu + um)*W3 + 
#     Aging.Prop*c(0,W3[1:(N.ages-1)]) -
#     theta1*W3 -
#     c(rep(0,10), (v_mat * 0.5) * period.birth.rate[1]*sum(States),rep(0,2)) - # W3 --> Vmat 
#     c(rep(0,10), (v_mat * 0.5) * period.birth.rate[1]*sum(States),rep(0,2))   # W3 --> WU3
#   
#   
#   dy[,'I4'] <- sigma3*lambda*S3 +
#     sigma3*lambda*W3 - 
#     gamma4*I4 - 
#     (mu + um)*I4 + 
#     Aging.Prop*c(0,I4[1:(N.ages-1)]) 
#   
#   
#   
#   # ------------ immunized arm (during the month of one-time vax)  ------------------
#   dy[,'MVmat'] <- 
#     - (mu + um) * MVmat +
#     Aging.Prop*c(0, MVmat[1:(N.ages-1)]) + 
#     v_mat * period.birth.rate*sum(States) * 0.5 - # half born to the vaccinated mothers, entering immunized arm
#     omega * MVmat
# 
#   dy[,'SVmat'] <- 
#     MVmat * omega - 
#     (mu + um) * SVmat +
#     Aging.Prop*c(0,SVmat[1:(N.ages-1)]) -
#     omega_mat * SVmat - 
#     RR * lambda * SVmat # SVmat --> IV1
#   
#   
#   dy[,'SV0'] <- theta1*WV1  -
#     lambda*SV0 -
#     (mu + um)*SV0 + 
#     Aging.Prop*c(0,SV0[1:(N.ages-1)]) +
#     SVmat *omega_mat 
#   
#   dy[,'IV1'] <-   lambda*SV0 - 
#     (gamma1 + mu + um)*IV1 + 
#     Aging.Prop*c(0,IV1[1:(N.ages-1)]) +
#     RR * lambda * SVmat
#   
#   dy[,'SV1'] <- gamma1*IV1 - 
#     sigma1*lambda*SV1 - 
#     (mu+um)*SV1 + 
#     Aging.Prop*c(0,SV1[1:(N.ages-1)]) +
#     theta1*WV2-
#     theta1*SV1 
#   
#   dy[,'WV1'] <-  theta1*SV1-
#     sigma1*lambda*WV1 -
#     (mu + um)*WV1 + 
#     Aging.Prop*c(0,WV1[1:(N.ages-1)])-
#     theta1*WV1 
#   
#   dy[,'IV2'] <- sigma1*lambda*SV1 - 
#     gamma2*IV2-
#     (mu + um)*IV2 + 
#     Aging.Prop*c(0,IV2[1:(N.ages-1)]) +
#     sigma1*lambda*WV1
#   
#   dy[,'SV2'] <- gamma2*IV2 - 
#     sigma2*lambda*SV2 -
#     (mu+um)*SV2 + 
#     Aging.Prop*c(0,SV2[1:(N.ages-1)])+
#     theta1*WV3 -
#     theta1*SV2 
#   
#   dy[,'WV2'] <-  theta1*SV2-
#     sigma2*lambda*WV2 -
#     (mu + um)*WV2 + 
#     Aging.Prop*c(0,WV2[1:(N.ages-1)])-
#     theta1*WV2 
#   
#   dy[,'IV3'] <- sigma2*lambda*SV2 -
#     (gamma3 + mu+um)*IV3 +  
#     Aging.Prop*c(0,IV3[1:(N.ages-1)]) +
#     sigma2*lambda*WV2
#   
#   dy[,'SV3'] <- gamma3*IV3 +  
#     gamma4*IV4 -
#     sigma3*lambda*SV3 -
#     (mu + um)*SV3 + 
#     Aging.Prop*c(0,SV3[1:(N.ages-1)]) -
#     theta1*SV3 +
#     Vmat * omega_mat
#   
#   dy[,'WV3'] <-  theta1*SV3 -
#     sigma3*lambda*WV3 -
#     (mu + um)*WV3 + 
#     Aging.Prop*c(0,WV3[1:(N.ages-1)]) -
#     theta1*WV3 
#   
#   dy[,'IV4'] <- sigma3*lambda*SV3 +
#     sigma3*lambda*WV3 - 
#     gamma4*IV4 - 
#     (mu + um)*IV4 + 
#     Aging.Prop*c(0,IV4[1:(N.ages-1)])  +
#     rrM * sigma3 * lambda * Vmat
#   
#   dy[,'Vmat'] <- 
#     Aging.Prop*c(0,Vmat[1:(N.ages-1)]) - 
#     (mu + um)*Vmat -
#     omega_mat * Vmat -
#     rrM * sigma3 * lambda * Vmat +
#     c(rep(0,10), (v_mat * 0.5) * period.birth.rate[1]*sum(States),rep(0,2)) + # W3 --> Vmat
#     c(rep(0,10), (v_mat * 0.5) * period.birth.rate[1]*sum(States),rep(0,2))   # S3 --> Vmat
# 
#   
#   
#   # ------------ unimmunized arm (during the month of one-time vax) ------------------
#   dy[,'MU'] <- v_mat * period.birth.rate*sum(States) * 0.5 -  # half of newborns (born to unvaccinated mothers) entering unimmunized arm
#     (omega+(mu+um))*MU + 
#     Aging.Prop*c(0,MU[1:(N.ages-1)])    
#   
#   dy[,'SU0'] <- theta1*WU1 +
#     omega*MU -
#     lambda*SU0 -
#     (mu + um)*SU0 + 
#     Aging.Prop*c(0,SU0[1:(N.ages-1)]) 
# 
#   dy[,'IU1'] <-   lambda*SU0 - 
#     (gamma1 + mu + um)*IU1 + 
#     Aging.Prop*c(0,IU1[1:(N.ages-1)]) 
#   
#   dy[,'SU1'] <- gamma1*IU1 - 
#     sigma1*lambda*SU1 - 
#     (mu+um)*SU1 + 
#     Aging.Prop*c(0,SU1[1:(N.ages-1)]) +
#     theta1*WU2-
#     theta1*SU1 
#   
#   dy[,'WU1'] <-  theta1*SU1-
#     sigma1*lambda*WU1 -
#     (mu + um)*WU1 + 
#     Aging.Prop*c(0,WU1[1:(N.ages-1)])-
#     theta1*WU1 
#   
#   dy[,'IU2'] <- sigma1*lambda*SU1 - 
#     gamma2*IU2-
#     (mu + um)*IU2 + 
#     Aging.Prop*c(0,IU2[1:(N.ages-1)]) +
#     sigma1*lambda*WU1
#   
#   dy[,'SU2'] <- gamma2*IU2 - 
#     sigma2*lambda*SU2 -
#     (mu+um)*SU2 + 
#     Aging.Prop*c(0,SU2[1:(N.ages-1)])+
#     theta1*WU3 -
#     theta1*SU2 
#   
#   dy[,'WU2'] <-  theta1*SU2-
#     sigma2*lambda*WU2 -
#     (mu + um)*WU2 + 
#     Aging.Prop*c(0,WU2[1:(N.ages-1)])-
#     theta1*WU2 
#   
#   dy[,'IU3'] <- sigma2*lambda*SU2 -
#     (gamma3 + mu+um)*IU3 +  
#     Aging.Prop*c(0,IU3[1:(N.ages-1)]) +
#     sigma2*lambda*WU2
#   
#   dy[,'SU3'] <- gamma3*IU3 +  
#     gamma4*IU4 -
#     sigma3*lambda*SU3 -
#     (mu + um)*SU3 + 
#     Aging.Prop*c(0,SU3[1:(N.ages-1)]) -
#     theta1*SU3  +
#     c(rep(0,10), (v_mat * 0.5) * period.birth.rate[1]*sum(States),rep(0,2)) # S3 --> SU3
#   
#   dy[,'WU3'] <-  theta1*SU3 -
#     sigma3*lambda*WU3 -
#     (mu + um)*WU3 + 
#     Aging.Prop*c(0,WU3[1:(N.ages-1)]) -
#     theta1*WU3 +
#     c(rep(0,10), (v_mat * 0.5) * period.birth.rate[1]*sum(States),rep(0,2)) # W3 --> WU3
#   
#   
#   dy[,'IU4'] <- sigma3*lambda*SU3 +
#     sigma3*lambda*WU3 - 
#     gamma4*IU4 - 
#     (mu + um)*IU4 + 
#     Aging.Prop*c(0,IU4[1:(N.ages-1)]) 
#   
#   
#   derivs <- as.vector(dy)
#   res <- list(derivs)  
#   
#   return(res)
#   
# }
