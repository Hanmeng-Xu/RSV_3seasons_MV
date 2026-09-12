

# ==============================================================================
# separately estimate waning for nirsevimab and MV
# ==============================================================================
mod_wane_adjRSVinc_decrease_separate <- "
model{
  for(j in 1:N_records){
    case_status[j] ~ dbin(pi[j], 1)
    logit(pi[j]) <-
      int +
      inprod(vax_time_group[j, 1:N_groups], beta[1:N_groups]) +
      inprod(M_cf[j, 1:N_cf], beta_cf[1:N_cf]) +
      beta_rsv_activity * rsv_activity[j] +
      beta_bw * birthwt[j] +
      beta_gestage * gestage[j] 
      #beta_age * age_at_test[j]  # only included for nir waning model
  }

  # waning -
  d[1] ~ dnorm(0, 0.0001)
  for(m in 2:N_groups){
    d[m] ~ dnorm(0, tau_d) T(0,)
  }
  beta[1] <- d[1]
  for(m in 2:N_groups){
    beta[m] <- beta[m-1] + d[m]
  }


  # ---- Other priors ----
  int               ~ dnorm(0, 1e-04)
  tau_d             ~ dgamma(0.01, 0.01)
  beta_rsv_activity ~ dnorm(0, 1e-04)
  beta_bw           ~ dnorm(0, 1e-04)
  beta_age          ~ dnorm(0, 1e-04)
  beta_gestage      ~ dnorm(0, 1e-04)
  for(k in 1:N_cf){
    beta_cf[k] ~ dnorm(0, 1e-04)
  }
  
}
"


# ==============================================================================
# spline model: waning of nir & MV estimated separately
# ==============================================================================
mod_wane_spline_separate <- "
model{
  for(j in 1:N_records){
    case_status[j] ~ dbin(pi[j], 1)
    logit(pi[j]) <-
      int +
      inprod(vax_time_group[j, 1:N_groups], beta[1:N_groups]) +
      inprod(M_cf[j, 1:N_cf], beta_cf[1:N_cf]) +
      beta_rsv_activity * rsv_activity[j] +
      beta_bw * birthwt[j] +
      beta_gestage * gestage[j] 
      # beta_age * age_at_test[j]  # only included for nir waning model
  }

  # ---- spline ----
  for(m in 1:N_groups){
    beta[m] <- inprod(B[m, 1:K], gamma[1:K])
  }
  gamma[1] ~ dnorm(0, 1e-04)
  for(k in 2:K){
    gamma[k] ~ dnorm(gamma[k-1], tau_spl)
  }

  # ---- Other priors ----
  int               ~ dnorm(0, 1e-04)
  tau_spl           ~ dgamma(0.01, 0.01)
  beta_rsv_activity ~ dnorm(0, 1e-04)
  beta_bw           ~ dnorm(0, 1e-04)
  beta_age          ~ dnorm(0, 1e-04)
  beta_gestage      ~ dnorm(0, 1e-04)
  for(k in 1:N_cf){
    beta_cf[k] ~ dnorm(0, 1e-04)
  }
  
}
"


# ==============================================================================
# separatelt estimate VE of nir and MV, AR(1) smoothing on each waning curve
# ==============================================================================
mod_wane_adjRSVinc_notrendAR_separate <- "
model{
  for(j in 1:N_records){
    case_status[j] ~ dbin(pi[j], 1)
    logit(pi[j]) <-
      int +
      inprod(vax_time_group[j, 1:N_groups], beta[1:N_groups]) +
      inprod(M_cf[j, 1:N_cf], beta_cf[1:N_cf]) +
      beta_rsv_activity * rsv_activity[j] +
      beta_bw * birthwt[j] +
      beta_gestage * gestage[j] 
      # beta_age * age_at_test[j]  # only included for nir waning model
  }

   # ---- waning: AR(1) ----
  beta[1] ~ dnorm(0, tau_ar * (1 - rho^2))   # stationary marginal
  for(m in 2:N_groups){
    beta[m] ~ dnorm(rho * beta[m-1], tau_ar)
  }

  # ---- Other priors ----
  int               ~ dnorm(0, 1e-04)
  rho               ~ dunif(-1, 1)
  tau_ar            ~ dgamma(0.01, 0.01)
  beta_rsv_activity ~ dnorm(0, 1e-04)
  beta_bw           ~ dnorm(0, 1e-04)
  beta_age          ~ dnorm(0, 1e-04)
  beta_gestage      ~ dnorm(0, 1e-04)
  for(k in 1:N_cf){
    beta_cf[k] ~ dnorm(0, 1e-04)
  }
  
}
"



# ==============================================================================
# separately estiamte VE for nir and MV
# AR(1) + linear structure for stabilization, on each waning curve
# ==============================================================================
mod_wane_adjRSVinc_notrendAR_wLinear_seperate <- "
model{
  for(j in 1:N_records){
    case_status[j] ~ dbin(pi[j], 1)
    logit(pi[j]) <-
      int +
      inprod(vax_time_group[j, 1:N_groups], beta[1:N_groups]) +
      inprod(M_cf[j, 1:N_cf], beta_cf[1:N_cf]) +
      beta_rsv_activity * rsv_activity[j] +
      beta_bw * birthwt[j] +
      beta_gestage * gestage[j] 
      # beta_age * age_at_test[j]  # only included for nir waning model
  }

  # ---- waning: linear trend + AR(1) deviations ----
  beta[1] ~ dnorm( mu0 + mu1*0 , tau_ar * (1 - rho^2))
  for(m in 2:N_groups){
    beta[m] ~ dnorm( mu0 + mu1*(m-1) + rho * beta[m-1], tau_ar)
  }

  # ---- Other priors ----
  int               ~ dnorm(0, 1e-04)
  # rho             ~ dunif(-1, 1)
  rho               ~ dunif(0, 1) # suppresses the sawtooth.
  tau_ar            ~ dgamma(0.01, 0.01)
  beta_rsv_activity ~ dnorm(0, 1e-04)
  beta_bw           ~ dnorm(0, 1e-04)
  beta_age          ~ dnorm(0, 1e-04)
  beta_gestage      ~ dnorm(0, 1e-04)
  for(k in 1:N_cf){
    beta_cf[k] ~ dnorm(0, 1e-04)
  }

  mu0 ~ dnorm(0, 1e-04)
  mu1 ~ dnorm(0, 1e-04)
}
"


################################################################################
# age as modifier, modify binary MV effect
################################################################################
mod_wane_agemodi_decrease <- "
model{
  for(j in 1:N_records){
    case_status[j] ~ dbin(pi[j], 1)
    logit(pi[j]) <-
      # int + # as we are not using a reference group, if we have a int here, the int will be not estimable
      inprod(age_group[j, 1:N_groups], alpha[1:N_groups]) +  # this is just like a confounder for categorical age variable (or can understand as baseline unvaxed )
      inprod(age_group[j, 1:N_groups], beta[1:N_groups]) * mv[j] +
      inprod(M_cf[j, 1:N_cf], beta_cf[1:N_cf]) + 
      beta_rsv_activity * rsv_activity[j] + 
      beta_bw * birthwt[j] +
      beta_gestage * gestage[j] 
      
  }

  # ---- waning ----
  d[1] ~ dnorm(0, 0.0001)
  for(m in 2:N_groups){
    d[m] ~ dnorm(0, tau_d) T(0,)
  }
  beta[1] <- d[1]
  for(m in 2:N_groups){
    beta[m] <- beta[m-1] + d[m]
  }
  
  # ---- age main effect, ALL infants, bin 1 fixed at 0 as reference ----
  alpha[1] <- 0
  for(m in 2:N_groups){
    alpha[m] ~ dnorm(alpha[m-1], tau_a)
  }



  # ---- Other priors ----
  int               ~ dnorm(0, 1e-04)
  tau_d             ~ dgamma(0.01, 0.01)
  tau_a             ~ dgamma(0.01, 0.01)
  beta_rsv_activity ~ dnorm(0, 1e-04)
  beta_bw           ~ dnorm(0, 1e-04)
  beta_gestage      ~ dnorm(0, 1e-04)
  for(k in 1:N_cf){
    beta_cf[k] ~ dnorm(0, 1e-04)
  }
  
}
"



mod_wane_agemodi_spline <- "
model{
  for(j in 1:N_records){
    case_status[j] ~ dbin(pi[j], 1)
    logit(pi[j]) <-
      # int + # as we are not using a reference group, if we have a int here, the int will be not estimable
      inprod(age_group[j, 1:N_groups], alpha[1:N_groups]) +  # this is just like a confounder for categorical age variable (or can understand as baseline for unvaxed )
      inprod(age_group[j, 1:N_groups], beta[1:N_groups]) * mv[j] +
      inprod(M_cf[j, 1:N_cf], beta_cf[1:N_cf]) + 
      beta_rsv_activity * rsv_activity[j] +
      beta_bw * birthwt[j] +
      beta_gestage * gestage[j]
      
  }

  # spline structure 
  for(m in 1:N_groups){
    beta[m] <- inprod(B[m, 1:K], gamma[1:K])
  }
  gamma[1] ~ dnorm(0, 1e-04)
  for(k in 2:K){
    gamma[k] ~ dnorm(gamma[k-1], tau_spl) # first order RW on the gamma
  }
  # 1st-order RW: the next slope (coeff) should be similar to the previous one
  # 2nd-order RW: the change in slope should be similar
  
  
  # age main effect, ALL infants, bin 1 fixed at 0 as reference 
  alpha[1] <- 0
  for(m in 2:N_groups){
    alpha[m] ~ dnorm(alpha[m-1], tau_a)
  }



  # priors
  int               ~ dnorm(0, 1e-04)
  tau_spl             ~ dgamma(0.01, 0.01)
  tau_a  ~ dgamma(0.01, 0.01)
  beta_rsv_activity ~ dnorm(0, 1e-04)
  beta_bw ~ dnorm(0, 1e-04)
  beta_gestage ~ dnorm(0, 1e-04)
  for(k in 1:N_cf){
    beta_cf[k] ~ dnorm(0, 1e-04)
  }
  
}
"







mod_wane_agemodi_notrendAR <- "
model{
  for(j in 1:N_records){
    case_status[j] ~ dbin(pi[j], 1)
    logit(pi[j]) <-
      # int + # as we are not using a reference group, if we have a int here, the int will be not estimable
      inprod(age_group[j, 1:N_groups], alpha[1:N_groups]) +  # this is just like a confounder for categorical age variable (or can understand as baseline unvaxed )
      inprod(age_group[j, 1:N_groups], beta[1:N_groups]) * mv[j] +
      inprod(M_cf[j, 1:N_cf], beta_cf[1:N_cf]) + 
      beta_rsv_activity * rsv_activity[j] +
      beta_bw * birthwt[j] +
      beta_gestage * gestage[j]
  }
  
  # ---- age main effect, ALL infants, bin 1 fixed at 0 as reference ----
  alpha[1] <- 0
  for(m in 2:N_groups){
    alpha[m] ~ dnorm(alpha[m-1], tau_a)
  }
 

  # ---- waning: AR(1) ----
  beta[1] ~ dnorm(0, tau_ar * (1 - rho^2))   # stationary marginal
  for(m in 2:N_groups){
    beta[m] ~ dnorm(rho * beta[m-1], tau_ar)
  }

  # ---- Other priors ----
  int               ~ dnorm(0, 1e-04)
  rho               ~ dunif(-1, 1)
  tau_a            ~ dgamma(0.01, 0.01)
  tau_ar            ~ dgamma(0.01, 0.01)
  beta_rsv_activity ~ dnorm(0, 1e-04)
  beta_bw           ~ dnorm(0, 1e-04)
  beta_gestage      ~ dnorm(0, 1e-04)
  for(k in 1:N_cf){
    beta_cf[k] ~ dnorm(0, 1e-04)
  }

}
"






mod_wane_agemodi_notrendAR_wLinear <- "
model{
  for(j in 1:N_records){
    case_status[j] ~ dbin(pi[j], 1)
    logit(pi[j]) <-
      # int + # as we are not using a reference group, if we have a int here, the int will be not estimable
      inprod(age_group[j, 1:N_groups], alpha[1:N_groups]) +  # this is just like a confounder for categorical age variable (or can understand as baseline unvaxed )
      inprod(age_group[j, 1:N_groups], beta[1:N_groups]) * mv[j] +
      inprod(M_cf[j, 1:N_cf], beta_cf[1:N_cf]) + 
      beta_rsv_activity * rsv_activity[j] +
      beta_bw * birthwt[j] +
      beta_gestage * gestage[j]
  }
  
  # ---- age main effect, ALL infants, bin 1 fixed at 0 as reference ----
  alpha[1] <- 0
  for(m in 2:N_groups){
    alpha[m] ~ dnorm(alpha[m-1], tau_a)
  }
 

 
  
  # ---- waning: linear trend + AR(1) deviations ----
  beta[1] ~ dnorm( mu0 + mu1*0 , tau_ar * (1 - rho^2))
  # beta[1] ~ dnorm( mu1*0 , tau_ar * (1 - rho^2))
  for(m in 2:N_groups){
    beta[m] ~ dnorm( mu0 + mu1*(m-1) + rho * beta[m-1], tau_ar)
    # beta[m] ~ dnorm( mu1*(m-1) + rho * beta[m-1], tau_ar)
  }

  # ---- Other priors ----
  int               ~ dnorm(0, 1e-04)
  rho               ~ dunif(0, 1) # suppresses the sawtooth.
  tau_ar            ~ dgamma(0.01, 0.01)
  tau_a            ~ dgamma(0.01, 0.01)
  beta_rsv_activity ~ dnorm(0, 1e-04)
  beta_bw           ~ dnorm(0, 1e-04)
  beta_gestage      ~ dnorm(0, 1e-04)
  for(k in 1:N_cf){
    beta_cf[k] ~ dnorm(0, 1e-04)
  }

  mu0 ~ dnorm(0, 1e-04)
  mu1 ~ dnorm(0, 1e-04)

}
"


