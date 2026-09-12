library(dplyr)
library(ggplot2)
library(deSolve)
library(lubridate)
library(tidyr)
library(lhs)
library(MASS)
library(cowplot)

source("code/fit_durMat/model_cc_mat.R") 


# -------------------------------------------------------------------------------
tofitSES <- "medium" # using transmission parameters from "unstratified", "low", "medium", "high"
# -------------------------------------------------------------------------------

if(tofitSES == "unstratified"){
  birth <-  readRDS(paste0("../rsv_CEA/data/transmission_model/BrNY_pred.rds"))
  yinit <- readRDS(paste0("../rsv_CEA/data/transmission_model/yinit.rds")) 
}else{
  birth <- readRDS(paste0("../rsv_CEA/data/transmission_model/BrNY_pred_CDC_wonder_bySES.rds"))[[paste0("BrNY.", tofitSES, "SES")]]
  yinit <-  readRDS(paste0("../rsv_CEA/data/transmission_model/yinit.NY.byses.1981.rds"))[[paste0("yinit.", tofitSES)]]
}

# birth <-  readRDS(paste0("data_and_parms/NY/BrNY_pred.rds"))
# yinit <- readRDS(paste0("data_and_parms/NY/yinit.rds")) 
contactUSAinfant <- readRDS( '../rsv_CEA/data/transmission_model/contactmatrix.rds') # contact matrix
c2 <- contactUSAinfant


# names of Age groups 
N_ages <- nrow(yinit) 
agenames <- c("<2m","2-3m","4-5m","6-7m","8-9m","10-11m","1Y","2-4Y","5-9Y","10-19Y","20-39Y","40-59Y","60Y+")
agenames_children <- c("<2m","2-3m","4-5m","6-7m","8-9m","10-11m","1Y","2-4Y")
al <- N_ages

# initialize pop matrix
yinit <- cbind(yinit, 
               W1 = rep(0,13), W2 = rep(0,13), W3 = rep(0,13), 
               MVmat = rep(0,13), SVmat = rep(0, 13), Vmat = rep(0,13), 
               WV1 = rep(0,13), WV2 = rep(0,13), WV3 = rep(0,13), 
               SV0 = rep(0,13), SV1 = rep(0,13), SV2 = rep(0,13), SV3 = rep(0,13), 
               IV1 = rep(0,13), IV2 = rep(0,13), IV3 = rep(0,13), IV4 = rep(0,13), 
               MU = rep(0,13),
               WU1 = rep(0,13), WU2 = rep(0,13), WU3 = rep(0,13), 
               SU0 = rep(0,13), SU1 = rep(0,13), SU2 = rep(0,13), SU3 = rep(0,13), 
               IU1 = rep(0,13), IU2 = rep(0,13), IU3 = rep(0,13), IU4 = rep(0,13))


rownames(yinit) <- agenames
yinit.vector <- as.vector(yinit)

name.array <- array(NA, dim=dim(yinit))
for(i in 1:dim(name.array)[1]){
  for(j in 1:dim(name.array)[2]){
    name.array[i,j] <- paste(dimnames(yinit)[[1]][i],dimnames(yinit)[[2]][j]  )
  }
}

name.vector <- as.vector(name.array)
names(yinit.vector) <- name.vector



#Relative infectiousness for 2nd and subsequent infections
rho1 = 0.75
rho2 = 0.51

# duration of infectiousness (months)
dur.days1 <- 10 #days
dur.days2 <- 7 #days
dur.days3 <- 5 #days

WidthAgeClassMonth = c(rep(2,times=6), 12,12*3,  60, 120, 240, 240, 240)  #Aging rate=1/width age class (months) 


#um : calibrated this parameter so we can reproduce the population growth
# recalibrated for NY
um = 0.000015


#Birth rate (births/person/YEAR
#Matrix: T rows, N_ages columns; columns 2:N_ages all 0s
PerCapitaBirthsYear=birth

#Relaive risk of infection following 1st, 2nd, 3rd+ infections
sigma1=0.76
sigma2=0.6
sigma3=0.4

#Relaive risk of hospitalizations following 1st, 2nd infections
hosp1= c(0.082251950, 0.047607003, 0.026990367, 0.016398471, 0.012376537, 0.012007695, 0.009609524,0.006999893, 0.001000000, 0.001000000, 0.001000000, 0.001000000, 0.001000000)
hosp2= 0.4 * hosp1
hosp3= rep(0,13)


length.step = 30.44 #days per month

# recovery rate following te 1st, 2nd, 3rd infection
gamma1= 1/(dur.days1/length.step)  #converts 1/days to 1/lenth.step
gamma2= 1/(dur.days2/length.step)  
gamma3= 1/(dur.days3/length.step)  


# recording fraction of hospitalizations (from burden project) by age 
# report_ratio <- c(0.76, 0.68, 0.62, 0.68, 0.67, 0.57, 0.50, 0.34, 0.12, 0.059, 0.023, 0.04, 0.04)


# Estimated transmission parameters from model
fitted_param_set <- readRDS(paste0("../rsv_CEA/data/fitted_parms/parm_STAN_10000LHS.rds")) %>% 
  filter(ses == tofitSES) %>%
  dplyr::select(beta, b1, phi, omega, theta1)

beta.median <- median(fitted_param_set$beta) # for test run and test fitting
b1.median <- median(fitted_param_set$b1)
phi.median <- median(fitted_param_set$phi)
omega.median <- median(fitted_param_set$omega)
theta1.median <- median(fitted_param_set$theta1)


# relative risk of infection post vaccination among infants and moms
# (if fixed: relative_risk_infants <- 1 - 0.849  # relative_risk_moms <- 1 - 0.9, ref: npj paper)
RR_rrM <- readRDS("../rsv_CEA/data/risk/RR_rrM.rds")
RR.median <- median(RR_rrM$RR_infants)
rrM.median <- median(RR_rrM$rrM_moms)




# time setup
start_time = 1 
tmax = (2035 - 1981) * 12 # 1981-7-1 ~ 2035-6-1
run_times <- seq(start_time, tmax, by = 1) 


parmset <- list(
  WidthAgeClassMonth=WidthAgeClassMonth,
  um=um,
  rho1=rho1,
  rho2=rho2,
  dur.days1=dur.days1,
  dur.days2=dur.days2,
  dur.days3=dur.days3,
  q=1,
  AllowWaning = 'Yes',
  c2=contactUSAinfant,
  sigma1=sigma1,
  sigma2=sigma2,
  sigma3=sigma3,
  gamma1=gamma1,
  gamma2=gamma2,
  gamma3=gamma3,
  gamma4=gamma3,
  time.step = 'month',
  yinit.matrix = yinit,
  PerCapitaBirthsYear = birth
)



# ------------------------------------------------------------------------------
# test run 
# ------------------------------------------------------------------------------
# preset duration of protection from mv
dur_mat_test <- 200

parms.test <- c(
  parmset, 
  list(
    RR = c(rep(1, 504), rep(RR.median, (tmax-504))), # relative risk of infection for infants (SVmat) born to vaccianted mothers, the reduced risk maintains after the vax month
    rrM = c(rep(1, 504), rep(rrM.median, (tmax-504))), # relative risk of infection for vaccinated mothers (Vmat), the reducted risk maintains after the vax months
    baseline.txn.rate = beta.median,
    b1 = b1.median,
    phi = phi.median,
    DurationMatImmunityDays = 1/(omega.median) * length.step, 
    DurationAbrImmunityDays = dur_mat_test,
    dur.immunity1 = (1/theta1.median) * length.step,
    WidthAgeClassMonth = WidthAgeClassMonth,
    cover_mat  = 1,
    # birth dose: one-time administration in October 2023
    VacPro  = c(rep(0,504), # 1981-7-1 ~ 2023-6-1: no vaccine available
                       c(0,0,0,1,0,0,0,0,0,0,0,0),
                       rep(rep(0,12), (tmax - 504)/12-1 ))
  ) 
  
)

# run ode 
results <- ode(y = yinit.vector, 
               t = run_times, 
               func = model_cc_mat, # the function is the same as fixwab model
               parms = parms.test)



# extract compartments 
St <- results[,-1]
M <- results[,grep(' M', colnames(results))]
I1 <- results[,grep(' I1', colnames(results))]
I2 <- results[,grep(' I2', colnames(results))]
I3 <- results[,grep(' I3', colnames(results))]
I4 <- results[,grep(' I4', colnames(results))]
S1 <- results[,grep(' S1', colnames(results))]
S2 <- results[,grep(' S2', colnames(results))]
S3 <- results[,grep(' S3', colnames(results))]
S0 <- results[,grep(' S0', colnames(results))]
W1 <- results[,grep(' W1', colnames(results))]
W2 <- results[,grep(' W2', colnames(results))]
W3 <- results[,grep(' W3', colnames(results))]

IV1 <- results[,grep(' IV1', colnames(results))]
IV2 <- results[,grep(' IV2', colnames(results))]
IV3 <- results[,grep(' IV3', colnames(results))]
IV4 <- results[,grep(' IV4', colnames(results))]
SV1 <- results[,grep(' SV1', colnames(results))]
SV2 <- results[,grep(' SV2', colnames(results))]
SV3 <- results[,grep(' SV3', colnames(results))]
SV0 <- results[,grep(' SV0', colnames(results))]
WV1 <- results[,grep(' WV1', colnames(results))]
WV2 <- results[,grep(' WV2', colnames(results))]
WV3 <- results[,grep(' WV3', colnames(results))]
MVmat <- results[,grep(' MVmat', colnames(results))]
SVmat <- results[,grep(' SVmat', colnames(results))]
Vmat <- results[,grep(' Vmat', colnames(results))]

MU <- results[,grep(' MU', colnames(results))]
IU1 <- results[,grep(' IU1', colnames(results))]
IU2 <- results[,grep(' IU2', colnames(results))]
IU3 <- results[,grep(' IU3', colnames(results))]
IU4 <- results[,grep(' IU4', colnames(results))]
SU1 <- results[,grep(' SU1', colnames(results))]
SU2 <- results[,grep(' SU2', colnames(results))]
SU3 <- results[,grep(' SU3', colnames(results))]
SU0 <- results[,grep(' SU0', colnames(results))]
WU1 <- results[,grep(' WU1', colnames(results))]
WU2 <- results[,grep(' WU2', colnames(results))]
WU3 <- results[,grep(' WU3', colnames(results))]


lambda1 <- matrix(0, nrow = tmax, ncol = N_ages) 
beta <-  parms.test$baseline.txn.rate / (parms.test$dur.days1 / 30.44) * parms.test$c2

for (t in 1:tmax) {
  # force of infection
  lambda1[t, ] <- as.vector((1 + parms.test$b1 * cos(2 * pi * (t - parms.test$phi * 12) / 12)) *
                              ( ( I1[t, ] + parms.test$rho1 * I2[t, ] + parms.test$rho2 * I3[t, ] + parms.test$rho2 * I4[t, ] +
                                    IV1[t, ] + parms.test$rho1 * IV2[t, ] + parms.test$rho2 * IV3[t, ] + parms.test$rho2 * IV4[t, ] + 
                                    IU1[t, ] + parms.test$rho1 * IU2[t, ] + parms.test$rho2 * IU3[t, ] + parms.test$rho2 * IU4[t, ]
                              ) %*% beta ) / sum(St[t, ]))
}

# incidence over time in the immunized and unimmunized arm
inc_immu = matrix(0, nrow = tmax, ncol = al) 
inc_unimmu = matrix(0, nrow = tmax, ncol = al) 
for(i in 1:al){
  inc_immu[, i] =  
    lambda1[,i] * SVmat[,i] * RR.median +
    lambda1[,i] * SV0[,i] +  
    lambda1[,i] * (SV1[,i] + WV1[,i]) * sigma1 +
    lambda1[,i] * (SV2[,i] + WV2[,i]) * sigma2 +
    lambda1[,i] * (SV3[,i] + WV3[,i]) * sigma3 
  inc_unimmu[, i] =  
    lambda1[,i] * SU0[,i] +  
    lambda1[,i] * (SU1[,i] + WU1[,i]) * sigma1 +
    lambda1[,i] * (SU2[,i] + WU2[,i]) * sigma2 +
    lambda1[,i] * (SU3[,i] + WU3[,i]) * sigma3 
}

# VE over time (waning VE) = 1- inc_immu / inc_unimmu (among under 12 months age groups)
VE = 1 - rowSums(inc_immu[, 1:6])/ rowSums(inc_unimmu[, 1:6]) # VE since birth of infants (from vaccination to mother ~ 2-6 weeks)
plot(VE[(509 + 1):(520 + 1)], type = "l") # lagged 4 months (to account for the interval between vax to delivery)




# ------------------------------------------------------------------------------
# fit to the case control data  
# ------------------------------------------------------------------------------
# model_cc <- "Decrease trend imposed" # decrease model has the lowest DIC
model_cc <- "B-spline"


ve.cc <- readRDS("code/fit_durMat/data/ve_wane_sep_mv_nir.rds") %>% # filter(endpoint == "MA RSV infection") %>% 
  mutate(VE_median = VE_median / 100, 
         VE_lb = VE_lb / 100, 
         VE_ub = VE_ub / 100) %>% 
  filter(product == "RSVpreF") %>% 
  filter(model == model_cc)

# plot the estimated waning VE estiamted from the cc study
ve.cc %>% 
  filter(month <= 12) %>%
  ggplot(aes(x = month, y = VE_median)) +
  geom_line() +
  theme_bw() +
  scale_x_continuous(breaks = 1:12, labels = as.character(1:12)) +
  ggtitle("Effectiveness of RSVpreF (against MA RSV infection) estimated from the the case-control study")



# ── Helper: run ODE and return 12-month model VE ─────────────────────────────
run_model_VE <- function(dur_mat) {
  
  parms.test <- c(
    parmset, 
    list(
      RR = c(rep(1, 504), rep( RR.median , (tmax-504))),  
      rrM = c(rep(1, 504), rep(rrM.median, (tmax-504))), 
      baseline.txn.rate = beta.median,
      b1 = b1.median,
      phi = phi.median,
      DurationMatImmunityDays = 1/(omega.median) * length.step, 
      DurationAbrImmunityDays = dur_mat,
      dur.immunity1 = (1/theta1.median) * length.step,
      WidthAgeClassMonth = WidthAgeClassMonth,
      cover_mat  = 1,
      # birth dose: one-time administration in October 2023
      VacPro  = c(rep(0,504), # 1981-7-1 ~ 2023-6-1: no vaccine available
                  c(0,0,0,1,0,0,0,0,0,0,0,0),
                  rep(rep(0,12), (tmax - 504)/12-1 ))
    )
  )
  
  results <- ode(y     = yinit.vector,
                 t     = run_times,
                 func  = model_cc_mat,
                 parms = parms.test)
  
  # ── Extract compartments ───────────────────────────────────────────────────
  St  <- results[, -1]
  I1  <- results[, grep(' I1',  colnames(results))]
  I2  <- results[, grep(' I2',  colnames(results))]
  I3  <- results[, grep(' I3',  colnames(results))]
  I4  <- results[, grep(' I4',  colnames(results))]
  IV1 <- results[, grep(' IV1', colnames(results))]
  IV2 <- results[, grep(' IV2', colnames(results))]
  IV3 <- results[, grep(' IV3', colnames(results))]
  IV4 <- results[, grep(' IV4', colnames(results))]
  SV0 <- results[, grep(' SV0', colnames(results))]
  SV1 <- results[, grep(' SV1', colnames(results))]
  SV2 <- results[, grep(' SV2', colnames(results))]
  SV3 <- results[, grep(' SV3', colnames(results))]
  WV1 <- results[, grep(' WV1', colnames(results))]
  WV2 <- results[, grep(' WV2', colnames(results))]
  WV3 <- results[, grep(' WV3', colnames(results))]
  
  IU1 <- results[, grep(' IU1', colnames(results))]
  IU2 <- results[, grep(' IU2', colnames(results))]
  IU3 <- results[, grep(' IU3', colnames(results))]
  IU4 <- results[, grep(' IU4', colnames(results))]
  SU0 <- results[, grep(' SU0', colnames(results))]
  SU1 <- results[, grep(' SU1', colnames(results))]
  SU2 <- results[, grep(' SU2', colnames(results))]
  SU3 <- results[, grep(' SU3', colnames(results))]
  WU1 <- results[, grep(' WU1', colnames(results))]
  WU2 <- results[, grep(' WU2', colnames(results))]
  WU3 <- results[, grep(' WU3', colnames(results))]
  
  MVmat <- results[,grep(' MVmat', colnames(results))]
  SVmat <- results[,grep(' SVmat', colnames(results))]
  Vmat <-  results[,grep(' Vmat', colnames(results))]
  
  # ── Force of infection ─────────────────────────────────────────────────────
  beta    <- parms.test$baseline.txn.rate / (parms.test$dur.days1 / 30.44) * parms.test$c2
  lambda1 <- matrix(0, nrow = tmax, ncol = N_ages)
  
  for (t in 1:tmax) {
    lambda1[t, ] <- as.vector(
      (1 + parms.test$b1 * cos(2 * pi * (t - parms.test$phi * 12) / 12)) *
        ((I1[t,]  + parms.test$rho1 * I2[t,]  + parms.test$rho2 * I3[t,]  + parms.test$rho2 * I4[t,]  +
            IV1[t,] + parms.test$rho1 * IV2[t,] + parms.test$rho2 * IV3[t,] + parms.test$rho2 * IV4[t,] +
            IU1[t,] + parms.test$rho1 * IU2[t,] + parms.test$rho2 * IU3[t,] + parms.test$rho2 * IU4[t,]
        ) %*% beta) / sum(St[t, ])
    )
  }
  
  # ── Incidence & VE ─────────────────────────────────────────────────────────
  inc_immu   <- matrix(0, nrow = tmax, ncol = al)
  inc_unimmu <- matrix(0, nrow = tmax, ncol = al)
  
  for (i in 1:al) {
    inc_immu[, i] =  
      lambda1[,i] * SVmat[,i] * RR.median +
      lambda1[,i] * SV0[,i] +  
      lambda1[,i] * (SV1[,i] + WV1[,i]) * sigma1 +
      lambda1[,i] * (SV2[,i] + WV2[,i]) * sigma2 +
      lambda1[,i] * (SV3[,i] + WV3[,i]) * sigma3
    inc_unimmu[, i] =  
      lambda1[,i] * SU0[,i] +  
      lambda1[,i] * (SU1[,i] + WU1[,i]) * sigma1 +
      lambda1[,i] * (SU2[,i] + WU2[,i]) * sigma2 +
      lambda1[,i] * (SU3[,i] + WU3[,i]) * sigma3
  }
  
  VE <- 1 - rowSums(inc_immu[,1:6]) / rowSums(inc_unimmu[,1:6]) # only look at infant age groups 
  
  # Return only the 16-month post-vaccination window (Oct 2023 = timestep 509)
  month_idx <- 508 + 1 + 2 + (1:6) # "+ 1": time lag between vax to delivery # "+ 2": from month 3 of the waning curve (as the first 2 months is having increasing trend in spline model)
  return(VE[month_idx])
}


# ── SSE objective ─────────────────────────────────────────────────────────────
sse_objective <- function(log_dur_mat) { # sum of squared errors
  dur_mat  <- exp(log_dur_mat)
  VE_model <- run_model_VE(dur_mat)
  sse      <- sum((VE_model - ve.cc$VE_median[3:8])^2, na.rm = TRUE) ## start fitting from month 3
  cat(sprintf("  dur_mat=%.1f  SSE=%.6f\n", dur_mat, sse))
  return(sse)
}


opt <- optim(
    par      = log(72),          # initial guess: 72 days
    fn       = sse_objective,
    method   = "Brent",          # fast & reliable for 1-D
    lower    = log(1),           # dur_mab >= 1 day
    upper    = log(365),         # dur_mab <= 365 days
    control  = list(maxit = 200)
)



# ------------------------------------------------------------------------------
# profile sampling to estimate uncertainty of dur_mat
# ------------------------------------------------------------------------------
# 0. Observed target + cached model evaluations ─────────────────────────────
ve_obs <- ve.cc %>% arrange(month) %>% slice(3:8) %>% pull(VE_median) ## start fitting from month 3
stopifnot(length(ve_obs) == 6, all(!is.na(ve_obs)))
n <- sum(!is.na(ve_obs))

.cache <- new.env(hash = TRUE)
ve_cached <- function(dur_mat) {
  key <- sprintf("%.4f", dur_mat)
  if (!is.null(.cache[[key]])) return(.cache[[key]])
  out <- tryCatch(run_model_VE(dur_mat), error = function(e) rep(NA_real_, 6))
  assign(key, out, envir = .cache)
  out
}

sse_of <- function(dur_mat) {
  pred <- ve_cached(dur_mat)
  if (all(is.na(pred))) return(NA_real_)
  sum((pred - ve_obs)^2, na.rm = TRUE)
}

# 1. Point estimate 
sse_objective <- function(log_dur_mat) {
  d <- exp(log_dur_mat); s <- sse_of(d)
  cat(sprintf("  dur_mat=%7.1f  SSE=%.6f\n", d, s))
  if (!is.finite(s)) return(1e10)
  s
}

opt <- optim(par = log(72), fn = sse_objective, method = "Brent",
             lower = log(1), upper = log(365), control = list(maxit = 200))

best_dur_mat <- exp(opt$par) # dur_mat = 122 days (based on ve_cc of decrease model)
sse_min      <- opt$value

# 2. Profile: coarse scan, then refine only where it matters 
dev_of <- function(d) n * log(sse_of(d) / sse_min)

coarse <- seq(20, 365, by = 10)
dev_c  <- vapply(coarse, dev_of, numeric(1))

# widen to deviance < 8 so the 3.84 crossing is safely bracketed
keep <- coarse[is.finite(dev_c) & dev_c < 8]
fine <- seq(max(20, min(keep) - 10), min(365, max(keep) + 10), by = 1)

prof <- data.frame(dur_mat = fine)
prof$deviance <- vapply(fine, dev_of, numeric(1))
prof$logLik   <- -0.5 * prof$deviance          # up to an additive constant

# 3. CI by root-finding (robust to a non-monotone profile)
cross <- function(thresh, side) {
  br <- if (side == "left") {
    prof[prof$dur_mat <= best_dur_mat, ][order(-prof$dur_mat[prof$dur_mat <= best_dur_mat]), ]
  } else {
    prof[prof$dur_mat >= best_dur_mat, ][order(prof$dur_mat[prof$dur_mat >= best_dur_mat]), ]
  }
  k <- which(br$deviance > thresh)[1]
  if (is.na(k) || k == 1) return(NA_real_)      # never crosses inside the grid
  uniroot(function(d) dev_of(d) - thresh,
          interval = sort(c(br$dur_mat[k - 1], br$dur_mat[k])),
          tol = 0.1)$root
}

thresh_chisq <- qchisq(0.95, df = 1)                          # 3.841
thresh_F     <- n * log(1 + qf(0.95, 1, n - 1) / (n - 1))     # small-sample correction (fitted to 6 VE data points)

ci <- rbind(
  chisq = c(cross(thresh_chisq, "left"), cross(thresh_chisq, "right")),
  F     = c(cross(thresh_F,     "left"), cross(thresh_F,     "right"))
)
colnames(ci) <- c("lower_95", "upper_95")
print(round(cbind(MLE = best_dur_mat, ci), 1))

# 4. Look at it
plot(prof$dur_mat, prof$deviance, type = "l",
     xlab = "dur_mat (days)", ylab = "profile deviance")
abline(h = c(thresh_chisq, thresh_F), lty = c(2, 3), col = c("grey40", "firebrick"))
abline(v = best_dur_mat, col = "steelblue")

# 5. Samples for downstream propagation 
rel_L <- exp(prof$logLik - max(prof$logLik))
dx    <- c(diff(prof$dur_mat), tail(diff(prof$dur_mat), 1))
mass  <- rel_L * dx; mass <- mass / sum(mass)

set.seed(1)
dur_mat_samples <- approx(cumsum(mass), prof$dur_mat, xout = runif(10000), rule = 2)$y
quantile(dur_mat_samples, c(0.025, 0.5, 0.975)) # 127 days (87, 205)
hist(dur_mat_samples, breaks = 60, xlab = "dur_mat (days)", main = "")

# save the dur_mat samples 
# saveRDS(dur_mat_samples, paste0("code/fit_durMat/data/dur_mat_samples_fitTO_", tofitSES, "SES.rds"))













# ------------------------------------------------------------------------------
# take a look at the dur_mat when fitting to different SES
# ------------------------------------------------------------------------------
# read in dur_mat_samples generated from fitting to different ses 
dur_mat_samples_bySES <- rbind(data.frame(dur_mat = readRDS("code/fit_durMat/data/dur_mat_samples_fitTO_highSES.rds"), ses = "high"),
      data.frame(dur_mat = readRDS("code/fit_durMat/data/dur_mat_samples_fitTO_mediumSES.rds"), ses = "medium"),
      data.frame(dur_mat = readRDS("code/fit_durMat/data/dur_mat_samples_fitTO_lowSES.rds"), ses = "low"),
      data.frame(dur_mat = readRDS("code/fit_durMat/data/dur_mat_samples_fitTO_unstratifiedSES.rds"), ses = "unstratified")
      )

dur_mat_samples_bySES$ses <- factor(
  dur_mat_samples_bySES$ses,
  levels = c("low", "medium", "high", "unstratified")   
)

ses_summ <- dur_mat_samples_bySES %>%
  group_by(ses) %>%
  summarise(
    med = median(dur_mat),
    lo  = quantile(dur_mat, 0.025),
    hi  = quantile(dur_mat, 0.975),
    .groups = "drop"
  ) %>%
  mutate(lab = sprintf("%.0f (%.0f, %.0f)", med, lo, hi))

y_lab <- max(dur_mat_samples_bySES$dur_mat)

plt.durMat.byses <- 
  ggplot(dur_mat_samples_bySES, aes(x = ses, y = dur_mat, fill = ses)) +
  geom_violin(trim = FALSE, alpha = 0.6, colour = NA) +
  geom_pointrange(
    data = ses_summ,
    aes(x = ses, y = med, ymin = lo, ymax = hi),
    inherit.aes = FALSE, colour = "black", linewidth = 0.5, size = 0.35
  ) +
  geom_text(
    data = ses_summ,
    aes(x = ses, y = y_lab, label = lab),
    inherit.aes = FALSE, size = 3.2, vjust = -0.6
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(
    x = "SES",
    y = "Duration of protection from\nmaternal RSVpreF vaccine (days)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

plt.durMat.byses

# ggsave(plt.durMat.byses, filename = "code/fit_durMat/data/durMat_bySES.pdf", width =6, height = 4)

# remove unstratified ses 
plt.durMat.byses.rmunstra <- 
  ggplot(dur_mat_samples_bySES %>% filter(ses != "unstratified"), aes(x = ses, y = dur_mat, fill = ses)) +
  geom_violin(trim = FALSE, alpha = 0.6, colour = NA) +
  geom_pointrange(
    data = ses_summ %>% filter(ses != "unstratified"),
    aes(x = ses, y = med, ymin = lo, ymax = hi),
    inherit.aes = FALSE, colour = "black", linewidth = 0.5, size = 0.35
  ) +
  geom_text(
    data = ses_summ%>% filter(ses != "unstratified"),
    aes(x = ses, y = y_lab, label = lab),
    inherit.aes = FALSE, size = 3.2, vjust = -0.6
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(
    x = "SES",
    y = "Duration of protection from\nmaternal RSVpreF vaccine (days)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

plt.durMat.byses.rmunstra

# ggsave(plt.durMat.byses.rmunstra, filename = "code/fit_durMat/data/durMat_bySES_rmunstra.pdf", width =6, height = 4)



# ------------------------------------------------------------------------------
# sample from the three SES-specific dur_mat --> create a new dur_mat
# ------------------------------------------------------------------------------
set.seed(1)
dur_mat_pooledSES <- sample((dur_mat_samples_bySES %>% filter(ses != "unstratified"))$dur_mat, 10000, replace = TRUE)
hist(dur_mat_pooledSES, breaks = 50)
summary(dur_mat_pooledSES)
print(quantile(dur_mat_pooledSES, c(0.025, 0.5, 0.975)) )


# plot the by ses and pooled dur_mat -------------
dur_mat_samples_bySES_addPooled <- rbind(dur_mat_samples_bySES,
      data.frame(dur_mat = dur_mat_pooledSES, ses = "pooled")) %>%
  mutate(ses = factor(ses, levels = c("low", "medium", "high", "unstratified", "pooled")))

ses_summ <- dur_mat_samples_bySES_addPooled %>%
  group_by(ses) %>%
  summarise(
    med = median(dur_mat),
    lo  = quantile(dur_mat, 0.025),
    hi  = quantile(dur_mat, 0.975),
    .groups = "drop"
  ) %>%
  mutate(lab = sprintf("%.0f (%.0f, %.0f)", med, lo, hi))

y_lab <- max(dur_mat_samples_bySES_addPooled$dur_mat)

plt.durMat.byses.addPooled <- 
  ggplot(dur_mat_samples_bySES_addPooled %>% filter(ses != "unstratified"), aes(x = ses, y = dur_mat, fill = ses)) +
  geom_violin(trim = FALSE, alpha = 0.6, colour = NA) +
  geom_pointrange(
    data = ses_summ%>% filter(ses != "unstratified"),
    aes(x = ses, y = med, ymin = lo, ymax = hi),
    inherit.aes = FALSE, colour = "black", linewidth = 0.5, size = 0.35
  ) +
  geom_text(
    data = ses_summ%>% filter(ses != "unstratified"),
    aes(x = ses, y = y_lab, label = lab),
    inherit.aes = FALSE, size = 3.2, vjust = -0.6
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(
    x = "SES",
    y = "Duration of protection from\nmaternal RSVpreF vaccine (days)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

plt.durMat.byses.addPooled

# ggsave(plt.durMat.byses.addPooled, filename = "code/fit_durMat/data/durMat_bySES_addPooled.pdf", width =6, height = 4)

# save dur_mat samples (ses = low, medium, high, unstratified, pooled)
# saveRDS(dur_mat_samples_bySES_addPooled, "code/fit_durMat/data/dur_mat_samples_combSES.rds")










# ------------------------------------------------------------------------------
# also get the samples of relative risk of infection for vaxed infants (RR)
# ------------------------------------------------------------------------------
# the RR previously used in the npj paper was from literature
# since we have results from case-control data, use 1- initial VE as the new input RR
RR_median <- 1- ve.cc$VE_median[1] # 0.37 
RR_ub <- 1- ve.cc$VE_lb[1] # 0.73
RR_lb <- 1- ve.cc$VE_ub[1] # 0.15

sigma_L <- (log(RR_median) - log(RR_lb)) / qnorm(0.975)   # 0.4606
sigma_U <- (log(RR_ub) - log(RR_median)) / qnorm(0.975)   # 0.3467

set.seed(1)
n_s  <- 10000
side <- rbinom(n_s, 1, 0.5)
z    <- abs(rnorm(n_s))
RR_samples <- exp(log(RR_median) + ifelse(side == 1, z * sigma_U, -z * sigma_L))

quantile(RR_samples, c(0.025, 0.5, 0.975))   # 0.150 / 0.370 / 0.730
# saveRDS(RR_samples, "code/fit_durMat/data/RR_samples_mat_decrease.rds")



