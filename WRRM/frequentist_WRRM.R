# install/load libraries 
if(!("pacman" %in% installed.packages()[,"Package"])) install.packages("pacman")
library(pacman)
p_load(geepack, geometry, dplyr, MASS, data.table, broom, msm, Rlab, tidyverse, utils) 


### Data generation settings
# Load in data generation function
source("DataGeneration.R")

## USER SETTINGS:
# specify what method you want to generate data for: BJSM or WRRM
method <- "WRRM"


# specify the sub-scenario you want to run
# a
# b
# c
# d

subscenario <- "a"

# Load in data generation parameters depending on sub-scenario you specify

if (subscenario == "a") {
  
  source("SubScenario_a.R")
  
} else if (subscenario == "b") {
  
  source("SubScenario_b.R")
  
} else if (subscenario == "c") {
  
  source("SubScenario_c.R")
  
} else if (subscenario == "d") {
  
  source("SubScenario_d.R")
  
} else if (subscenario == "e") {
  
  source("SubScenario_e.R")
}
# specify what scenario you want to run (1,2,3)
scenario <- 1

if (scenario == 1){
  pNP1=0.50 #  desired proportion of individuals expressing No Preference in stage 1
  pNP2=0.50 # desired proportion of patients expressing No Preference in stage 2 (among non-responders)
  
} else if (scenario == 2){
  pNP1=1/3
  pNP2=1/3
  
}  else if (scenario == 3){
  pNP1=0.5
  pNP2=1/3
  
}

# Specify number of subjects in trial
N <- 1000

# Specify theta targets
pTheta_A=0.4 # desired proportion of individuals expressing preference for treatment A among those with a preference in stage 1
pTheta_C=0.4 # desired proportion of individuals expressing preference for treatment C among those with a preference in stage 2 (among non-responders)


### TRUE DTRS ###

expected_pref <- c()  # expected DTR response rates from our simulated data
expected_pref[1] <- pi_A * pA0A + (1 - pi_A) * pi_AC  #AAC00
expected_pref[2] <- pi_A * pA0A + (1 - pi_A) * pi_AD  #AAD00
expected_pref[3] <- pi_B * pB0B + (1 - pi_B) * pi_BC  #BBC00
expected_pref[4] <- pi_B * pB0B + (1 - pi_B) * pi_BD  #BBD00
expected_pref[5] <- pi_A * pA0A + (1 - pi_A) * pA0C1 #AAC01
expected_pref[6] <- pi_A * pA0A + (1 - pi_A) * pA0D1 #AAD01
expected_pref[7] <- pi_B * pB0B + (1 - pi_B) * pB0C1 #BBC01
expected_pref[8] <- pi_B * pB0B + (1 - pi_B) * pB0D1 #BBD01
expected_pref[9] <- pi_A1 * pA1A + (1 - pi_A1) * pA1C0 #AAC10
expected_pref[10] <- pi_A1 * pA1A + (1 - pi_A1) * pA1D0 #AAD10
expected_pref[11] <- pi_B1 * pB1B + (1 - pi_B1) * pB1C0 #BBC10
expected_pref[12] <- pi_B1 * pB1B + (1 - pi_B1) * pB1D0 #BBD10
expected_pref[13] <- pi_A1 * pA1A + (1 - pi_A1) * pA1C1 #AAC11
expected_pref[14] <- pi_A1 * pA1A + (1 - pi_A1) * pA1D1 #AAD11
expected_pref[15] <- pi_B1 * pB1B + (1 - pi_B1) * pB1C1 #BBC11
expected_pref[16] <- pi_B1 * pB1B + (1 - pi_B1) * pB1D1 #BBD11

# Define contrast matrix for DTR calculation
contrast_dtr <- matrix(c(1,0,0,0,1,0,1,0,1, #AC00
                         1,0,0,0,1,0,-1,0,-1, #AD00
                         1,0,0,0,-1,0,1,0,-1, #BC00
                         1,0,0,0,-1,0,-1,0,1, #BD00
                         0,1,0,0,1,0,0,1,1, #AC01
                         0,1,0,0,1,0,0,-1,-1, #AD01
                         0,1,0,0,-1,0,0,1,-1, #BC01
                         0,1,0,0,-1,0,0,-1,1, #BD01
                         0,0,1,0,0,1,1,0,1, #AC10
                         0,0,1,0,0,1,-1,0,-1, #AD10
                         0,0,1,0,0,-1,1,0,-1, #BC10
                         0,0,1,0,0,-1,-1,0,1, #BD10
                         0,0,0,1,0,1,0,1,1, #AC11
                         0,0,0,1,0,1,0,-1,-1, #AD11
                         0,0,0,1,0,-1,0,1,-1, #BC11
                         0,0,0,1,0,-1,0,-1,1 #BD11
),nrow = 16, ncol=9, byrow = TRUE)

true_DTR_mat <- matrix(c(expected_pref), ncol = 1)
rownames(true_DTR_mat) <- c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")


### SIMULATION ###
# find number of sims needed to get 500 given the above settings 
source("nsim_toget500_funs.R")
iterations_needed <- count_iterations()
n.sim <- iterations_needed
num_skip <- rep(0, n.sim) # number simulations skipped

## Store WRRM PRPP-SMART results 
DTR_hat <- matrix(NA,nrow = 16, ncol = n.sim) # matrix to store preference DTR estimate per simulation 
parameter_hat <- matrix(NA, nrow=9, ncol=n.sim) # matrix to store parameter estimates per simulation 
variance_param_hat <- matrix(NA, nrow=9, ncol=n.sim) # matrix to store sandwich variance of parameter estimates per simulation
variance_dtr_hat <- matrix(NA, nrow=16, ncol=n.sim) # matrix to store delta method DTR variance estimates per simulation
ci_hat <- matrix(NA, nrow=16, ncol=n.sim) # matrix to store whether ci covers truth per simulation  
n.DTR <- matrix(NA,nrow = 16, ncol = n.sim) # matrix to store sample size for each DTR path per simulation

# simulation start
for (i in 1:n.sim){
  set.seed(i+100000)
  
  # Generate data
  data <- generate_data(method = method, N=N, pNP1=pNP1, pTheta_A=pTheta_A, pNP2=pNP2, pTheta_C=pTheta_C, pi_A=pi_A, pi_B=pi_B, pi_A1=pi_A1, pi_B1=pi_B1, pi_AC=pi_AC, pi_AD=pi_AD, pi_BC=pi_BC, pi_BD=pi_BD, pA0A=pA0A, pB0B=pB0B, pA1A=pA1A, pB1B=pB1B, pA1C1=pA1C1, pA1C0=pA1C0, pA1D1=pA1D1, pA1D0=pA1D0, pA0C1=pA0C1, pA0D1=pA0D1, pB1C1=pB1C1, pB1C0=pB1C0, pB1D1=pB1D1, pB1D0=pB1D0, pB0C1=pB0C1, pB0D1=pB0D1)   
  df <- data[[2]] # data for prpp_smart full analysis replicated and has calculated weights
  
  df <-df[order(df$id),] # sort data by id for gee statement
  
  # check to make sure at least three subjects per treatment path in PRPP-SMART data if not skip simulation
  trialpath_df <- data[[1]] %>% dplyr::group_by(Treatment_Path) %>% dplyr::count() %>% dplyr::filter(n >= 3)
  if (nrow(trialpath_df) < 20) {num_skip[i]=1 ; next} # check to make sure data in each pathway
  
  ### WRRM PRPP-SMART ANALYSIS ###
  
  # fit 
  gee.fit <- geeglm(Y ~ xalpha00 + xalpha01 + xalpha10 + xalpha11 + xbeta0 + xbeta1 + xtheta0 + xtheta1 + xgamma-1, id=id, weights=w, family = "binomial", corstr = "independence", data = df) # warning is normal since using non-integer weights
  
  ## parameter estimates
  parameter_hat[,i] <- gee.fit$coefficients
  rownames(parameter_hat) <- c("alpha00", "alpha01", "alpha10", "alpha11", "beta0", "beta1","theta0", "theta1","gamma")
  
  ## DTR estimates 
  for(j in 1:16){
    DTR_hat[j,i] = exp((contrast_dtr[j,]%*%parameter_hat[,i]))/(1 + exp((contrast_dtr[j,]%*%parameter_hat[,i])))
    
  }
  rownames(DTR_hat) <- c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")
  
  
  ## Parameter Variances 
  
  # robust variances of each parameter 
  variance_param_hat[,i] <- diag(gee.fit$geese$vbeta)
  rownames(variance_param_hat) <- c("alpha00", "alpha01", "alpha10", "alpha11", "beta0", "beta1","theta0", "theta1","gamma")
  
  ## Delta method for DTR variances
  
  a00 <- parameter_hat[1,i]
  a01 <- parameter_hat[2,i]
  a10 <- parameter_hat[3,i]
  a11 <- parameter_hat[4,i]
  b0 <- parameter_hat[5,i]
  b1 <- parameter_hat[6,i]
  t0 <- parameter_hat[7,i]
  t1 <- parameter_hat[8,i]
  g <- parameter_hat[9,i]
  
  # Extract variance covariance matrix from geeglm
  varcov <- gee.fit$geese$vbeta
  
  # delta method for dtr variances
  for (d in 1:16){
    
    if (d == 1){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 + x2 + x3 + x4))), c(a00, b0, t0, g), varcov[c(1,5,7,9), c(1,5,7,9)]))^2
    }
    
    if (d == 2){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 + x2 - x3 - x4))), c(a00, b0, t0, g), varcov[c(1,5,7,9), c(1,5,7,9)]))^2
      
    }
    
    if (d == 3){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 - x2 + x3 - x4))), c(a00, b0, t0, g), varcov[c(1,5,7,9), c(1,5,7,9)]))^2
      
    }
    
    if (d == 4){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 - x2 - x3 + x4))), c(a00, b0, t0, g), varcov[c(1,5,7,9), c(1,5,7,9)]))^2
      
    }
    
    if (d == 5){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 + x2 + x3 + x4))), c(a01, b0, t1, g), varcov[c(2,5,8,9), c(2,5,8,9)]))^2
      
    }
    
    if (d == 6){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 + x2 - x3 - x4))), c(a01, b0, t1, g), varcov[c(2,5,8,9), c(2,5,8,9)]))^2
      
    }
    
    if (d == 7){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 - x2 + x3 - x4))), c(a01, b0, t1, g), varcov[c(2,5,8,9), c(2,5,8,9)]))^2
      
    }
    
    if (d == 8){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 - x2 - x3 + x4))), c(a01, b0, t1, g), varcov[c(2,5,8,9), c(2,5,8,9)]))^2
      
    }
    
    if (d == 9){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 + x2 + x3 + x4))), c(a10, b1, t0, g), varcov[c(3,6,7,9), c(3,6,7,9)]))^2
      
    }
    
    if (d == 10){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 + x2 - x3 - x4))), c(a10, b1, t0, g), varcov[c(3,6,7,9), c(3,6,7,9)]))^2
      
    }
    
    if (d == 11){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 - x2 + x3 - x4))), c(a10, b1, t0, g), varcov[c(3,6,7,9), c(3,6,7,9)]))^2
      
    }
    
    if (d == 12){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 - x2 - x3 + x4))), c(a10, b1, t0, g), varcov[c(3,6,7,9), c(3,6,7,9)]))^2
      
    }
    
    if (d == 13){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 + x2 + x3 + x4))), c(a11, b1, t1, g), varcov[c(4,6,8,9), c(4,6,8,9)]))^2
      
    }
    
    if (d == 14){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 + x2 - x3 - x4))), c(a11, b1, t1, g), varcov[c(4,6,8,9), c(4,6,8,9)]))^2
      
    }
    
    if (d == 15){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 - x2 + x3 - x4))), c(a11, b1, t1, g), varcov[c(4,6,8,9), c(4,6,8,9)]))^2
      
    }
    
    if (d == 16){
      variance_dtr_hat[d,i] <- (deltamethod(~ 1/(1+exp(-(x1 - x2 - x3 + x4))), c(a11, b1, t1, g), varcov[c(4,6,8,9), c(4,6,8,9)]))^2
      
    }
  }
  
  rownames(variance_dtr_hat) <- c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")
  
  
  ## calculate CI coverage of DTRs
  
  lwr_dtr <- c()
  upper_dtr <- c()
  
  for(d in 1:16){
    lwr_dtr[d] <- DTR_hat[d,i] - qnorm((1+0.95)/2)*sqrt(variance_dtr_hat[d,i])
    upper_dtr[d] <- DTR_hat[d,i] + qnorm((1+0.95)/2)*sqrt(variance_dtr_hat[d,i])
    
  }
  dtr_ci_mat <- cbind(lwr_dtr, upper_dtr)
  rownames(dtr_ci_mat) <-  c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")
  
  param_in_ci <- data.table::between(true_DTR_mat[,1], dtr_ci_mat[, 1], dtr_ci_mat[, 2])
  ci_hat[,i] <- param_in_ci
  rownames(ci_hat) <- c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")
  
}

######################### EVALUATION ############################

### Save Settings
# set date and time for file saving 
st<-format(Sys.time(), "%Y_%m_%d_%H_%M")

# Define the folder path where you want your results saved
#folder_path <- "path to your folder where you want to store frequentist results"

# Example: 
folder_path <- paste0("/Users/mariwank/Downloads/SimResults/WRRM/Scenario", scenario, "/", subscenario)


# Create the folder if it doesn't exist
if (!file.exists(folder_path)) {
  dir.create(folder_path, recursive = TRUE)
  cat("Folder created at", folder_path, "\n")
} else {
  cat("Folder already exists at", folder_path, "\n")
}

### DTR Full Model Results 
DTR = c("AAC00","AAD00","BBC00","BBD00","AAC01","AAD01","BBC01","BBD01",
        "AAC10","AAD10","BBC10","BBD10","AAC11","AAD11","BBC11","BBD11")
DTR_hat_avg = c()
DTR_sd_hat = c()
DTR_avg_sd = c()
DTR_avg_n = c()

for(i in 1:16){
  DTR_hat_avg[i] <- round(mean(DTR_hat[i,], na.rm = TRUE),4)
  DTR_sd_hat[i] <- round(sd(DTR_hat[i,], na.rm = TRUE),4) # sd of mean DTR estimates over 500 sims (sd of 500 DTR estimates)
  DTR_avg_sd[i] <- round(sqrt(mean(variance_dtr_hat[i,], na.rm = TRUE)),4) # mean of variance of the DTR estimates over 500 sims (mean of 500 variances then square root to get sd) sqrt(mean(var(dtr hat))) equivalent to mean(sqrt(var(dtr hat)))
  DTR_avg_n[i] <- round(mean(n.DTR[i,], na.rm = TRUE),1)
  
}

# calculate bias
DTR_bias = round(DTR_hat_avg - expected_pref,4)

# calculate rMSE
rMSE_DTR <- sqrt(DTR_sd_hat^2 + DTR_bias^2)


DTR_results_tbl <- data.frame(DTR=DTR,True_DTR=expected_pref,
                              DTR_Hat_Avg=DTR_hat_avg, 
                              Bias=DTR_bias, 
                              SE = DTR_sd_hat, 
                              Avg_se = DTR_avg_sd,
                              rMSE = rMSE_DTR) 

# Define the file name
file_name <- paste0("DTR_Results_", st, ".csv")

# Create the file path
file_path <- file.path(folder_path, file_name)

# Write the data to a CSV file
write.csv(DTR_results_tbl, file = file_path, row.names = FALSE)
cat("CSV file saved to", file_path, "\n")

## nominal coverage 
ci_coverage = c()

for(i in 1:nrow(ci_hat)){
  ci_coverage[i] <- mean(ci_hat[i,], na.rm = TRUE)
  
}

ci_results_tbl <- data.frame(DTR, 
                             CI_Coverage=ci_coverage)

# Define the file name
file_name <- paste0("DTR_CI_Results_", st, ".csv")

# Create the file path
file_path <- file.path(folder_path, file_name)

# Write the data to a CSV file
write.csv(ci_results_tbl, file = file_path, row.names = FALSE)
cat("CSV file saved to", file_path, "\n")


## number simulations skipped due to positivity assumption
num_skip_total <- sum(num_skip)

mean_skip <- mean(num_skip)
skip_df <- data.frame(num_sim = n.sim, 
                      num_skip_total = num_skip_total, 
                      avg_num_skip = mean_skip) 

# Define the file name
file_name <- paste0("Num_Sim_Skip_Summary_", st, ".csv")

# Create the file path
file_path <- file.path(folder_path, file_name)

# Write the data to a CSV file
write.csv(skip_df, file = file_path, row.names = FALSE)
cat("CSV file saved to", file_path, "\n")



