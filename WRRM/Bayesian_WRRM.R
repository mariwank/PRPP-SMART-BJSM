# install/load libraries 
library(pacman)
p_load(survey, rstan, tidyverse, utils, matrixStats, Rlab) 

### STAN FUNCTIONS/CODE
# Cluster 
source("cs_samp_func.R") 
source("cs_samp_brms_func.R") 
source("grad_par_func.R") 
source("row_subset_func.R") 
source("DEadj_func.R") 

# load stan models
mod_dm <- stan_model('bayes_prpp_smart.stan', auto_write = TRUE) # prpp-smart stan

### STAN SETTINGS

MCMC_SAMPLE <- 5500
BURN.IN <- 500
n_MCMC_chain <- 1
n_MCMC_thin <- 1

## Prior distributions hyperparameter specifications 
alpha_mu00 <- 0 
alpha_mu01 <- 0
alpha_mu10 <- 0
alpha_mu11 <- 0
beta0_mu <- 0
beta1_mu <- 0
theta0_mu <- 0
theta1_mu <- 0
sigma2 <- 1
sigma <- sqrt(sigma2) 
gamma_mu <- 0  


#### Data generation settings
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
# e

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


scenario_vec <- c(1,2,3)

for (z in scenario_vec) {

    # specify what scenario you want to run (1,2,3)
    scenario <- z
    
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
  
    
    ### SIMULATION ###
    # find number of sims needed to get 500 given the above settings 
    source("nsim_toget500_funs.R")
    iterations_needed <- count_iterations()
    n.sim <- iterations_needed
    num_skip <- rep(0, n.sim) # number simulations skipped
    
    ## Store Full WRRM PRPP-SMART results 
    n.DTR <- matrix(NA,nrow = 16, ncol = n.sim) # matrix to store sample size for each DTR path per simulation
    DTR_hat_org = matrix(NA,nrow = 16, ncol = n.sim) # matrix to store preference DTR estimate per simulation from original MCMC sample
    DTR_hat_adj = matrix(NA,nrow = 16, ncol = n.sim) # matrix to store preference DTR estimate per simulation from adjusted MCMC sample
    DTR_var_org = matrix(NA,nrow = 16, ncol = n.sim) # matrix to store preference DTR variance estimate per simulation from original MCMC sample
    DTR_var_adj = matrix(NA,nrow = 16, ncol = n.sim) # matrix to store preference DTR variance estimate per simulation from adjusted MCMC sample
    parameter_hat_org = matrix(NA, nrow=9, ncol=n.sim) # matrix to store parameter estimates per simulation from original MCMC sample
    parameter_hat_adj = matrix(NA, nrow=9, ncol=n.sim) # matrix to store parameter estimates per simulation from adjusted MCMC sample
    parameter_var_org = matrix(NA, nrow=9, ncol=n.sim) # matrix to store parameter variance estimates per simulation from original MCMC sample
    parameter_var_adj = matrix(NA, nrow=9, ncol=n.sim) # matrix to store parameter variance estimates per simulation from adjusted MCMC sample
    ci_hat_org = matrix(NA, nrow=16, ncol=n.sim) # matrix to store whether ci covers truth per simulation from original MCMC sample 
    ci_hat_adj = matrix(NA, nrow=16, ncol=n.sim) # matrix to store whether ci covers truth per simulation from adjusted MCMC sample

    # simulation start
    for (i in 1:n.sim){
      set.seed(i+100000)
      
      # Generate data
      data <- generate_data(method=method, N=N, pNP1=pNP1, pTheta_A=pTheta_A, pNP2=pNP2, pTheta_C=pTheta_C, pi_A=pi_A, pi_B=pi_B, pi_A1=pi_A1, pi_B1=pi_B1, pi_AC=pi_AC, pi_AD=pi_AD, pi_BC=pi_BC, pi_BD=pi_BD, pA0A=pA0A, pB0B=pB0B, pA1A=pA1A, pB1B=pB1B, pA1C1=pA1C1, pA1C0=pA1C0, pA1D1=pA1D1, pA1D0=pA1D0, pA0C1=pA0C1, pA0D1=pA0D1, pB1C1=pB1C1, pB1C0=pB1C0, pB1D1=pB1D1, pB1D0=pB1D0, pB0C1=pB0C1, pB0D1=pB0D1)   
      df <- data[[2]] # data for prpp_smart full analysis replicated and has calculated weights
      
      df$norm_w <- df$w/mean(df$w) # normalize weights to replicated sample size
      
      # check to make sure at least three subjects per treatment path in PRPP-SMART data if not skip simulation
      trialpath_df <- data[[1]] %>% dplyr::group_by(Treatment_Path) %>% dplyr::count() %>% dplyr::filter(n >= 3)
      if (nrow(trialpath_df) < 20) {num_skip[i]=1 ; next} # check to make sure data in each pathway 
      
      ### PRPP-SMART ANALYSIS ###
      
      # create survey design object
      dclus1<-svydesign(id=~id, weights=~norm_w, data=df)
      # create replicated survey design object
      rclus1 <- survey::as.svrepdesign(design = dclus1, type = "JK1") # jacknife1 method 
      
      #Set the Data for Stan
      dfcs <- rclus1$variables
      weights_use <- rclus1$pweights
      data_stan <- with(dfcs, list(y=Y, xalpha00 = xalpha00, xalpha01 = xalpha01, xalpha10=xalpha10, xalpha11=xalpha11, xbeta0 = xbeta0, xbeta1 = xbeta1, xtheta0 = xtheta0, xtheta1 = xtheta1, xgamma=xgamma, N=nrow(dfcs), alpha_mu00=alpha_mu00, alpha_mu01=alpha_mu01, alpha_mu10=alpha_mu10, alpha_mu11=alpha_mu11, beta0_mu=beta0_mu, beta1_mu=beta1_mu, theta0_mu=theta0_mu, theta1_mu=theta1_mu, sigma=sigma, gamma_mu=gamma_mu, weights=weights_use))
      ctrl_stan<-list("chains"=n_MCMC_chain,"iter"=MCMC_SAMPLE,"warmup"=BURN.IN,"thin"=n_MCMC_thin)
      pars_stan <- c("alpha00","alpha01", "alpha10", "alpha11","beta0", "beta1", "theta0", "theta1", "gamma","DTR") # parameters to monitor 
      
      
      mod1 <- cs_sampling(svydes = rclus1, mod_stan = mod_dm,
                          data_stan = data_stan, par_stan = pars_stan, ctrl_stan = ctrl_stan, rep_design = TRUE) # if supplying rep survey object then set rep_design = TRUE if supply srvey design object set rep_design = FALSE
      
      
      # estimates from original sample
      DTR_hat_org[,i] = unname(colMeans(mod1$sampled_parms[,10:25]))
      rownames(DTR_hat_org) <- c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")
      
      parameter_hat_org[,i] = unname(colMeans(mod1$sampled_parms[,1:9]))
      rownames(parameter_hat_org) <- c("alpha00", "alpha01", "alpha10", "alpha11", "beta0", "beta1", "theta0", "theta1", "gamma")
      
      DTR_var_org[,i] = unname(colVars(mod1$sampled_parms[,10:25]))
      rownames(DTR_var_org) <- c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")
      
      parameter_var_org[,i] = unname(colVars(mod1$sampled_parms[,1:9]))
      rownames(parameter_var_org) <- c("alpha00", "alpha01", "alpha10", "alpha11", "beta0", "beta1", "theta0", "theta1", "gamma")
      
      # estimates from adjusted sample
      DTR_hat_adj[,i] = unname(colMeans(mod1$adjusted_parms[,10:25]))
      rownames(DTR_hat_adj) <- c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")
      
      parameter_hat_adj[,i] = unname(colMeans(mod1$adjusted_parms[,1:9]))
      rownames(parameter_hat_adj) <- c("alpha00", "alpha01", "alpha10", "alpha11", "beta0", "beta1", "theta0", "theta1", "gamma")
      
      DTR_var_adj[,i] = unname(colVars(mod1$adjusted_parms[,10:25]))
      rownames(DTR_var_adj) <- c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")
      
      parameter_var_adj[,i] = unname(colVars(mod1$adjusted_parms[,1:9]))
      rownames(parameter_var_adj) <- c("alpha00", "alpha01", "alpha10", "alpha11", "beta0", "beta1", "theta0", "theta1", "gamma")
      
      # check nominal coverage - original sample
      dtr1_quantiles <- c(quantile(mod1$sampled_parms[,10], .025), quantile(mod1$sampled_parms[,10], .975))
      dtr2_quantiles <- c(quantile(mod1$sampled_parms[,11], .025), quantile(mod1$sampled_parms[,11], .975))
      dtr3_quantiles <- c(quantile(mod1$sampled_parms[,12], .025), quantile(mod1$sampled_parms[,12], .975))
      dtr4_quantiles <- c(quantile(mod1$sampled_parms[,13], .025), quantile(mod1$sampled_parms[,13], .975))
      dtr5_quantiles <- c(quantile(mod1$sampled_parms[,14], .025), quantile(mod1$sampled_parms[,14], .975))
      dtr6_quantiles <- c(quantile(mod1$sampled_parms[,15], .025), quantile(mod1$sampled_parms[,15], .975))
      dtr7_quantiles <- c(quantile(mod1$sampled_parms[,16], .025), quantile(mod1$sampled_parms[,16], .975))
      dtr8_quantiles <- c(quantile(mod1$sampled_parms[,17], .025), quantile(mod1$sampled_parms[,17], .975))
      dtr9_quantiles <- c(quantile(mod1$sampled_parms[,18], .025), quantile(mod1$sampled_parms[,18], .975))
      dtr10_quantiles <- c(quantile(mod1$sampled_parms[,19], .025), quantile(mod1$sampled_parms[,19], .975))
      dtr11_quantiles <- c(quantile(mod1$sampled_parms[,20], .025), quantile(mod1$sampled_parms[,20], .975))
      dtr12_quantiles <- c(quantile(mod1$sampled_parms[,21], .025), quantile(mod1$sampled_parms[,21], .975))
      dtr13_quantiles <- c(quantile(mod1$sampled_parms[,22], .025), quantile(mod1$sampled_parms[,22], .975))
      dtr14_quantiles <- c(quantile(mod1$sampled_parms[,23], .025), quantile(mod1$sampled_parms[,23], .975))
      dtr15_quantiles <- c(quantile(mod1$sampled_parms[,24], .025), quantile(mod1$sampled_parms[,24], .975))
      dtr16_quantiles <- c(quantile(mod1$sampled_parms[,25], .025), quantile(mod1$sampled_parms[,25], .975))
      
      quantile_mat <- rbind(dtr1_quantiles, dtr2_quantiles, dtr3_quantiles, dtr4_quantiles, dtr5_quantiles, dtr6_quantiles, dtr7_quantiles, dtr8_quantiles, dtr9_quantiles, dtr10_quantiles, dtr11_quantiles, dtr12_quantiles, dtr13_quantiles, dtr14_quantiles, dtr15_quantiles, dtr16_quantiles)
      param_in_ci <- data.table::between(true_DTR_mat[,1], quantile_mat[, 1], quantile_mat[, 2])
      ci_hat_org[,i] <- param_in_ci
      rownames(ci_hat_org) <- c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")
      
      # check nominal coverage - adjusted sample
      dtr1_quantiles2 <- c(quantile(mod1$adjusted_parms[,10], .025), quantile(mod1$adjusted_parms[,10], .975))
      dtr2_quantiles2 <- c(quantile(mod1$adjusted_parms[,11], .025), quantile(mod1$adjusted_parms[,11], .975))
      dtr3_quantiles2 <- c(quantile(mod1$adjusted_parms[,12], .025), quantile(mod1$adjusted_parms[,12], .975))
      dtr4_quantiles2 <- c(quantile(mod1$adjusted_parms[,13], .025), quantile(mod1$adjusted_parms[,13], .975))
      dtr5_quantiles2 <- c(quantile(mod1$adjusted_parms[,14], .025), quantile(mod1$adjusted_parms[,14], .975))
      dtr6_quantiles2 <- c(quantile(mod1$adjusted_parms[,15], .025), quantile(mod1$adjusted_parms[,15], .975))
      dtr7_quantiles2 <- c(quantile(mod1$adjusted_parms[,16], .025), quantile(mod1$adjusted_parms[,16], .975))
      dtr8_quantiles2 <- c(quantile(mod1$adjusted_parms[,17], .025), quantile(mod1$adjusted_parms[,17], .975))
      dtr9_quantiles2 <- c(quantile(mod1$adjusted_parms[,18], .025), quantile(mod1$adjusted_parms[,18], .975))
      dtr10_quantiles2 <- c(quantile(mod1$adjusted_parms[,19], .025), quantile(mod1$adjusted_parms[,19], .975))
      dtr11_quantiles2 <- c(quantile(mod1$adjusted_parms[,20], .025), quantile(mod1$adjusted_parms[,20], .975))
      dtr12_quantiles2 <- c(quantile(mod1$adjusted_parms[,21], .025), quantile(mod1$adjusted_parms[,21], .975))
      dtr13_quantiles2 <- c(quantile(mod1$adjusted_parms[,22], .025), quantile(mod1$adjusted_parms[,22], .975))
      dtr14_quantiles2 <- c(quantile(mod1$adjusted_parms[,23], .025), quantile(mod1$adjusted_parms[,23], .975))
      dtr15_quantiles2 <- c(quantile(mod1$adjusted_parms[,24], .025), quantile(mod1$adjusted_parms[,24], .975))
      dtr16_quantiles2 <- c(quantile(mod1$adjusted_parms[,25], .025), quantile(mod1$adjusted_parms[,25], .975))
      
      quantile_mat2 <- rbind(dtr1_quantiles2, dtr2_quantiles2, dtr3_quantiles2, dtr4_quantiles2, dtr5_quantiles2, dtr6_quantiles2, dtr7_quantiles2, dtr8_quantiles2, dtr9_quantiles2, dtr10_quantiles2, dtr11_quantiles2, dtr12_quantiles2, dtr13_quantiles2, dtr14_quantiles2, dtr15_quantiles2, dtr16_quantiles2)
      param_in_ci2 <- data.table::between(true_DTR_mat[,1], quantile_mat2[, 1], quantile_mat2[, 2])
      ci_hat_adj[,i] <- param_in_ci2
      rownames(ci_hat_adj) <- c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")
      
    }
    
    
    ######################### EVALUATION ############################
    
    ### Save Settings
    # set date and time for file saving 
    st<-format(Sys.time(), "%Y_%m_%d_%H_%M")
    
    # Define the folder path where you want your results saved
    folder_path <- paste0("/Users/mariwank/Downloads/SimResults/WRRM/Bayesian/Scenario", scenario, "/", subscenario)
    
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
    
    DTR_hat_avg_org = c()
    DTR_hat_avg_adj = c()
    DTR_sd_hat_org = c()
    DTR_sd_hat_adj = c()
    DTR_avg_sd_org = c()
    DTR_avg_sd_adj = c()
    DTR_avg_n = c()
    
    #DTR_avg_n = c()
    for(i in 1:16){
      DTR_hat_avg_org[i] <- round(mean(DTR_hat_org[i,], na.rm = TRUE),4)
      DTR_hat_avg_adj[i] <- round(mean(DTR_hat_adj[i,], na.rm = TRUE),4)
      DTR_sd_hat_org[i] <- round(sd(DTR_hat_org[i,], na.rm = TRUE),4) # sd of mean DTR estimates over 500 sims (sd of 500 posterior means)
      DTR_sd_hat_adj[i] <- round(sd(DTR_hat_adj[i,], na.rm = TRUE),4) # sd of mean DTR estimates over 500 sims (sd of 500 posterior means)
      DTR_avg_sd_org[i] <- round(sqrt(mean(DTR_var_org[i,], na.rm = TRUE)),4) # mean of variance of the DTR estimates over 500 sims (mean of 500 posterior variances then square root to get sd)
      DTR_avg_sd_adj[i] <- round(sqrt(mean(DTR_var_adj[i,], na.rm = TRUE)),4) # mean of variance of the DTR estimates over 500 sims (mean of 500 posterior variances then square root to get sd)
      
      DTR_avg_n[i] <- round(mean(n.DTR[i,], na.rm = TRUE),1)
    }
    
    # calculate bias
    DTR_bias_org = round(DTR_hat_avg_org - expected_pref,4)
    DTR_bias_adj = round(DTR_hat_avg_adj - expected_pref,4)
    
    # calculate rMSE
    rMSE_DTR_org <- sqrt(DTR_sd_hat_org^2 + DTR_bias_org^2)
    rMSE_DTR_adj <- sqrt(DTR_sd_hat_adj^2 + DTR_bias_adj^2)
    
    DTR_results_tbl <- data.frame(DTR=DTR,
                                  True_DTR=expected_pref,
                                  DTR_Hat_Avg_org=DTR_hat_avg_org, 
                                  DTR_Hat_Avg_adj=DTR_hat_avg_adj,
                                  Bias_org=DTR_bias_org, 
                                  Bias_adj=DTR_bias_adj,
                                  Sd_org = DTR_sd_hat_org, 
                                  Sd_adj = DTR_sd_hat_adj, 
                                  Avg_sd_org = DTR_avg_sd_org, 
                                  Avg_sd_adj = DTR_avg_sd_adj, 
                                  rMSE_org = rMSE_DTR_org,
                                  rMSE_adj = rMSE_DTR_adj) 
    
    
    # Define the file name
    file_name <- paste0("DTR_Results_", st, ".csv")
    
    # Create the file path
    file_path <- file.path(folder_path, file_name)
    
    # Write the data to a CSV file
    write.csv(DTR_results_tbl, file = file_path, row.names = FALSE)
    cat("CSV file saved to", file_path, "\n")
    
    
    ## nominal coverage 
    ci_coverage_org = c()
    ci_coverage_adj = c()
    
    for(i in 1:nrow(ci_hat_org)){
      ci_coverage_org[i] <- round(mean(ci_hat_org[i,], na.rm = TRUE),2)
      ci_coverage_adj[i] <- round(mean(ci_hat_adj[i,], na.rm = TRUE),2)
      
    }
    
    ci_results_tbl <- data.frame(DTR, 
                                 CI_Coverage_org=ci_coverage_org, 
                                 CI_Coverage_adj=ci_coverage_adj)
    
    
    
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
    

}
