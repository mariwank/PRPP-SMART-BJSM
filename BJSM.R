# install/load libraries 
if(!("pacman" %in% installed.packages()[,"Package"])) install.packages("pacman")
library(pacman)
p_load(R2jags, coda, tidyverse, utils) 


### Data generation settings
# Load in data generation function
source("DataGeneration.R")

## USER SETTINGS:
# specify what method you want to generate data for: BJSM or WRRM
method <- "BJSM"


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
  
} 

# specify the prior setting (i.e., hyperparameter values) for BJSM 
# 1
# 2
# 3
# 4

prior_setting <- 1

if (prior_setting == 1) {
  
  pi_prior.a <- c(rep(1/3,6)) # pi_J, pi_jk ~ Beta(a,b)
  pi_prior.b <- c(rep(1/3,6)) # pi_J, pi_jk ~ Beta(a,b)
  beta1_prior.r <- 2          # Beta1j ~ Gamma(r, mu) 
  beta1_prior.mu <- 2         # Beta1j ~ Gamma(r, mu)
  alphaP_prior.r <- 2         
  alphaP_prior.mu <- 2
  
  
} else if (prior_setting == 2) {
  
  pi_prior.a <- c(rep(1/3,6))  # pi_J, pi_jk ~ Beta(a,b)
  pi_prior.b <- c(rep(1/3,6))  # pi_J, pi_jk ~ Beta(a,b)
  beta1_prior.r <- 1           # Beta1j ~ Gamma(r, mu)
  beta1_prior.mu <- 1          # Beta1j ~ Gamma(r, mu)
  alphaP_prior.r <- 1          # Alpha1j = Alpha2k ~ Gamma(r, mu) 
  alphaP_prior.mu <- 1
  
} else if (prior_setting == 3) {
  
  pi_prior.a <- c(rep(1/3,6))  # pi_J, pi_jk ~ Beta(a,b)
  pi_prior.b <- c(rep(1/3,6))  # pi_J, pi_jk ~ Beta(a,b)
  beta1_prior.r <- 1/2         # Beta1j ~ Gamma(r, mu)    
  beta1_prior.mu <- 1/2        # Beta1j ~ Gamma(r, mu)
  alphaP_prior.r <- 1/2        # Alpha1j = Alpha2k ~ Gamma(r, mu) 
  alphaP_prior.mu <- 1/2       # Alpha1j = Alpha2k ~ Gamma(r, mu) 
  
} else if (prior_setting == 4) {
  
  pi_prior.a <- c(rep(1,6))  # pi_J, pi_jk ~ Beta(a,b)
  pi_prior.b <- c(rep(1,6))  # pi_J, pi_jk ~ Beta(a,b)
  beta1_prior.r <- 1/2       # Beta1j ~ Gamma(r, mu)  
  beta1_prior.mu <- 1/2      # Beta1j ~ Gamma(r, mu)
  alphaP_prior.r <- 1/2      # Alpha1j = Alpha2k ~ Gamma(r, mu) 
  alphaP_prior.mu <- 1/2     # Alpha1j = Alpha2k ~ Gamma(r, mu) 
  
} 


# Specify number of subjects in trial
N <- 1000

# Specify theta targets
pTheta_A=0.4 # desired proportion of individuals expressing preference for treatment A among those with a preference in stage 1
pTheta_C=0.4 # desired proportion of individuals expressing preference for treatment C among those with a preference in stage 2 (among non-responders)


### CALCULATE TRUE DTRS ###
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


# Create to check nominal coverage
true_piDTR_mat <- matrix(c(pi_A, pi_B, pi_AC, pi_AD, pi_BC, pi_BD, expected_pref, alphaP_link, beta1_link))
rownames(true_piDTR_mat) <- c("PiA", "PiB", "PiAC", "PiAD", "PiBC", "PiBD", "AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11", "alphaP1A", "alphaP1B", "alphaP2C", "alphaP2D", "Beta1A", "Beta1B")


scenario_vec <- c(1,2,3) # run each sub-scenario across three primary preference scenarios

for (s in scenario_vec) {

   
    scenario <- s
    
    
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
    
    # find number of sims needed to get 500 given the above settings 
    source("nsim_toget500_funs.R")
    iterations_needed <- count_iterations()
    
    # BJSM Set-Up
    
    ##### Estimation #####
    NUM_PATHS <- 6 # number of treatment paths independent of preference in SMART (AA, BB, AC, AD, BC, BD) OR max(treatment_stageII)
    MCMC_SAMPLE <- 5500
    BURN.IN <- 500
    n_MCMC_chain <- 1
    n_MCMC_thin <- 1
    saved_parms = c('pi','beta', 'alphaP','DTR') # parameters to monitor 
    

    # Jags model
    bjsm_model = function()
    { 
      for (i in 1:N){   # N is total sample size
        # likelihood
        Y1[i]~dbern(pi_1[i])
        Y2[i]~dbern(pi_2[i])
        # explaining
        pi_1[i] <- ifelse(preference_stageI[i] == 0, pi[treatment_stageI[i]], # A0/B0
                          ifelse(preference_stageI[i] == 1 && treatment_stageI[i] == 1, alphaP[1] * pi[treatment_stageI[i]], alphaP[2] * pi[treatment_stageI[i]])) # A1/B1
        
        
        pi_2[i] <- ifelse(Y1[i] == 1 && preference_stageI[i] == 0, 
                          pi[treatment_stageII[i]] * beta[treatment_stageI[i]], # A0A/B0B
                          ifelse(Y1[i] == 1 && preference_stageI[i] == 1 && treatment_stageI[i] == 1,  alphaP[1] * pi[treatment_stageII[i]] * beta[treatment_stageI[i]], # A1A
                                 ifelse(Y1[i] == 1 && preference_stageI[i] == 1 && treatment_stageI[i] == 2, alphaP[2] * pi[treatment_stageII[i]] * beta[treatment_stageI[i]], # B1B
                                        ifelse(Y1[i] == 0 &&  preference_stageI[i] == 0 && preference_stageII[i] == 0, pi[treatment_stageII[i]], # A0C0/A0D0/B0C0/B0D0
                                               ifelse(Y1[i] == 0 &&  preference_stageI[i] == 0 && preference_stageII[i] == 1 && treatmentCD_stageII[i] == 3, alphaP[3] * pi[treatment_stageII[i]], # A0C1/B0C1
                                                      ifelse(Y1[i] == 0 &&  preference_stageI[i] == 0 && preference_stageII[i] == 1 && treatmentCD_stageII[i] == 4, alphaP[4] * pi[treatment_stageII[i]], # A0D1/B0D1
                                                             ifelse(Y1[i] == 0 &&  preference_stageI[i] == 1 && preference_stageII[i] == 0 && treatment_stageI[i] == 1, alphaP[1] * pi[treatment_stageII[i]], # A1C0/A1D0
                                                                    ifelse(Y1[i] == 0 &&  preference_stageI[i] == 1 && preference_stageII[i] == 0 && treatment_stageI[i] == 2, alphaP[2] * pi[treatment_stageII[i]], # B1C0/B1D0 
                                                                           ifelse(Y1[i] == 0 &&  preference_stageI[i] == 1 && preference_stageII[i] == 1 && treatment_stageI[i] == 1 && treatmentCD_stageII[i] == 3, alphaP[1] * alphaP[3] * pi[treatment_stageII[i]], # A1C1
                                                                                ifelse(Y1[i] == 0 &&  preference_stageI[i] == 1 && preference_stageII[i] == 1 && treatment_stageI[i] == 1 && treatmentCD_stageII[i] == 4, alphaP[1] * alphaP[4] * pi[treatment_stageII[i]], # A1D1
                                                                                    ifelse(Y1[i] == 0 &&  preference_stageI[i] == 1 && preference_stageII[i] == 1 && treatment_stageI[i] == 2 && treatmentCD_stageII[i] == 3, alphaP[2] * alphaP[3] * pi[treatment_stageII[i]], alphaP[2] * alphaP[4] * pi[treatment_stageII[i]]))))))))))) # B1C1 and B1D1
        
        
        
        
      }
      
      alphaP[1] ~ dgamma(alphaP_prior_a,alphaP_prior_b) # alpha1A
      alphaP[2] ~ dgamma(alphaP_prior_a,alphaP_prior_b) # alpha1B
      alphaP[3] ~ dgamma(alphaP_prior_a,alphaP_prior_b) # alpha2C
      alphaP[4] ~ dgamma(alphaP_prior_a,alphaP_prior_b) # alpha2D
      beta[1] ~ dgamma(beta1_prior_a,beta1_prior_b)     # beta1A
      beta[2] ~ dgamma(beta1_prior_a,beta1_prior_b)     # beta1B
      
      for (j in 1:num_paths){
        pi[j] ~ dbeta(pi_prior_a[j],pi_prior_b[j]); T(0.001,0.999)
      }
      
      # DTR
      DTR[1] = pi[1] * (pi[1]*beta[1]) + (1 - pi[1]) * pi[3] ## AAC00
      DTR[2] = pi[1] * (pi[1]*beta[1]) + (1 - pi[1]) * pi[4] ## AAD00
      DTR[3] = pi[2] * (pi[2]*beta[2]) + (1 - pi[2]) * pi[5] ## BBC00
      DTR[4] = pi[2] * (pi[2]*beta[2]) + (1 - pi[2]) * pi[6] ## BBD00
      DTR[5] = pi[1] * (pi[1]*beta[1]) + (1 - pi[1]) * (alphaP[3]*pi[3]) ## AAC01
      DTR[6] = pi[1] * (pi[1]*beta[1]) + (1 - pi[1]) * (alphaP[4]*pi[4]) ## AAD01
      DTR[7] = pi[2] * (pi[2]*beta[2]) + (1 - pi[2]) * (alphaP[3]*pi[5]) ## BBC01
      DTR[8] = pi[2] * (pi[2]*beta[2]) + (1 - pi[2]) * (alphaP[4]*pi[6])  ## BBD01
      DTR[9] = (alphaP[1]*pi[1]) * (pi[1]*alphaP[1]*beta[1]) + (1 - alphaP[1]*pi[1]) * (alphaP[1]*pi[3]) ## AAC10
      DTR[10] = (alphaP[1]*pi[1]) * (pi[1]*alphaP[1]*beta[1]) + (1 - alphaP[1]*pi[1]) * (alphaP[1]*pi[4]) ## AAD10
      DTR[11] = (alphaP[2]*pi[2]) * (pi[2]*alphaP[2]*beta[2]) + (1 - alphaP[2]*pi[2]) * (alphaP[2]*pi[5]) ## BBC10
      DTR[12] = (alphaP[2]*pi[2]) * (pi[2]*alphaP[2]*beta[2]) + (1 - alphaP[2]*pi[2]) * (alphaP[2]*pi[6]) ## BBD10
      DTR[13] = (alphaP[1]*pi[1]) * (pi[1]*alphaP[1]*beta[1])  + (1 - alphaP[1]*pi[1]) * (alphaP[1]*alphaP[3]*pi[3]) ## AAC11
      DTR[14] = (alphaP[1]*pi[1]) * (pi[1]*alphaP[1]*beta[1])  + (1 - alphaP[1]*pi[1]) * (alphaP[1]*alphaP[4]*pi[4]) ## AAD11
      DTR[15] =  (alphaP[2]*pi[2]) * (pi[2]*alphaP[2]*beta[2])  + (1 - alphaP[2]*pi[2]) * (alphaP[2]*alphaP[3]*pi[5]) ## BBC11
      DTR[16] = (alphaP[2]*pi[2]) * (pi[2]*alphaP[2]*beta[2])  + (1 - alphaP[2]*pi[2]) * (alphaP[2]*alphaP[4]*pi[6]) ## BBD11
      
    }
    
    ### SIMULATION ###
    n.sim <- iterations_needed
    num_skip <- rep(0, n.sim) # number simulations skipped
    
    # store results
    DTR_hat_bjsm <- c()            # store posterior means 
    DTR_hat_bjsm_var <- c()        # store posterior variances
    pi_hat_bjsm <- c()             # store posterior means
    pi_hat_bjsm_var <- c()         # store posterior variances
    beta_hat_bjsm <- c()           # store posterior means
    beta_hat_bjsm_var <- c()       # store posterior variances
    alpha_hat_bjsm <- c()          # store posterior means
    alpha_hat_bjsm_var <- c()      # store posterior variances
    ci_hat <- c()
    
    for (i in 1:n.sim) {
      set.seed(i+100000)
      data <- generate_data(method=method, N=N, pNP1=pNP1, pTheta_A=pTheta_A, pNP2=pNP2, pTheta_C=pTheta_C, pi_A=pi_A, pi_B=pi_B, pi_A1=pi_A1, pi_B1=pi_B1, pi_AC=pi_AC, pi_AD=pi_AD, pi_BC=pi_BC, pi_BD=pi_BD, pA0A=pA0A, pB0B=pB0B, pA1A=pA1A, pB1B=pB1B, pA1C1=pA1C1, pA1C0=pA1C0, pA1D1=pA1D1, pA1D0=pA1D0, pA0C1=pA0C1, pA0D1=pA0D1, pB1C1=pB1C1, pB1C0=pB1C0, pB1D1=pB1D1, pB1D0=pB1D0, pB0C1=pB0C1, pB0D1=pB0D1)   
      
      bjsm_df <- data[[2]] # BJSM formatted data
      
      # check to make sure simulated data has at least 3 subjects per treatment path
      trialpath_df <- data[[1]] %>% dplyr::group_by(Treatment_Path) %>% dplyr::count() %>% dplyr::filter(n >= 3)
      if (nrow(trialpath_df) < 20) {num_skip[i]=1 ; next} # check to make sure data in each pathway 
      
      data_list <- list(N = nrow(bjsm_df),
                        num_paths = NUM_PATHS, 
                        Y1 = bjsm_df$response_stageI,
                        Y2 = bjsm_df$response_stageII,
                        treatment_stageI = bjsm_df$treatment_stageI,
                        treatment_stageII = bjsm_df$treatment_stageII,
                        preference_stageI = bjsm_df$preference_stageI,
                        preference_stageII = bjsm_df$preference_stageII,
                        treatmentCD_stageII = bjsm_df$treatmentCD_stageII,
                        pi_prior_a = pi_prior.a, 
                        pi_prior_b = pi_prior.b, 
                        alphaP_prior_a = alphaP_prior.r, 
                        alphaP_prior_b = alphaP_prior.mu, 
                        beta1_prior_a = beta1_prior.r,   
                        beta1_prior_b = beta1_prior.mu  
      )
      
      jags.out <- jags(data=data_list, model.file=bjsm_model,parameters.to.save=saved_parms,n.chains=n_MCMC_chain,
                       n.thin=n_MCMC_thin, n.iter=MCMC_SAMPLE,n.burnin=BURN.IN, progress.bar = "none", quiet=TRUE)
      
      
      store <- print(jags.out)
      mcmc.out <- summary(coda::as.mcmc(jags.out))
      
      
      # pi estimates posterior mean
      pi_hat_bjsm <- rbind(pi_hat_bjsm, store$mean$pi) # 6 columns n.sim rows "PiA", "PiB", "PiAC", "PiAD", "PiBC", "PiBD"
      
      # pi estimates posterior variance
      pi_hat_bjsm_var <- rbind(pi_hat_bjsm_var, store$sd$pi^2) # 6 columns n.sim rows "PiA", "PiB", "PiAC", "PiAD", "PiBC", "PiBD"
      
      # DTR estimates posterior mean
      DTR_hat_bjsm <- rbind(DTR_hat_bjsm, store$mean$DTR) # 16 columsn n.sim rows "AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11"
      
      # DTR estimates posterior variance
      DTR_hat_bjsm_var <- rbind(DTR_hat_bjsm_var, store$sd$DTR^2) # 16 columsn n.sim rows "AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11"
      
      # beta estimates posterior mean
      beta_hat_bjsm <- rbind(beta_hat_bjsm, store$mean$beta) # 2 columsn n.sim rows "Beta1A", "Beta1B"
      
      # beta estimates posterior variances
      beta_hat_bjsm_var <- rbind(beta_hat_bjsm_var, store$sd$beta^2) # 2 columsn n.sim rows "Beta1A", "Beta1B"
      
      # alpha estimates posterior mean
      alpha_hat_bjsm <- rbind(alpha_hat_bjsm, store$mean$alphaP) # 4 columsn n.sim rows "Alpha1A", "Alpha1B", "Alpha2C", "Alpha2D"
      
      # alpha estimates posterior variances
      alpha_hat_bjsm_var <- rbind(alpha_hat_bjsm_var, store$sd$alphaP^2) # 4 columsn n.sim rows "Alpha1A", "Alpha1B", "Alpha2C", "Alpha2D"
      
      # Coverage 
      quantile_mat <- mcmc.out$quantiles[c(1:6,8:29),c(1,5)] 
      # reorder quantile mat to be in the same order as true_piDTR_mat
      rowname_order <- c("pi[1]", "pi[2]", "pi[3]", "pi[4]", "pi[5]", "pi[6]", "DTR[1]", "DTR[2]", "DTR[3]", "DTR[4]", "DTR[5]", "DTR[6]", "DTR[7]", "DTR[8]", "DTR[9]", "DTR[10]", "DTR[11]", "DTR[12]", "DTR[13]", "DTR[14]", "DTR[15]", "DTR[16]", "alphaP[1]", "alphaP[2]", "alphaP[3]", "alphaP[4]", "beta[1]", "beta[2]")
      quantile_mat <- quantile_mat[rowname_order,,drop=FALSE]
      param_in_ci <- data.table::between(true_piDTR_mat[,1], quantile_mat[, 1], quantile_mat[, 2])
      ci_hat <- rbind(ci_hat, param_in_ci) 
      
    }
    
    ######################### EVALUATION ############################
    
    ### Save Settings
    # set date and time for file saving 
    st<-format(Sys.time(), "%Y_%m_%d_%H_%M")
    
    # Define the folder path where you want your results saved
    # folder_path <- "path to your folder where you want to store results"
    
    # Example: 
    folder_path <- paste0("/Users/mariwank/Downloads/SimResults/Scenario", scenario, "/", subscenario, "/PriorSetting", prior_setting)
    
    # Create the folder if it doesn't exist
    if (!file.exists(folder_path)) {
      dir.create(folder_path, recursive = TRUE)
      cat("Folder created at", folder_path, "\n")
    } else {
      cat("Folder already exists at", folder_path, "\n")
    }
    
    ## Linkage Parameters (alpha/beta)
    beta_alpha_effect_output <- data.frame(link_param = c("Beta1A", "Beta1B", "AlphaP1A", "AlphaP1B", "AlphaP2C", "AlphaP2D"),
                                           true_param = c(beta1_link[1], beta1_link[2], alphaP_link[1], alphaP_link[2], alphaP_link[3], alphaP_link[4]),
                                           param_hat = c(apply(beta_hat_bjsm,2,mean, na.rm = TRUE), apply(alpha_hat_bjsm,2,mean, na.rm = TRUE)),
                                           sd_param_hat = c(apply(beta_hat_bjsm,2,sd, na.rm = TRUE), apply(alpha_hat_bjsm,2,sd, na.rm = TRUE)),
                                           avg_sd_param_hat = c(sqrt(apply(beta_hat_bjsm_var,2,mean, na.rm = TRUE)), sqrt(apply(alpha_hat_bjsm_var,2,mean, na.rm = TRUE))))
    beta_alpha_effect_output$bias <- beta_alpha_effect_output$param_hat - beta_alpha_effect_output$true_param
    beta_alpha_effect_output$rMSE <- sqrt(beta_alpha_effect_output$sd_param_hat^2 + beta_alpha_effect_output$bias^2)
    
    # Define the file name
    file_name <- paste0("LinkageParam_results_", st, ".csv")
    
    # Create the file path
    file_path <- file.path(folder_path, file_name)
    
    # Write the data to a CSV file
    write.csv(beta_alpha_effect_output, file = file_path, row.names = FALSE)
    cat("CSV file saved to", file_path, "\n")
    
    ## Indifference Main Effects
    trt_effect_output <- data.frame(pi = c("PiA", "PiB", "PiAC", "PiAD", "PiBC", "PiBD"),
                                    true_pi = c(pi_A,pi_B,pi_AC,pi_AD,pi_BC,pi_BD),
                                    pi_hat = apply(pi_hat_bjsm,2,mean, na.rm = TRUE),
                                    sd_pi_hat = apply(pi_hat_bjsm,2,sd, na.rm = TRUE),
                                    avg_sd_pi_hat = sqrt(apply(pi_hat_bjsm_var,2,mean, na.rm = TRUE)))
    trt_effect_output$bias <- trt_effect_output$pi_hat - trt_effect_output$true_pi
    trt_effect_output$rMSE <- sqrt(trt_effect_output$sd_pi_hat^2 + trt_effect_output$bias^2)
    
    # Define the file name
    file_name <- paste0("MainEffectIndifference_results_", st, ".csv")
    
    # Create the file path
    file_path <- file.path(folder_path, file_name)
    
    # Write the data to a CSV file
    write.csv(trt_effect_output, file = file_path, row.names = FALSE)
    cat("CSV file saved to", file_path, "\n")
    
    ## DTR
    dtr_effect_output <- data.frame(dtr = c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11"),
                                    true_dtr = expected_pref,
                                    DTR_hat = apply(DTR_hat_bjsm,2, mean, na.rm = TRUE),
                                    sd_DTR_hat = apply(DTR_hat_bjsm,2,sd, na.rm = TRUE),
                                    avg_sd_DTR_hat = sqrt(apply(DTR_hat_bjsm_var,2, mean, na.rm = TRUE)))
    dtr_effect_output$bias <- dtr_effect_output$DTR_hat - dtr_effect_output$true_dtr
    dtr_effect_output$rMSE <- sqrt(dtr_effect_output$sd_DTR_hat^2 + dtr_effect_output$bias^2)
    
    # Define the file name
    file_name <- paste0("DTR_results_", st, ".csv")
    
    # Create the file path
    file_path <- file.path(folder_path, file_name)
    
    # Write the data to a CSV file
    write.csv(dtr_effect_output, file = file_path, row.names = FALSE)
    cat("CSV file saved to", file_path, "\n")
    
    
    ## coverage 
    params <- c("PiA", "PiB", "PiAC", "PiAD", "PiBC", "PiBD", "AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11", "alphaP1A",  "alphaP1B", "alphaP2C", "alphaP2D", "Beta1A", "Beta1B")
    coverage_output <- data.frame(Parameter = params, 
                                  Coverage = apply(ci_hat, 2, mean, na.rm = TRUE))
    
    # Define the file name
    file_name <- paste0("CI_results_", st, ".csv")
    
    # Create the file path
    file_path <- file.path(folder_path, file_name)
    
    # Write the data to a CSV file
    write.csv(coverage_output, file = file_path, row.names = FALSE)
    cat("CSV file saved to", file_path, "\n")
    
    ## Sims skipped
    num_skip_total <- sum(num_skip)
    mean_skip <- mean(num_skip)
    skip_df <- data.frame(num_sim = n.sim, num_skip_total = num_skip_total, avg_num_skip = mean_skip)
    
    # Define the file name
    file_name <- paste0("Num_Sim_Skip_Summary_", st, ".csv")
    
    # Create the file path
    file_path <- file.path(folder_path, file_name)
    
    # Write the data to a CSV file
    write.csv(skip_df, file = file_path, row.names = FALSE)
    cat("CSV file saved to", file_path, "\n")
  
}
