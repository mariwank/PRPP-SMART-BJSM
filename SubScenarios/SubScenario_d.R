# Sub-Scenario (d)

# Data gen parameters 

beta1_link <- c(1.01,0.88)                          # Beta1A (responders to A) Beta1B (resopnders to B)       
alphaP_link <- c(1.30, 1.39, 1.40, 2.15)     # alpha_p1A (1st stage preference, T1=A), alpha_p1B (1st stage preference, T1=B), alpha_p2C (2nd stage preference, T2=C), alpha_p2D (2nd stage preference, T2=D)

# First stage treatment response rates
pi_A <- 0.6                                          # stage 1 response rate to randomize A: Pr(R=1|T1=A,P1=0)
pi_B <- 0.45                                         # stage 1 response rate to randomize B: Pr(R=1|T1=B,P1=0)
pi_A1 <- pi_A*alphaP_link[1]                         # stage 1 response rate to prefer A: Pr(R=1|T1=A,P1=1)
pi_B1 <- pi_B*alphaP_link[2]                         # stage 1 response rate to prefer B: Pr(R=1|T1=B,P1=1)

# Second stage treatment response rates
pi_AC <- 0.475                                        # Second stage response rate of non-responders to randomized A who receive randomized C in the second stage: Pr(Y=1|T1=A,P1=0,NR,P2=0,T2=C)        
pi_AD <- 0.225                                         # Second stage response rate of non-responders to randomized A who receive randomized D in the second stage: Pr(Y=1|T1=A,P1=0,NR,P2=0,T2=D)
pi_BC <- 0.317                                         # Second stage response rate of non-responders to randomized B who receive randomized C in the second stage: Pr(Y=1|T1=B,P1=0,NR,P2=0,T2=C)
pi_BD <- 0.1                                         # Second stage response rate of non-responders to randomized B who receive randomized D in the second stage: Pr(Y=1|T1=B,P1=0,NR,P2=0,T2=D)
pA0A <- pi_A * beta1_link[1]                         # Second stage response rate of responders to randomized A
pB0B <- pi_B * beta1_link[2]                         # Second stage response rate of responders to randomized B
pA1A <- pi_A * alphaP_link[1] * beta1_link[1]        # Second stage response rate of responders to preferred A
pB1B <- pi_B * alphaP_link[2] * beta1_link[2]        # Second stage response rate of responders to preferred B
pA0C1 <- alphaP_link[3] * pi_AC                      # Second stage response rate of non-responders to randomized A who receive preferred C in the second stage
pA0D1 <- alphaP_link[4] * pi_AD                      # Second stage response rate of non-responders to randomized A who receive preferred D in the second stage
pA1C0 <- alphaP_link[1] * pi_AC                      # Second stage response rate of non-responders to preferred A who receive randomized C in the second stage
pA1D0 <- alphaP_link[1] * pi_AD                      # Second stage response rate of non-responders to preferred A who receive randomized D in the second stage
pA1C1 <- alphaP_link[1] * alphaP_link[3] * pi_AC     # Second stage response rate of non-responders to preferred A who receive preferred C in the second stage
pA1D1 <- alphaP_link[1] * alphaP_link[4] * pi_AD     # Second stage response rate of non-responders to preferred A who receive preferred D in the second stage
pB0C1 <- alphaP_link[3] * pi_BC                      # Second stage response rate of non-responders to randomized B who receive preferred C in the second stage
pB0D1 <- alphaP_link[4] * pi_BD                      # Second stage response rate of non-responders to randomized B who receive preferred D in the second stage
pB1C0 <- alphaP_link[2] * pi_BC                      # Second stage response rate of non-responders to preferred B who receive randomized C in the second stage
pB1D0 <- alphaP_link[2] * pi_BD                      # Second stage response rate of non-responders to preferred B who receive randomized D in the second stage
pB1C1 <- alphaP_link[2] * alphaP_link[3] * pi_BC     # Second stage response rate of non-responders to preferred B who receive preferred C in the second stage
pB1D1 <- alphaP_link[2] * alphaP_link[4] * pi_BD     # Second stage response rate of non-responders to preferred B who receive preferred D in the second stage


