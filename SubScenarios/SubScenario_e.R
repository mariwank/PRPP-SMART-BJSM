# Sub-Scenario (e) 
# additive preference effect

# additive effect 
preferAEff <- 0.02
preferBEff <- 0.1
preferCEff <- 0.03
preferDEff <- 0.01

# First stage treatment response rates
pi_A <- 0.6       # stage 1 response rate to randomize A: Pr(R=1|T1=A,P1=0)
pi_B <- 0.45      # stage 1 response rate to randomize B: Pr(R=1|T1=B,P1=0)
pi_A1 <- pi_A +  preferAEff    # stage 1 response rate to prefer A: Pr(R=1|T1=A,P1=1)
pi_B1 <- pi_B +  preferBEff    # stage 1 response rate to prefer B: Pr(R=1|T1=B,P1=1)


pi_AC <- 0.7   # Second stage response rate of non-responders to randomized A who receive randomized C in the second stage: Pr(Y=1|T1=A,P1=0,NR,P2=0,T2=C)        
pi_AD <- 0.55   # Second stage response rate of non-responders to randomized A who receive randomized D in the second stage: Pr(Y=1|T1=A,P1=0,NR,P2=0,T2=D)
pi_BC <- 0.35    # Second stage response rate of non-responders to randomized B who receive randomized C in the second stage: Pr(Y=1|T1=B,P1=0,NR,P2=0,T2=C)
pi_BD <- 0.2   # Second stage response rate of non-responders to randomized B who receive randomized D in the second stage: Pr(Y=1|T1=B,P1=0,NR,P2=0,T2=D)

pA0A <- pi_A       # Second stage response rate of responders to randomized A
pB0B <- pi_B       # Second stage response rate of responders to randomized B
pA1A <- pi_A1       # Second stage response rate of responders to preferred A
pB1B <- pi_B1      # Second stage response rate of responders to preferred B
pA0C1 <- pi_AC + preferCEff   # Second stage response rate of non-responders to randomized A who receive preferred C in the second stage
pA0D1 <- pi_AD + preferDEff   # Second stage response rate of non-responders to randomized A who receive preferred D in the second stage
pA1C0 <- pi_AC + preferAEff     # Second stage response rate of non-responders to preferred A who receive randomized C in the second stage
pA1D0 <- pi_AD + preferAEff  # Second stage response rate of non-responders to preferred A who receive randomized D in the second stage
pA1C1 <- pi_AC +  preferAEff + preferCEff # Second stage response rate of non-responders to preferred A who receive preferred C in the second stage
pA1D1 <- pi_AD +  preferAEff + preferDEff   # Second stage response rate of non-responders to preferred A who receive preferred D in the second stage
pB0C1 <- pi_BC + preferCEff    # Second stage response rate of non-responders to randomized B who receive preferred C in the second stage
pB0D1 <- pi_BD + preferDEff  # Second stage response rate of non-responders to randomized B who receive preferred D in the second stage
pB1C0 <- pi_BC + preferBEff   # Second stage response rate of non-responders to preferred B who receive randomized C in the second stage
pB1D0 <- pi_BD + preferBEff   # Second stage response rate of non-responders to preferred B who receive randomized D in the second stage
pB1C1 <- pi_BC +  preferBEff + preferCEff  # Second stage response rate of non-responders to preferred B who receive preferred C in the second stage
pB1D1 <- pi_BD +  preferBEff + preferDEff  # Second stage response rate of non-responders to preferred B who receive preferred D in the second stage

