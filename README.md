# PRPP-SMART BJSM
Companion code for "Bayesian Estimation of Dynamic Treatment Regimens in a Partially Randomized Patient Preference, Sequential, Multiple-Assignment, Randomized Trial".

## Usage Note
The code in this repository is in active development. To view or use stable code, see the appropriate releases:
- [v1.0.0](../../releases/tag/v1.0.0): Release accompanying initial submission to _Statistics in Biopharmaceutical Research_.

## File Descriptions
- [DataGeneration.R](DataGeneration.R): Main function used to generate data from a two-stage PRPP-SMART with binary end-of-stage outcome.
- [BJSM.R](BJSM.R): Simulation code for Bayesian Joint Stage Model (BJSM) to estimate embedded dynamic treatment regimes (DTRs) in PRPP-SMART under all scenarios/sub-scenarios considered in the manuscript. Note, simulations are reccomened to be run on a high-performance cluster rather than a laptop.
- [nsim_toget_500.R](nsim_toget_500.R): Functions used to determine the number of simulations needed to run per scenario/sub-scenario in order to achieve 500 total simulations. 

## Folder Descriptions
- The [SubScenarios](SubScenarios)  folder contains the data generation parameters, including linkage parameter values and treatment outcome rates, that define each sub-scenario analyzed in the manuscript "Bayesian Estimation of Dynamic Treatment Regimens in a Partially Randomized Patient Preference, Sequential, Multiple-Assignment, Randomized Trial." Each sub-scenario is specified in its own .R file. For a detailed description of the sub-scenarios, refer to the manuscript.


