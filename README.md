# PRPP-SMART BJSM
Companion code for "Bayesian Estimation of Dynamic Treatment Regimens in a Partially Randomized Patient Preference, Sequential, Multiple-Assignment, Randomized Trial".

## Usage Note
The code in this repository is in active development. To view or use stable code, see the appropriate releases:
- [v1.0.0](../../releases/tag/v1.0.0): Release accompanying initial submission to _Statistics in Biopharmaceutical Research_.

## File Descriptions
- [DataGeneration.R](DataGeneration.R): Main function used to generate data from a two-stage PRPP-SMART with binary end-of-stage outcome.
- [BJSM.R](BJSM.R): Simulation code for the Bayesian Joint Stage Model (BJSM) to estimate embedded dynamic treatment regimes (DTRs) in PRPP-SMART under all scenarios and sub-scenarios where BJSM assumptions are satisfied, as described in the manuscript. Simulations are recommended to be run on a high-performance computing cluster rather than a laptop.
- [BJSM_AssumptionsViolated.R](BJSM_AssumptionsViolated.R): Simulation code for the Bayesian Joint Stage Model (BJSM) to estimate embedded dynamic treatment regimes (DTRs) in PRPP-SMART under scenarios where the assumptions of BJSM are **violated**. While the BJSM model code remains the same across [BJSM.R](BJSM.R) and [BJSM_AssumptionsViolated.R](BJSM_AssumptionsViolated.R), the evaluation and output differ. Specifically, when BJSM assumptions are violated, linkage parameters cannot be evaluated. To maintain clean code organization, sub-scenarios where BJSM assumptions are satisfied and violated are handled in separate files. Simulations are recommended to be run on a high-performance computing cluster rather than a laptop.
- [nsim_toget_500.R](nsim_toget_500.R): Functions used to determine the number of simulations needed to run per scenario/sub-scenario in order to achieve 500 total simulations. 

## Folder Descriptions
- The [SubScenarios](SubScenarios)  folder contains the data generation parameters, including linkage parameter values and treatment outcome rates, that define each sub-scenario analyzed in the manuscript "Bayesian Estimation of Dynamic Treatment Regimens in a Partially Randomized Patient Preference, Sequential, Multiple-Assignment, Randomized Trial." Each sub-scenario is specified in its own .R file. For a detailed description of the sub-scenarios, refer to the manuscript.


