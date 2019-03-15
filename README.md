# How to run AWPE survey simulations

The simulations simulate the outcome of surveys monitoring reproductive success of American White Pelicans (AWPE).  
Simulations results and analyses form the basis of the report _Protocol for Monitoring Reproductive Success of  
Western American White Pelicans: Phase 1: Evaluating Study Design_, by N. Nur, C. Moulton, G. Block, and  
L. Salas. (2018. Final Report to United States Fish & Wildlife Service, Inventory & Monitoring Initiative,  
Pacific Southwest Region. Point Blue Conservation Science, 3820 Cypress Dr., Petaluma, CA, 94954, USA.)  

These surveys estimate the number of active nests (N), number of fledglings (F) and ratio of fledglings to nests (R).  
Thus, we simulate true N, true R (number of fledglings per nest), and thus true F, which is the product of N x R.  
We then simulate estimates of N, F, and R and evaluate precision in estimating R under different sampling designs.

The simulations use a set of R objects we constructed for general purposes of power analysis simulations. We are  
not including the code for those objects here because the license for this repository does not meet the requirements  
for us to publish these objects. In this repository we include the following:
1. The complete logic of our approach, including the code that uses the power analysis objects
2. The simulation definition files (.yaml files in the SimulationDefinitions folder) that specify what to simulate and how
3. The resulting synthetic data from a run of 1,000 simulations (.RData files in the SimulationData folder)
4. The code scripts that process these synthetic data and generates the results (see 3_compileResults.R) used in  
the Stata analyses and presented in the report (see "To run simulations, extract results and summarize" below)
5. An additional script file that provides further details (visualizeR_seR.R)

## To run simulations, extract results, and summarize

The script 1_runSurveySimulations.R contains all the code to run simulations. It sources the power analysis objects and  
then processes the simulation scenarios specified in the yaml files found in the SimulationDefinitions folder. Note that  
the names of these scenarios refer to names of the .yaml files. Refer to the How_To file in the SimulationDefinitions folder  
to understand how these simulations are defined in each of the .yaml files.

Once the simulations are completed, results are extracted with the script 2_processSimulations.R. This script uses the  
instructions in the .yaml files in the folder AnalysisDefinitions to understand what information to extract from the  
simulated data. The user is welcomed to try running this code with the simulation data in the SimulationData folder.  
Just make sure to adjust the paths to files as needed.

Once the results are extracted, these are summarized by the script 3_compileResults.R, and exported to Stata using  
a comma-separated-values file as output. 

### About this work

The simulations and calculations were developed by Nadav Nur and Leo Salas, Point Blue Conservation Science, per contract  
with the U.S. Fish and Wildlife Service - Region 8, on behalf of the American White Pelican Reproductive Success Protocol  
project team. 
  
For any further inquiries, contact Dr. Nadav Nur at nnur@pointblue.org
