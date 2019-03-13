# How to run AWPE survey simulations

The simulations intend to resemble the outcome of field surveys of American White Pelicans (AWPE). These surveys  
estimate the number of adults (N) and ratio of chicks to adults (R). Thus, we simulate N, and we simulate fecundity  
(the number of chicks per adult - F), from which we calculate a simulated R.

The simulations use a set of R objects we constructed for general purposes of power analysis simulations. We are  
not including the code for those objects here because the license for this repository does not meet the requirements  
for us to publish these objects. However, we are including here:
1. The complete logic of our approach, including the code that uses the power analysis objects
2. The simulation definition files (.yaml files in the SimulationDefinitions folder) that specify what to simulate and how
3. The resulting synthetic data from a run of 1,000 simulations (.RData files in the SimulationData folder)
4. The code scripts that process these synthetic data and generate the results used in the Stata analyses and presented in the report
5. An additional script file that shows justifications for some of our decisions (visualizeR_seR.R)

## To run simulations, extract results, and summarize

The script 1_runSurveySimulations.R contains all the code to run simulations. It sources the power analysis objects and  
then processes the simulation scenarios specified therein. Note that the names of these scenarios refer to names of .yaml  
files in the SimulationDefinitions folder. Refer to the How_To file in that folder to understand how these simulations  
are defined in each of the .yaml files.

Once the simulations complete, results are extracted with the script 2_processSimulations.R. This script uses the  
instructions in the .yaml files in the folder AnalysisDefinitions to understand what information to extract from the  
simulated data. The user is welcomed to try running this code with the simulation data in the SimulationData folder.  
Just make sure to adjust the paths to files as needed.

Once the results are extracted, these are summarized by the script compileResults.R, and exported to Stata using  
a comma-separated file as output. 

[Best to read this file by copying its content in a Markdown editor, for example: https://dillinger.io/]

