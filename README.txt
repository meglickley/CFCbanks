Lickley et al. (in Revisions) Quantifying contributions of chlorofluorocarbon banks to emissions and impacts on the ozone layer and climate

Code for generating figures in main text. 

Contained in Code Directory

estimating_RFandDE.m: matlab script to that simulates the bottom up accounting of banks by equipment type and generates joint prior distributions of RF (release fraction) and DE (direct emissions)

MainScript_multi.m: matlab script that runs the BPE model. 

simulation_model.m: matlab function that is called from MainScript_multi.m

Bank_Emiss.m:  matlab function that is called from simulation_model.m

Make_Figures.m: matlab script that makes the figures after the prior and posterior samples have been generated for each scenario and gas. 

SPARC_lifetimes.mat: file containing SPARC atmospheric lifetimes for each gas. 

Figures_data: a folder containing data from previous published assessments that are used to create figures in the main text  

boundedline.m: a function used to create the figures. 

CFC11, CFC12, CFC113 folders: contains Input and Output Folders.  Each Input folder for each gas contains the following: 

- a5_na5.mat: a file containing production data used in MainScript_multi.m
- AFEAS_cfcXXproduction.mat: a file containing AFEAS production data by equipment type
- WMO2002.mat: a file containing production data used in MainScript_multi.m
- wmo2018.mat: a file containing atmospheric concentrations 


%%%%%%%%%%%%%%%%%%  INSTRUCTION FOR RUNNING THE BPE MODEL %%%%%%%%%%%%%%%%%%%%%%

Step 1:  Open estimating_RFandDE.m script.  Change HomeDir to specify location of file. Run this script.  This will generate joint priors for RF and DE that is then stored in the Input Folder. 

Step 2:  Open MainScript_multi.m.  Change HomeDir to specify location of file.  Run the following scenarios by changing the variables in the ‘Specify Scenarios’ Section: 
		- meanSPARC = 1, CFC11 = 1, all other variables = 0
		- meanSPARC = 1, CFC11 = 1, fugitive_emissions = 1, all other variables = 0
		- constant_LT = 1, CFC11 = 1, all other variables = 0
		- constant_LT_scen2 = 1, CFC11 = 1, all other variables = 0
		- meanSPARC = 1, CFC12 = 1, all other variables = 0
		- constant_LT = 1, CFC12 = 1, all other variables = 0
		- meanSPARC = 1, CFC113 = 1, all other variables = 0
		- constant_LT = 1, CFC113 = 1, all other variables = 0

Step 3:  Run the Make_Figures.m script to generate figures and populate table 1 values 


