# Community_Gap-Filling
Source code and input files for reproducing the results presented in "A gap-filling algorithm for prediction of metabolic interactions in microbial communities".

## Input Files
The folder input files contains all the necessary files for the use of the community gap-filling method.
* The folder input files/database contains the database file (BiGG.mat).
* The folder input files/toy_E.coli_community contains the models for the *E. coli* glucose utilizer (model1.mat) and the *E. coli* acetate utilizer (model2.mat), as well as the file with the permitted exchanges for the toy E. coli commuinity (community_media.txt). The folder also contains the used model for the *E. coli* glycerol utilizer (model3.mat).
* The folder input files/gut-microbiome_community contains the models for *B. adolescentis* ATCC 15703 (model1.mat) and *F. prausnitzii* A2-165 (model2.mat), as well as the file with the permitted exchanges for the commuinity (community_media.txt).
* The folder input files/ACT-3_community contains the models for *Dehalobacter* sp. CF (model1.mat) and *Bacteroidales* sp. CF50 (model2.mat), as well as the file with the permitted exchanges for the commuinity (community_media.txt).

## Source Code 
The script main.m calls the function make_member.m in order to make the necessary community members that are merged into a community model. It also uses a file with the community media composition in order to constrain the exchange fluxes of the community model. Then, the script creates an MILP from the community and finds multiple alternative solutions of the problem. The script gives as an output an .xlsx file with the reaction fluxes calculated for the best solution, an .xlsx file with the fluxes of the exchange reactions for the 10 best solutions, and an .xlsx file with the fluxes of the reactions that have been added from the database to the community model for the 10 best solutions. The inputs of main.m are defined in the following lines:
* Line 6: Give the database in a .mat file.
* Lines 7 and 8: Give the metabolic reconstructions for the microorganisms of the community in .mat files.
* Line 9: Give the media composition for the community in a .txt file.
* Lines 23 and 24: Set the lower bound for the flux through the biomass reaction of each organism.
* Lines 105-113: Give the parameters for the calculation of multiple solutions of the MILP with the CPLEX procedure "populate".

All simulations were performed with CobraToolbox v.3, CPLEX 12.8, and MATLAB R2017b.
