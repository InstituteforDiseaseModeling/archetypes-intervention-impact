Code originally by Sam Bhatt, refactored by Amelia Bertozzi-Villa, for generating ITN coverage "cubes". 

* sam_orig is the code directly from sam, unchanged except for some comments. 
Script order:
1. Create_database_ITN access_covariates_means.R (generates the household level use and access metrics used by INLA)
2. acc_new_deviation.r (Collates covariates and runs INLA command for "access deviation"-- that is, deviation from the national mean of access)
3. gap_new_deviation.r (runs INLA command for net use gap)
4. ITN prediction.r (uses the above outputs to generate cubes)
note: in Sam's scripts, "/home/drive/" refers to the Z:/ driv (smb://map-fs1.ndph.ox.ac.uk/map_data) and /home/backup/ refers to smb://map-fs1.ndph.ox.ac.uk/map_pit/sam/bld1/bras2280/.

* run_on_gcloud is a minimal refactoring of sam's code to run on the google cloud framework-- mostly filepath updates and some function refactoring to get the process running

* amelia_refactor is a complete refactor of the process, but preserving all bugs such that it replicates the outputs of run_on_glcoud exactly. 

Having done this work within the malaria-atlas-project repository, I will now move amelia_refactor into its own repository for bug testing and further updates. 

Last Updated: ABV, 3 June 2019
