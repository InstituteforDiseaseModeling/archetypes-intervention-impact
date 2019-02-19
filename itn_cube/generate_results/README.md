Code originally by Sam Bhatt, mid-refactor by Amelia Bertozzi-Villa, for generating ITN coverage "cubes". 

sam_orig script order:

1. Create_database_ITN access_covariates_means.R (generates the household level use and access metrics used by INLA)
2. acc_new_deviation.r (Collates covariates and runs INLA command for "access deviation"-- that is, deviation from the national mean of access)
3. gap_new_deviation.r (runs INLA command for net use gap)
4. ITN prediction.r (uses the above outputs to generate cubes)

note: in Sam's scripts, "/home/drive/" refers to the Z:/ driv (smb://map-fs1.ndph.ox.ac.uk/map_data) and /home/backup/ refers to smb://map-fs1.ndph.ox.ac.uk/map_pit/sam/bld1/bras2280/.


amelia_refactor script order and dependencies:

create_database_refactored.r
	create_database_functions.r
acc_deviation_refactored.r
	acc_deviation_functions.r
	algorithm.V1.R
	INLAFunctions.R
[gap_new_devaition.r still to-be-commented/refactored from sam's]
itn_prediction_refactored.r
	itn_prediction_functions.r
scratch_itn_adjust.r (notes from my and sam's meeting in London)

Last Updated: ABV, 19 Feb 2019