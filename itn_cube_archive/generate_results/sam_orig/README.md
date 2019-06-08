Sam's original ITN cube code.


Script order:
1. Create_database_ITN access_covariates_means.R (generates the household level use and access metrics used by INLA)
2. acc_new_deviation.r (Collates covariates and runs INLA command for "access deviation"-- that is, deviation from the national mean of access)
3. gap_new_deviation.r (runs INLA command for net use gap)
4. ITN prediction.r (uses the above outputs to generate cubes)
note: in Sam's scripts, "/home/drive/" refers to the Z:/ driv (smb://map-fs1.ndph.ox.ac.uk/map_data) and /home/backup/ refers to smb://map-fs1.ndph.ox.ac.uk/map_pit/sam/bld1/bras2280/.
