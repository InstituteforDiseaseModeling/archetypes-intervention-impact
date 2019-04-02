library(data.table)

main_dir <- "/Users/bertozzivill/Desktop"

orig <- fread(file.path(main_dir, "users_amelia_itn_cube_create_database_output_ITN_final_clean_access_18March2019.csv"))
new <- fread(file.path(main_dir, "users_amelia_itn_cube_create_database_output_ITN_final_clean_access_18March2019_COMPARE.csv"))

sameness_cutoff <- 1e-10 #past this point, values are similar to a rounding error
 
if (nrow(orig)!=nrow(new)){
  print("Datasets have varying lengths!")
}else{
  for (this_name in names(orig)){
    mismatch <- orig[orig[[this_name]]!=new[[this_name]]]
    if (nrow(mismatch)>0){
      print(paste("discrepancy in", this_name))
      orig$new_val <- new[[this_name]]
      orig$diff <- abs(orig[[this_name]] - orig[["new_val"]])
      max_diff <- max(orig$diff)
      if (max_diff>sameness_cutoff){
        print(paste("columns differ by at most", max_diff, "-- Explore further"))
      }else{
        print(paste("columns differ by at most", max_diff, "-- No issue"))
      }
    }
  }
}


