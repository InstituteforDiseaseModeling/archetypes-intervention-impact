###############################################################################
## check_file_similarity.r
## Amelia Bertozzi-Villa
## April 2019
## 
## I'm in the process of refactoring Sam Bhatt's ITN Cube code. This script 
## compares original outputs to my new outputs to ensure the refactor is 
## producing the same results.
##############################################################################

library(data.table)

check_sameness <- function(orig, new, sameness_cutoff=1e-10){
  if (nrow(orig)!=nrow(new)){
    print("Datasets have varying lengths!")
  }else{
    for (this_name in names(orig)){
      
      print(paste("assessing column", this_name))
      
      mismatch <- copy(orig)
      setnames(mismatch, this_name, "old_val")
      mismatch$new_val <- new[[this_name]]
      
      mismatch <- mismatch[old_val!=new_val]
      #print(mismatch)

      if (nrow(mismatch)>0){
        print(paste("discrepancy in", this_name))
        
        mismatch$diff <- abs(mismatch$old_val - mismatch$new_val)
        max_diff <- max(mismatch$diff)

        if (max_diff>sameness_cutoff){
          print(paste("columns differ by at most", max_diff, "-- Explore further"))
          print(paste(nrow(mismatch[diff>sameness_cutoff]), "rows are above cutoff."))
          print(summary(mismatch))
        }else{
          print(paste("columns differ by at most", max_diff, "-- No issue"))
        }
      }
    }
  }
}

 



