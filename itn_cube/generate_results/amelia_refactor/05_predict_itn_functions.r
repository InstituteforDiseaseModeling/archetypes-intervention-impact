###############################################################################################################
## 05_predict_itn_functions.r
## Amelia Bertozzi-Villa
## May 2019
## 
## Functions to accompany 05_predict_itn.r-- prepping survey data to go into the ITN cube model. 
##############################################################################################################

Inv.IHS <- function(x, theta){  # reverse IHS transformation
  (1/theta)*sinh(theta * x)
}
