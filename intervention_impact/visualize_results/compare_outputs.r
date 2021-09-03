###############################################################################################################
## compare_outputs.r
## Amelia Bertozzi-Villa
## January 2020
## 
##############################################################################################################

library(lattice)
library(data.table)
library(raster)
library(rasterVis)
library(PNWColors)
library(MapSuite)
library(gridExtra)
library(RColorBrewer)

rm(list=ls())

func_dir <- file.path("~/repos/malaria-atlas-project/intervention_impact/visualize_results")
setwd(func_dir)
source("pr_to_r0.r")
source("map_ii_functions.r")

old_dir <- "20191008_replicate_megatrends/results/rasters"
new_dir <- "20191009_mega_era5_new_arch/results/megatrend_ii_rasters"
suffix <- ""
base_dir <- file.path(Sys.getenv("HOME"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/intervention_impact")

## color palettes and breaks
color_vals <- generate_full_pal()
prev_cols <- c(color_vals$pal[1:10], rev(brewer.pal(11, "RdYlBu")))
prev_breaks <- c(color_vals$breaks[1:9], seq(0.005, 1, length.out = 12))

level_palettes <- list(par=list(pal=wpal("seaside"),
                                breaks=NA,
                                zscaleLog=T,
                                name="Pop at Risk"),
                       pfpr=list(pal=prev_cols,
                                 breaks=prev_breaks,
                                 name="PfPR"),
                       repro_number=list(pal=c("#e0e0e0", terrain.colors(31)[1:30]),
                                                   breaks=c(0,1,seq(1.5, 5, length.out=25), seq(5.1, 80, length.out=4)),
                                                   zscaleLog=F,
                                         name="R0")
                     )

divpal <- c(pnw_palette("Lake", 60)[10:59],  rev(pnw_palette("Shuksan", 35))[3:17], rev(pnw_palette("Starfish", 75))[1:35])
diff_palettes <- list(par=list(pal=divpal,
                               breaks=(seq(-250000, 250000, 10000)),
                               zscaleLog=F),
                      pfpr=list(pal=divpal,
                                breaks=seq(-0.5, 0.5, 0.025),
                                zscaleLog=F),
                      repro_number=list(pal=divpal,
                                        breaks=seq(-0.5, 0.5, 0.025),
                                        zscaleLog=F)
                      )


### Intervention Labels  -----------------------------------------------------

ints <- fread(file.path(base_dir, new_dir, "../../input/interventions.csv"))
ints <- dcast.data.table(ints, int_id ~ int, value.var = "cov")

make_int_label <- function(row){
 
  if (!is.na(row$al_cm)){
    label <- paste0(row$itn, "/", row$irs, "/", row$al_cm, "% ITN/IRS/AL CM")
  }else{
    label <- paste0(row$itn, "/", row$irs, "/", row$dp_cm, "% ITN/IRS/DP CM")
  }
  
  other_labels <- c(mAb="mAb", pev="PEV", tbv="TBV", atsb="ATSB Kill Rate")
  
  for (var in names(other_labels) ){
    if (!is.na(row[[var]])){
      label <- paste0(label, " + ", row[[var]], "% ", other_labels[[var]])
    }
  }

  return(label)
  
}

labs <- sapply(1:nrow(ints), function(idx){
  return(make_int_label(ints[idx]))
})

labs <- c("Baseline", labs)





## plot

make_plots <- F

for (metric in names(level_palettes)){
  
  print(paste("Comparing for", metric))
  
  tifname <- paste0(metric, "_Africa.tif")
  this_level_pal <- level_palettes[[metric]]
  this_diff_pal <- diff_palettes[[metric]]
  
  if (!file.exists(file.path(base_dir, new_dir, suffix, tifname)) | !file.exists(file.path(base_dir, old_dir, suffix, tifname))){
    warning("One or more input files not found! Skipping this metric")
    next  
  }
  
  new <- brick(file.path(base_dir, new_dir, suffix, tifname))
  old <- brick(file.path(base_dir, old_dir, suffix, tifname))
  diff <- old-new
  
  if (metric=="par"){
    par_new <- cellStats(new, sum)/1000000
    par_old <- cellStats(old, sum)/1000000
    par_diff <- par_old-par_new
    perc_diff <- par_diff/par_old*100
    
    par_dt <- data.table(int_number = 0:(length(labs)-1),
                         int_name=labs,
                         par_old=round(par_old, 2),
                         par_new=round(par_new, 2),
                         par_diff=round(par_diff, 2),
                         par_perc_diff=round(perc_diff,2))
    
    par_dt[par_old==0 & par_diff==0, par_perc_diff:=0]
    
    par_dt_long <- melt(par_dt, id.vars="int_name", measure.vars = c("par_diff", "par_perc_diff"))
    par_dt_long[, variable:=ifelse(variable %like% "perc", "PAR Percent Difference (%)", "PAR Difference (Millions)")]
    
    par_dens_plot <- ggplot(par_dt_long, aes(x=value, color=variable, fill=variable)) +
                      geom_density(alpha=0.5) +
                      facet_grid(variable~.) +
                      theme_minimal() +
                      theme(legend.position = "none") +
                      labs(x="Value", 
                           y="Density")
    
    pdf(file.path(base_dir, new_dir, suffix, "comparison_par_density.pdf"), width=8, height=5)
      print(par_dens_plot)
    graphics.off()
    
    names(par_dt) <- c("Intervention", "Description", "Old PAR", "New PAR", "Difference (Millions)", "% Difference")
    write.csv(par_dt, file = file.path(base_dir, new_dir, suffix, "comparison_par.csv"), row.names = F)
  }
  
  
  if (make_plots){
    print("Plotting")
    nplots <- nlayers(new)
    pdf(file.path(base_dir, new_dir, suffix, paste0("comparison_", metric, ".pdf")), width=8, height=3)
    
    for (i in 1:nplots){
      
      if (metric=="par"){
        new_plot <- levelplot(new[[i]]+1, par.settings=rasterTheme(this_level_pal$pal),
                              zscaleLog=T, xlab=NULL, ylab=NULL, margin=F, scales=list(draw=F), main=paste("New:", round(par_new[[i]], 1), "M"), cex.main=0.75)
        
        old_plot <- levelplot(old[[i]]+1, par.settings=rasterTheme(this_level_pal$pal), 
                              zscaleLog=T, xlab=NULL, ylab=NULL, margin=F, scales=list(draw=F), main=paste("Old:", round(par_old[[i]], 1), "M"), cex.main=0.75)
        
        full_plot <- grid.arrange(new_plot, old_plot, ncol=2, top=paste0(this_level_pal$name, ": ", labs[i]))
        print(full_plot)
  
      }else{
        new_plot <- levelplot(new[[i]], par.settings=rasterTheme(this_level_pal$pal), at=this_level_pal$breaks,
                              zscaleLog=F, xlab=NULL, ylab=NULL, margin=F, scales=list(draw=F), main=paste0("New (PAR ", round(par_new[[i]], 1), " M)"), cex.main=0.75)
        
        old_plot <- levelplot(old[[i]], par.settings=rasterTheme(this_level_pal$pal), at=this_level_pal$breaks,
                              zscaleLog=F, xlab=NULL, ylab=NULL, margin=F, scales=list(draw=F), main=paste0("Old (PAR ", round(par_old[[i]], 1), " M)"), cex.main=0.75)
        
        diff_plot <- levelplot(diff[[i]], par.settings=rasterTheme(this_diff_pal$pal), at=this_diff_pal$breaks,
                               zscaleLog=F, xlab=NULL, ylab=NULL, margin=F, scales=list(draw=F), main=paste0("Difference: ", round(par_diff[[i]], 1), " M"), cex_main=0.75 )
        full_plot <- grid.arrange(new_plot, old_plot, diff_plot, ncol=3, top=paste0(this_level_pal$name, ": ", labs[i]))
        print(full_plot)
        
      }
    }
  }
  
  graphics.off()
  
}



