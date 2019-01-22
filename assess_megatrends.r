## -----------------------------------------------------------------------------------------------------------------
# Megatrends: Characterizing residual transmission
# assess_megatrends.r
# 
# Amelia Bertozzi-Villa, Institute for Disease Modeling, University of Oxford
# December 2018
# 
# As per the request of the WHO-SAGme, look at the megatrends + existing intervention forecasts and
# take a first pass at determining:
#   1. Whether pockets of residual transmission are determined by high initial transmission or something else;
#   2. If "something else", what mechanistic modeling scenarios result in higher residual transmission.
# 
## -----------------------------------------------------------------------------------------------------------------------

library(raster)
library(data.table)
library(stats)
library(ggplot2)
library(rasterVis)
library(gridExtra)

rm(list=ls())

root_dir <- ifelse(Sys.getenv("USERPROFILE")=="", Sys.getenv("HOME"))
base_dir <- file.path(root_dir, "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact")

main_dir <- file.path(base_dir, "/megatrends")
cluster_fname <- file.path(base_dir, "lookup_tables/interactions/africa_clusters_v4.tif")


anthro_endo_map <- data.table(cluster=1:6,
                              sitename=c("aba", "kananga", "kasama", "djibo", "gode", "moine"),
                              anthro=c(74.45, 65.02, 79.04, 76.6, 75, 75.78),
                              endo=c(80, 85, 80.38, 55.6, 50, 52.73))
anthro_endo_map[, human_indoor:= (anthro*endo)/100]

colors <- c("#00a08a", "#d71b5a", "#f2a200", "#f98400", "#902e57", "#5392c2")

# load in data, clip to africa
cluster_layer <- raster(cluster_fname)
megatrends_noints <- raster(file.path(main_dir, "pete_analysis", "actual_ssp2_base2000_2050.tif"))
megatrends_ints <- raster(file.path(main_dir, "pete_analysis", "actual_ssp2_base2016_2050.tif"))

int_type <- "itn_irs_act"

if (int_type=="itn_act"){
  interventions <- raster(file.path(main_dir, "abv_pfpr_africa_ITN.0.8..ACT.0.8..tif"))
  interventions_pete <- raster(file.path(main_dir, "pete_analysis", "actual_ssp2_2050_ITN80ACT80-14.tif"))
  int_label <- "ITN 80%, ACT 80%"
}else if (int_type=="itn_irs_act"){
  interventions <- raster(file.path(main_dir, "abv_pfpr_africa_ITN.0.8..IRS.0.8..ACT.0.8..tif"))
  interventions_pete <- raster(file.path(main_dir, "pete_analysis", "actual_ssp2_2050_ITN80IRS80ACT80-14.tif"))
  int_label <- "ITN 80%, IRS 80%, ACT 80%"
}else{
  stop(paste("unrecognized intervention type", int_type))
}

# ensure consistent extents in all rasters, and
# set near-zero areas to zero.
cutoff_pr <- 0.0001
cluster_layer <- crop(cluster_layer, interventions)
megatrends_orig <- copy(megatrends_noints)

rastlist <- list("megatrends_noints"=megatrends_noints, "megatrends_ints"=megatrends_ints,
                 "interventions"=interventions, "interventions_pete"=interventions_pete)

print("cropping and removing zeros")
for (rastname in names(rastlist)){
  rast <- rastlist[[rastname]]
  rast[rast<cutoff_pr] <- 0
  newrast <- crop(rast, cluster_layer)
  names(newrast) <- rastname
  assign(rastname, newrast)
}

# mask clusters to megatrends
cluster_layer <- raster::mask(cluster_layer, megatrends_noints)
names(cluster_layer) <- "cluster"

# my "interventions" dataset is unbounded by the "megatrends with interventions" dataset.
# I need to bound it to be comparable to Pete's 
# (and to explore the difference between bounded and unbounded residual transmission)

# TODO: find reason for differences between my and Pete's bounded rasters. 
bounded_interventions <-min(stack(interventions,megatrends_ints))
names(bounded_interventions) <- "bounded_interventions"
megatrends_diff <- interventions - bounded_interventions

# compare my bounded intervention to Pete's
abv_pg_diff <- bounded_interventions - interventions_pete
plot(abv_pg_diff>0.01 | abv_pg_diff<(-0.01), main=">1% difference between ABV and PG")


# spatially disaggregate pixels that track lookup table vs megatrends 
pdf(file.path(main_dir, "cluster_plots.pdf"), width=12, height=6)

for (cluster_idx in 1:6){
  this_color <- colors[[cluster_idx]]
  stacked_series <- stack(cluster_layer==cluster_idx,
                          cluster_layer==cluster_idx & bounded_interventions>cutoff_pr,
                          cluster_layer==cluster_idx & bounded_interventions>cutoff_pr & bounded_interventions<interventions)
  
  names(stacked_series) <- c("Full Cluster", "Residual Transmission", "Megatrends RT")
  cluster_plot <- levelplot(stacked_series, att="ID", col.regions=c("#A9A9A9", this_color),
                            xlab=NULL, ylab=NULL, scales=list(draw=F),
                            main = int_label, colorkey=F, margin=F)
  print(cluster_plot)
}
graphics.off()



# explore areas of residual transmission--------------------------


raster_to_dt <- function(rast){
  vals <- as.matrix(rast)
  vals <- data.table(id = which(!is.na(vals) & !is.infinite(vals)),
                     value = vals[!is.na(vals) & !is.infinite(vals)])
  vals[, type:= names(rast)]
  return(vals)
}

reduction_dt <- lapply(list(interventions, bounded_interventions, megatrends_ints, megatrends_noints, cluster_layer), raster_to_dt)
reduction_dt <- rbindlist(reduction_dt)
reduction_dt <- dcast.data.table(reduction_dt, id ~ type)
reduction_dt <- reduction_dt[complete.cases(reduction_dt)] 
reduction_dt <- reduction_dt[megatrends_noints>0]
reduction_dt <- merge(reduction_dt, anthro_endo_map, by="cluster", all=T)

reduction_dt[, cluster:=as.factor(cluster)]
reduction_dt[, eliminate:=floor(1-bounded_interventions)]
reduction_dt[, init_prev_class:= cut(megatrends_noints, breaks=seq(0,1,0.05),
                                     labels=c("0-5", "5-10", "10-15", "15-20", "20-25", "25-30", "30-35",
                                              "35-40", "40-45", "45-50", "50-55", "55-60", "60-65", "65-70",
                                              "70-75", "75-80", "80-85", "85-90", "90-95", "95-100"))]


summary_dt <- reduction_dt[,lapply(.SD, mean), by=list(sitename, cluster, init_prev_class)]
summary_dt <- merge(summary_dt, reduction_dt[, list(count=.N), by=list(cluster, init_prev_class)],
                    by=c("cluster", "init_prev_class"), all=T)
summary_dt[, count_class:=cut(count, breaks=quantile(count, probs=seq(0,1,0.2)),
                              labels=c("4-161", "162-1,743", "1,744-7,839", "7,840-14,003", "14,004-59,512"), include.lowest=T)]

elim_plot <- ggplot(summary_dt, aes(x=init_prev_class, y=human_indoor)) + 
            geom_point(aes(color=eliminate, size=count_class), shape=15) +
            scale_size_discrete(name="Pixel Count") + 
            scale_color_gradientn(colors=brewer.pal(7,"Spectral"), name="Proportion of\nPixels Eliminating") +
            theme(axis.text.x = element_text(angle=45, hjust=1)) +
            labs(x="Initial PfPR (%)",
                 y="Indoor Biting (%)",
                 title=paste("Probability of Elimination Under", int_label))

resid_plot <- ggplot(summary_dt, aes(x=init_prev_class, y=human_indoor)) + 
                geom_point(aes(color=bounded_interventions, size=count_class), shape=15) +
                scale_size_discrete(name="Pixel Count") + 
                scale_color_gradientn(colors=rev(brewer.pal(7,"Spectral")), name="Mean Final\nPrevalence") +
                theme(axis.text.x = element_text(angle=45, hjust=1)) +
                labs(x="Initial PfPR (%)",
                     y="Indoor Biting (%)",
                     title=paste("Residual Transmission Under", int_label))

pdf(file.path(main_dir, paste0("resid_by_endo_", int_type, ".pdf")), height=7, width=16)
  grid.arrange(grobs=list(elim_plot, resid_plot), layout_matrix=rbind(c(1,2)))
graphics.off()


png(file.path(main_dir, "megatrend_compare.png"), width=900, height=600)
ggplot(reduction_dt, aes(x=megatrends_noints, y=megatrends_ints)) +
  geom_point(alpha=0.25) +
  facet_wrap(~cluster) +
  geom_abline() +
  theme_minimal() +
  theme(legend.position="none") +
  xlim(0,1)  +
  labs(title="Megatrend (Base 2000) vs Megatrend (Base 2016), \n by Transmission Archetype",
       x="Megatrends, Base 2000",
       y="Megatrends, Base 2016")
graphics.off()

png(file.path(main_dir, "megatrend_int_impact.png"), width=900, height=600)
ggplot(reduction_dt[bounded_interventions>cutoff_pr], aes(x=megatrends_noints, y=bounded_interventions)) + 
        # geom_point(aes(y=megatrends_ints), alpha=0.1) + 
        geom_point(aes(color=cluster), alpha=0.5) + 
        scale_color_manual(values=colors) + 
        facet_wrap(~cluster) + 
        geom_abline() +
        theme_minimal() + 
        theme(legend.position="none") + 
        xlim(0,1)  + 
        labs(title="Megatrend (Base 2000) vs Megatrend Plus Interventions, \n by Transmission Archetype",
             x="Megatrend, Base 2000",
             y=paste("Megatrend Base 2000 +", int_label))
graphics.off()
