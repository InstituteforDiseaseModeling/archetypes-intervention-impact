library(data.table)
library(ggplot2)
library(gridExtra)
library(viridis)

rm(list=ls())

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
colors <- gg_color_hue(5)[2:5]

main_dir <- file.path(Sys.getenv("USERPROFILE"), 
                      "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/lookup_tables/interactions")

gates_dir <- file.path(main_dir, "gates_examples")
full_sweep_path <- file.path(main_dir, "version_2", "lookup_full_interactions_v2.csv")

get_smooth <- function(x, y){
  lo <- loess(y[y>0]~x[y>0])
  predictions <- c(y[y==0], predict(lo))
  return(pmax(predictions, rep(0, length(predictions))))
}


# setup: load all data, aggregate homo biting and corr net examples
full_sweep_data <- fread(full_sweep_path)

# homo biting
initial_homo <- fread(file.path(gates_dir, "initial_homo_biting.csv"))
final_homo <- fread(file.path(gates_dir, "itns_homo_biting.csv"))
homo_nets <- merge(initial_homo, final_homo, by=c("Site_Name", "Run_Number", "x_Temporary_Larval_Habitat"), all=T)
homo_nets[, Run_Number:=factor(Run_Number)]
homo_nets[, mean_initial:= mean(initial_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, ITN_Coverage)]
homo_nets[, mean_final:=mean(final_prev), by=list(Site_Name, x_Temporary_Larval_Habitat, ITN_Coverage)]
homo_nets[, net_type:="Homogeneous Biting"]
het_nets <- full_sweep_data[Site_Name=="kananga" & IRS_Coverage==0 & ACT_Coverage==0,
                            list(Site_Name, net_type="Heterogeneous Biting", Run_Number, x_Temporary_Larval_Habitat,
                                 ITN_Coverage, initial_prev, final_prev, mean_initial, mean_final)]
all_nets <- rbind(het_nets, homo_nets)
all_nets[, ITN_Coverage:=factor(ITN_Coverage)]
all_nets[, net_type:=factor(net_type, levels=c("Homogeneous Biting", "Heterogeneous Biting"))]

minmaxes <- all_nets[, list(mean_initial=unique(mean_initial),
                            min_final=min(final_prev),
                            max_final=max(final_prev)),
                     by=list(net_type, x_Temporary_Larval_Habitat, ITN_Coverage)]

minmaxes_smooth <- lapply(unique(minmaxes$net_type), function(net){
  sub_list <- lapply(unique(minmaxes$ITN_Coverage), function(cov){
    subset <- minmaxes[net_type==net & ITN_Coverage==cov]
    subset[, smooth_min:= get_smooth(mean_initial, min_final)]
    subset[, smooth_max:= get_smooth(mean_initial, max_final)]
  })
  sub_list <- rbindlist(sub_list)
})
minmaxes_smooth <- rbindlist(minmaxes_smooth)

summary_nets <- unique(all_nets[, list(net_type, ITN_Coverage, x_Temporary_Larval_Habitat, mean_initial, mean_final)])
summary_nets <- merge(summary_nets, 
                      minmaxes_smooth, 
                      by=c("net_type", "ITN_Coverage", "x_Temporary_Larval_Habitat", "mean_initial"))


# pdf(paste0(gates_dir, "/biting.pdf"))
ggplot(summary_nets[ITN_Coverage %in% c(0.2, 0.6)], aes(x=mean_initial, y=mean_final,
                                                    fill=interaction(net_type, ITN_Coverage),
                                                    group=interaction(net_type, ITN_Coverage))) +
  geom_abline(size=1.5)+
  # geom_ribbon(aes(ymin=min_final, ymax=max_final), alpha=0.75) + 
  geom_ribbon(aes(ymin=smooth_min, ymax=smooth_max), alpha=0.25) +
  geom_line(aes(color=interaction(net_type, ITN_Coverage)), size=1.5) +
  # scale_color_manual(values=brewer.pal(7, "RdPu")[3:7], name="ITN Coverage") +
  scale_color_brewer(type="qual", palette = "Paired") + 
  scale_fill_brewer(type="qual", palette = "Paired") + 
  theme_minimal() +
  labs(x="Initial PfPR",
       y="Final PfPR") 
# graphics.off()


# Main takeaways: djibo, kananga, karen at 80% each
takeaway_sites <- c("djibo", "kananga", "karen")
takeaway_labels <- c("SE Asia", "Sahel", "Central Africa")

takeaways <- full_sweep_data[Site_Name %in% takeaway_sites]
takeaways[, site:=factor(Site_Name, levels=takeaway_sites, labels=takeaway_labels)]


minmaxes <- takeaways[, list(mean_initial=unique(mean_initial),
                            min_final=min(final_prev),
                            max_final=max(final_prev)),
                     by=list(site, Intervention, x_Temporary_Larval_Habitat)]

minmaxes_smooth <- lapply(unique(minmaxes$site), function(site_name){
  sub_list <- lapply(unique(minmaxes$Intervention), function(int_name){
    subset <- minmaxes[site==site_name & Intervention==int_name]
    subset[, smooth_min:= get_smooth(mean_initial, min_final)]
    subset[, smooth_max:= get_smooth(mean_initial, max_final)]
  })
  sub_list <- rbindlist(sub_list)
})
minmaxes_smooth <- rbindlist(minmaxes_smooth)

summary_result <- unique(takeaways[, list(site, Intervention, x_Temporary_Larval_Habitat, mean_initial, mean_final)])
summary_result <- merge(summary_result, 
                      minmaxes_smooth, 
                      by=c("site", "Intervention", "x_Temporary_Larval_Habitat", "mean_initial"))

# 1: Differing residual transmission, but same overall story
# png(paste0(gates_dir, "/resid_trans.png"), width=1300, height=600, res=170)
res_trans_ints <- c("ITN 0.6; " , "IRS 0.6; ", "ITN 0.6; IRS 0.6; "
  )
int_count <- length(res_trans_ints)
subset <- summary_result[Intervention %in% res_trans_ints]
subset[, Intervention:=factor(Intervention, levels=res_trans_ints)]
ggplot(subset[site!="SE Asia"], aes(x=mean_initial, y=mean_final)) +
  geom_abline(size=1.5)+
  geom_ribbon(aes(ymin=smooth_min, ymax=smooth_max, fill=interaction(Intervention, site)), alpha=0.25) +
  geom_line(aes(color=interaction(Intervention, site)), size=2) +
  scale_color_manual(values=c(brewer.pal(5, "Reds")[2:(2+int_count-1)], brewer.pal(5, "Blues")[2:(2+int_count-1)])) + 
  scale_fill_manual(values=c(brewer.pal(5, "Reds")[2:(2+int_count-1)], brewer.pal(5, "Blues")[2:(2+int_count-1)])) + 
  theme_minimal() +
  labs(x="Initial PfPR",
       y="Final PfPR") 
# graphics.off()

# 2: Different relative effects in different places 
int_subset <- c("ACT 0.2; ", "ACT 0.6; ", "ITN 0.4; ACT 0.2; ", "IRS 0.4; ACT 0.2; ", "ITN 0.4; IRS 0.4; ACT 0.2; ")
subset <- summary_result[Intervention %in% int_subset]
subset[, Intervention:=factor(Intervention, levels=int_subset)]


# png(paste0(gates_dir, "/relative_imp_1.png"), width=1300, height=600, res=170)
ggplot(subset, aes(x=mean_initial, y=mean_final, color=Intervention)) +
  geom_abline(size=1.5, alpha=0.1)+
  geom_ribbon(aes(ymin=smooth_min, ymax=smooth_max, fill=Intervention), alpha=0.25) +
  geom_line(size=1.5) +
  scale_color_manual(values=brewer.pal(6, "PuRd")[2:6]) + 
  scale_fill_manual(values=brewer.pal(6, "PuRd")[2:6]) + 
  theme_minimal() +
  labs(x="Initial PfPR",
       y="Final PfPR") +
  facet_grid(~site)
# graphics.off()

