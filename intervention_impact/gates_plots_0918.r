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

pdf(paste0(gates_dir, "/biting.pdf"))
ggplot(all_nets[ITN_Coverage %in% c(0.2, 0.6)], aes(x=mean_initial, y=mean_final)) +
  geom_abline(size=1.5)+
  # geom_line(aes(group=interaction(ITN_Coverage, Run_Number, net_type)), alpha=0.25) + todo: geom_ribbons
  geom_line(aes(group=interaction(net_type, ITN_Coverage), color=interaction(net_type, ITN_Coverage)), size=1.5) +
  # scale_color_manual(values=brewer.pal(7, "RdPu")[3:7], name="ITN Coverage") +
  scale_color_brewer(type="qual", palette = "Paired") + 
  theme_minimal() +
  labs(x="Initial PfPR",
       y="Final PfPR") 
graphics.off()


# Main takeaways: djibo, kananga, karen at 80% each
takeaway_sites <- c("djibo", "kananga", "karen")
takeaway_labels <- c("SE Asia", "Sahel", "Central Africa")

takeaways <- full_sweep_data[Site_Name %in% takeaway_sites]
takeaways[, site:=factor(Site_Name, levels=takeaway_sites, labels=takeaway_labels)]

# 1: Differing residual transmission, but same overall story
# png(paste0(gates_dir, "/resid_trans.png"), width=1300, height=600, res=170)
res_trans_ints <- c("ITN 0.6; " #, "IRS 0.6; ", "ITN 0.6; IRS 0.6; "
                    )
int_count <- length(res_trans_ints)
subset <- takeaways[Intervention %in% res_trans_ints]
subset[, Intervention:=factor(Intervention, levels=res_trans_ints)]
ggplot(subset[Site_Name!="karen"], aes(x=mean_initial, y=mean_final)) +
  geom_abline(size=1.5)+
  # geom_line(aes(group=interaction(Intervention, Run_Number)), alpha=0.25) +
  geom_line(aes(color=interaction(Intervention, site)), size=2) +
  scale_color_manual(values=c(brewer.pal(5, "Reds")[2:(2+int_count-1)], brewer.pal(5, "Blues")[2:(2+int_count-1)])) + 
  theme_minimal() +
  labs(x="Initial PfPR",
       y="Final PfPR") 
# graphics.off()

# 2: Different relative effects in different places 
int_subset <- c("ACT 0.4; ", "ITN 0.2; ACT 0.4; ", "IRS 0.4; ACT 0.6; ", "ITN 0.6; IRS 0.6; ")
subset <- takeaways[Intervention %in% int_subset]
subset[, Intervention:=factor(Intervention, levels=int_subset)]

# png(paste0(gates_dir, "/relative_imp_1.png"), width=1300, height=600, res=170)
ggplot(subset, aes(x=mean_initial, y=mean_final)) +
  geom_abline(size=1.5, alpha=0.1)+
  # geom_line(aes(group=interaction(Intervention, Run_Number)), alpha=0.25) +
  geom_line(aes(color=interaction(Intervention)), size=1.5) +
  scale_color_manual(values=c("#252525", "#969696", "#dd3497", "#fa9fb5")) + 
  theme_minimal() +
  labs(x="Initial PfPR",
       y="Final PfPR") +
  facet_grid(~site)
# graphics.off()

int_subset <- c("ITN 0.8; ", "IRS 0.8; ", "ACT 0.8; ", "ITN 0.8; IRS 0.8; ")
subset <- takeaways[Intervention %in% int_subset]
subset[, Intervention:=factor(Intervention, levels=int_subset)]

png(paste0(gates_dir, "/relative_imp_2.png"), width=1300, height=600, res=170)
ggplot(subset, aes(x=mean_initial, y=mean_final)) +
  geom_abline(size=1.5, alpha=0.25)+
  # geom_line(data=takeaways[Intervention=="ITN 0.8; IRS 0.8; "], linetype=4, color="#fa9fb5", size=1.5) +
  geom_line(size=1.5, aes(color=Intervention, linetype=Intervention)) +
  scale_color_manual(values=c("#252525", "#969696", "#dd3497", "#fa9fb5")) +
  scale_linetype_manual(values=c(3,2,1,4)) +
  theme_minimal() +
  labs(x="Initial PfPR",
       y="Final PfPR",
       title="Varying Best Intervention in Sites") +
  facet_grid(~site)
graphics.off()