
rm(list=ls())

library(data.table)
library(ggplot2)
library(rstan)


input_dir <- "/Volumes/GoogleDrive/My Drive/itn_cube/create_database/input"

# p0 & p1-- from stock & flow, nat'l time series of p0=p(hh has >0 nets) and p1=avg # of nets
# 40 countries (list length), houshold size 1-10 (columns)
load(file.path(input_dir, "prop0prop1.rData"))

format_dataset <- function(subset, metric, iso3){
  
  subset$time <- as.numeric(rownames(subset))
  subset <- data.table(subset)
  subset[, iso3:=iso3]
  subset <- melt(subset, id.vars=c("time", "iso3"), variable.name="hh_size", value.name=metric)
  subset[, hh_size:=as.integer(hh_size)]
  
  return(subset)
}

stock_and_flow <- lapply(1:length(out), function(idx){
  iso3 <- Cout[[idx]]
  
  subset_list <- out[[idx]]
  p0 <- format_dataset(data.frame(subset_list[,,1]), metric="prob_hasnet", iso3 = iso3)
  p1 <- format_dataset(data.frame(subset_list[,,2]), metric="mean_nets_per_hh", iso3 = iso3)
  country_subset <- merge(p0, p1, by=c("iso3", "hh_size", "time"), all=T)
  
  return(country_subset)
})

stock_and_flow <- rbindlist(stock_and_flow)
stock_and_flow[, hh_size:=as.factor(hh_size)]

ggplot(stock_and_flow[hh_size==2], aes(x=time, y=1-prob_hasnet)) +
  geom_line() + 
  facet_wrap(~iso3) +
  labs(x="Time",
       y="P0 (Probability HH has Any Nets)",
       title="Stock and Flow Means, HH Size=2")







