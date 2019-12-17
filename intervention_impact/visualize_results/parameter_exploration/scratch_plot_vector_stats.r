
library(data.table)
library(ggplot2)
library(bit64)
library(plyr)

in_dir <- "/Users/bertozzivill/Downloads/ReportVectorStats.csv"

data <- fread(in_dir)
data[, NodeID:=as.character(NodeID)]

from <- as.character(c(1475161276, 1518218607, 1402613965, 1585194826, 1671175813, 1758473462, 2186687066, 820652505))
to <- c("Nigeria", "Cameroon", "Burkina", "DRC", "Mozambique", "Ethiopia", "Myanmar", "Ecuador")

data[, site:= mapvalues(NodeID, from, to)]

ggplot(data, aes(x=Time, y=VectorPopulation)) + 
  geom_line(aes(color=Species)) +
  facet_wrap(~site)