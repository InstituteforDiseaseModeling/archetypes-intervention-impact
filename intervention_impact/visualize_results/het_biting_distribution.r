library(data.table)
library(ggplot2)
library(gridExtra)
library(MASS)

rm(list=ls())
theme_set(theme_minimal(base_size = 18))
in_dir <- file.path(Sys.getenv("USERPROFILE"), 
                      "Dropbox (IDM)/Malaria Team Folder/data/Burkina/INDIE/Fig2a_Guelbeogo_elife2018.csv")

data <- fread(in_dir)
data <- data[!is.na(Peak), list(id=V1, bites=Peak)]

# counts of bites
ggplot(data, aes(x=id, y=bites)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  labs(title="Peak Season Bites by Individual")

# convert to a risk (ie proportion of bites)
data[, risk:= bites/sum(bites)]
ggplot(data, aes(x=id, y=risk)) +
  geom_bar(stat="identity")+
  theme_minimal() +
  labs(title="Fraction of Bites by Individual",
       y="Risk")

# convert to a risk relative to the mean(?)
data[, relrisk:= risk/mean(risk)]
ggplot(data, aes(x=id, y=relrisk)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept=1, size=1.5, color="red") +
  theme_minimal() +
  labs(title="Biting Risk by Individual (Relative to Mean)",
       y="Relative Risk")

# fit distribution
histval <- data$relrisk
dist_fit <- fitdistr(histval, densfun="exponential")

h<-hist(histval,breaks=30)
histo <- data.table(xhist=c(min(h$breaks),h$breaks),
                    yhist=c(0,h$density,0))

xfit <- seq(0,max(histval),length=40)
dist <- data.table(xfit=xfit,
                   yfit=dexp(xfit, rate=dist_fit$estimate))
out_dir <- file.path(Sys.getenv("USERPROFILE"), 
          "Dropbox (IDM)/Malaria Team Folder/projects/map_intervention_impact/lookup_tables/interactions/gates_examples/")
pdf(file.path(out_dir, "expbiting.pdf"), width=5, height=5)
ggplot(histo, aes(x=xhist, y=yhist)) +
  geom_bar(stat="identity", fill="#AF81AD") + 
  geom_line(data=dist, aes(x=xfit, y=yfit), color="#F1B657", size=1.5) +
  labs(x="Relative Risk",
       y="")
graphics.off()

