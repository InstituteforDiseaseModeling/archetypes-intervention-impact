library(data.table)
library(ggplot2)
library(gridExtra)
library(MASS)

rm(list=ls())

in_dir <- file.path(Sys.getenv("USERPROFILE"), 
                      "Dropbox (IDM)/Malaria Team Folder/data/Burkina/INDIE/Fig2a_Guelbeogo_elife2018.csv")

data <- fread(in_dir)
data <- data[!is.na(Peak), list(id=V1, bites=Peak)]

# counts of bites
ggplot(data, aes(x=id, y=bites)) +
  geom_bar(stat="identity")

# convert to a risk (ie proportion of bites)
data[, risk:= bites/sum(bites)]
ggplot(data, aes(x=id, y=risk)) +
  geom_bar(stat="identity")

# convert to a risk relative to the mean(?)
data[, relrisk:= risk/mean(risk)]
# ggplot(data, aes(x=id, y=relrisk)) +
#   geom_bar(stat="identity")

# fit distribution
histval <- data$relrisk
dist_fit <- fitdistr(histval, densfun="exponential")

h<-hist(histval,breaks=30)
histo <- data.table(xhist=c(min(h$breaks),h$breaks),
                    yhist=c(0,h$density,0))

xfit <- seq(0,max(histval),length=40)
dist <- data.table(xfit=xfit,
                   yfit=dexp(xfit, rate=dist_fit$estimate))

ggplot(histo, aes(x=xhist, y=yhist)) +
  geom_bar(stat="identity") + 
  geom_line(data=dist, aes(x=xfit, y=yfit), color="red", size=1) +
  theme_minimal() +
  labs(x="Relative Risk",
       y="Density",
       title="Data and Exponential Distirbution, Biting Risk")


