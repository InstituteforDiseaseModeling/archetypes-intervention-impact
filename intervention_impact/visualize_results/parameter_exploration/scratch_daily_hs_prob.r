library(data.table)
library(ggplot2)

daily_prob<-0.7
not_prob<- 1-daily_prob
data <- data.table(day=1:20)
data[, prob:= (1-not_prob^(day))]

ggplot(data, aes(x=day, y=prob)) +
  geom_point() +
  geom_hline(yintercept = 0.90)