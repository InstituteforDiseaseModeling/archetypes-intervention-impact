
library(data.table)
library(ggplot2)

itn_func <- function(x, B){ (B^0.4 -B^20) * x^(1.5*B)}
interact<-function(x,y) (y^0.4-y^20)*(x^(1.5*y))

itn_adjust <- data.table(expand.grid(seq(0,1,0.1), seq(0,1,0.01)))
names(itn_adjust) <- c("x", "B")
itn_adjust$result <- mapply(interact, itn_adjust$x, itn_adjust$B)
itn_adjust[, x:=as.factor(x)]
setnames(itn_adjust, "x", "Coverage")

ggplot(itn_adjust, aes(x=B, y=result, color=Coverage)) +
        geom_line() + 
        labs(x="Baseline Prevalence (B)",
             y="Post-adjustment value") +
        theme_minimal()
        
itn_adjust_2 <- data.table(expand.grid(seq(0,1,0.1), seq(0,1,0.01)))
names(itn_adjust_2) <- c("B", "x")
itn_adjust_2$result <- mapply(itn_func, itn_adjust_2$x, itn_adjust_2$B)
itn_adjust_2[, B:=as.factor(B)]
setnames(itn_adjust_2, "B", "Prevalence")

pdf("Desktop/itn.pdf", width=6, height=5)
ggplot(itn_adjust_2, aes(x=x, y=result, color=Prevalence)) +
        geom_abline()+
        geom_line() + 
        labs(x="Net Coverage",
             y="Post-adjustment value") +
        theme_minimal()
graphics.off()
