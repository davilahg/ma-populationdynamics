Seeds <- read.csv("nseeds.csv")
nplot <- nlevels(as.factor(as.character(Survival$plot)))
k.p.list <- list()
mat.p <- 5 # mature plot
Survival <- subset(Survival, !is.na(ln.h2))
