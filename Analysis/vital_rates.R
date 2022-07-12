# Population dynamics of Mimosa acantholoba
#
try(setwd("/Volumes/GoogleDrive/Mi\ unidad/Edgar/Trabajo/FC/Direcci�n\ de\ tesis/L\ 2017\ Gerardo\ D�vila"), silent = TRUE)
try(setwd("/Volumes/GoogleDrive/My\ Drive/Edgar/Trabajo/FC/Direcci�n\ de\ tesis/L\ 2017\ Gerardo\ D�vila"), silent = TRUE)
try(setwd('Insync/davilahg@ciencias.unam.mx/Google Drive/L 2017 Gerardo D�vila/'), silent = TRUE)
try(setwd('/home/gerardo/ma-populationdynamics'), silent = TRUE)
# load requiered libraries:
#install.packages(c('AICcmodavg', 'betareg', 'EnvStats', 'gamm4', 'ggplot2', 'mgcv', 'lme4', 'mgcv', 'reshape2', 'gridExtra', 'plot3D', 'fields'))
library(AICcmodavg)
library(betareg)
library(EnvStats)
library(gamm4)
library(ggplot2)
library(metR)
library(lme4)
library(mgcv)
library(reshape2)
library(gridExtra)
library(plot3D)
library(fields)
library(brms)
library(patchwork)
library(ggpubr)
library(dplyr)
library(MASS)
library(ggnewscale)

# read data
{
  data <- read.csv("./Data/cleaned_data/complete_data.csv")
  survival <- subset(data, !is.na(h1)) # s
  survival$sup <- as.factor(as.character(survival$sup))
  growth <- subset(data, !is.na(h1) & !is.na(h2) & sup == 1) # g
  data_r <- read.csv("./Data/cleaned_data/reproduction.csv")
  Age.mature <- 100
  reproduction <- transform(data_r, Age_init = 0, Fst_cens = 2003, Age = 0, ln.h1 = log(h1))
  for (i in 1:nrow(reproduction)) {
      if (reproduction$PlotName[i] == "DIV") {reproduction$Fst_cens[i] <- 2005}
  }
  for (i in 1:nrow(reproduction)) { # setting initial ages
      if (reproduction$PlotName[i] == "RIC") {reproduction$Age_init[i] <- 30}
      else if (reproduction$PlotName[i] == "ISP") {reproduction$Age_init[i] <- 23}
      else if (reproduction$PlotName[i] == "SEP") {reproduction$Age_init[i] <- 40}
      else if (reproduction$PlotName[i] == "ABE") {reproduction$Age_init[i] <- 10}
      else if (reproduction$PlotName[i] == "SEL") {reproduction$Age_init[i] <- Age.mature}
      else if (reproduction$PlotName[i] == "ISC") {reproduction$Age_init[i] <- 18}
      else if (reproduction$PlotName[i] == "BES") {reproduction$Age_init[i] <- 11}
      else if (reproduction$PlotName[i] == "TOA") {reproduction$Age_init[i] <- 0}
      else if (reproduction$PlotName[i] == "MAL") {reproduction$Age_init[i] <- 36}
      else if (reproduction$PlotName[i] == "HIL") {reproduction$Age_init[i] <- 5}
      else if (reproduction$PlotName[i] == "MAR") {reproduction$Age_init[i] <- 3}
      else if (reproduction$PlotName[i] == "FID") {reproduction$Age_init[i] <- 9}
      else if (reproduction$PlotName[i] == "ESS") {reproduction$Age_init[i] <- 16}
      else if (reproduction$PlotName[i] == "TOB") {reproduction$Age_init[i] <- 1}
      else if (reproduction$PlotName[i] == "DIA") {reproduction$Age_init[i] <- 7}
      else if (reproduction$PlotName[i] == "DIV") {reproduction$Age_init[i] <- 56}
  }
  for (i in 1:nrow(reproduction)) {
   reproduction$Age[i] <- reproduction$Census[i] - reproduction$Fst_cens[i] + reproduction$Age_init[i]
   if (reproduction$Age[i] > Age.mature) {
     reproduction$Age[i] <- Age.mature
   } # f1
  }
  fruitdata <- subset(reproduction, Rep > 0)
  fruitdata <- transform(fruitdata, TotalFruitset = TotalFruitset-1) # f2
  seeds <- read.csv('./Data/cleaned_data/seeds.csv') # f3
  seeds <- transform(seeds, ln.h1 = log(h1))
  obs_recruits <- read.csv('./Data/cleaned_data/observed_recruitment.csv')
  obs_recruits<- subset(obs_recruits, !is.na(Altura))
  obs_recruits.h <- obs_recruits$Altura
  understory <- read.csv('./Data/cleaned_data/understory.csv')
  understory <- subset(understory, !is.na(h1))
  understory.h <- understory$h1
  smalltrees.h <- c(obs_recruits.h, understory.h)
  smalltrees.h.l <- log(smalltrees.h)
  smalltrees.k <- kmeans(smalltrees.h.l, 6, iter.max = 15, nstart = 25)
}
# definitions
{
 m <- 100 # small-matrix dimension
 M <- m # big-matrix dimension
 min.lh1 <- min(survival$ln.h1, na.rm = TRUE)-0.0001
 min.lh2 <- min(survival$ln.h1, na.rm = TRUE)-0.0001
 max.lh1 <- max(survival$ln.h1, na.rm = TRUE)+0.0001
 max.lh2 <- max(survival$ln.h2, na.rm = TRUE)+0.0001
 E.pred <- seq(min(min.lh1, min.lh2), max(max.lh1, max.lh2), length.out = M+1)
 h.pred <- seq(min(survival$h1, na.rm = TRUE), max(survival$h1, na.rm = TRUE), length.out = m)
 h <- E.pred[2]-E.pred[1]
 X.pred <- (E.pred[2:(M+1)]+E.pred[1:M])/2
 t1.pred <- rep(0, length.out = M)
 sM <- M/m
 e.pred <- seq(min(X.pred), max(X.pred), length.out = m+1)
 x.pred <- (e.pred[2:(m+1)]+e.pred[1:m])/2
 Age.pred <- 1:Age.mature
 nplot <- nlevels(survival$PlotCode)
 ex.pr <- exp(x.pred)
 total_area <- 400 # m^2 per plot
 understory_area <- 32 # m^2 per plot
 total_plot_n <- 17
 understory_plot_n <- 5
 area_scale <- total_area / understory_area
}
# modelling vital rates
{
 #gam.s <- gamm4(sup~t2(ln.h1, Age, k = c(4,4)), random = ~(1+ln.h1+Age|Census)+(0++Age+ln.h1|PlotCode), data = survival, family = binomial)
 gam.s <- readRDS('./Data/r_data/survival_model.RDS')
 #gam.g <- gamm4(ln.h2~t2(ln.h1, Age, k = c(4,4)), random = ~(1+Age+ln.h1+Age:ln.h1|Census)+(1+Age+ln.h1+Age:ln.h1|PlotCode/UID), data = growth)
 gam.g <- readRDS('./Data/r_data/growth_model.RDS')
 gam.f1 <- gamm4(Rep~s(ln.h1, k = 4)+s(Age, k = 4), random = ~(0+Age+ln.h1|PlotCode), family = binomial, data = reproduction)
 gam.f2 <- readRDS('./Data/r_data/gam-f2.rds')
 gam.f3 <- readRDS('./Data/r_data/gam-f3.rds')
} 
# f4 data processing
{
  und_plots <- c("SC02", "SC04", "SC06", "SC07", "SC08", "SC10", "SC11", "SC12", "SC13", "SC14", "SC15", "SC16", "SC17", "2", "4", "6", "7", "8", "10", "11", "12", "13", "14", "15", "16", "17")
  und_surv_rows <- which(survival$PlotCode %in% und_plots) # get data from the main dataset which was collected in the same plots as understory data
  und_survival <- survival[und_surv_rows,] # 13 PARCELAS
  und_age_range <- sort(unique(und_survival$Age))  
  tot_seeds_t <- as.data.frame(matrix(NA, nrow = max(und_age_range), ncol = 6))
  names(tot_seeds_t) <- c('Age', 'SeedsUnd', 'Seeds1', 'Seeds2.5', 'Seeds5','Seeds')
  tot_seeds_t$Age <- 1:max(und_age_range)
  und_s_und <- und_survival[which(und_survival$maxDBH < 1),]
  und_s_1 <- und_survival[which(und_survival$maxDBH >= 1 & und_survival$maxDBH < 2.5),]
  und_s_2.5 <- und_survival[which(und_survival$maxDBH >= 2.5 & und_survival$maxDBH < 5),]
  und_s_5 <- und_survival[which(und_survival$maxDBH >= 5),]
  for (a in und_age_range) {
    n_t_rows <- which(und_s_und$Age == a) ## DB < 1
    n_t_data <- und_s_und[n_t_rows, ]
    n_t_data <- subset(n_t_data, sup == 1)
    if (nrow(n_t_data) > 0) {
      n_f1 <- predict(gam.f1$gam, newdata = data.frame(ln.h1 = n_t_data$ln.h1, Age = rep(a, nrow(n_t_data))), type = "response", re.form = NA)
      n_f2 <- fitted(gam.f2, newdata = data.frame(ln.h1 = n_t_data$ln.h1, Age = rep(a, nrow(n_t_data))), type = "response", re_formula = NA)[, 1]
      n_f3 <- fitted(gam.f3, newdata = data.frame(ln.h1 = n_t_data$ln.h1, Age = rep(a, nrow(n_t_data))), type = "response", re_formula = NA)[, 1]
      n_f123_t <- sum(n_f1*n_f2*n_f3)*area_scale # (12.5) 32m2 to 400m2
      tot_seeds_t$SeedsUnd[a] <- n_f123_t
    }
  }
  for (a in und_age_range) {
    n_t_rows <- which(und_s_1$Age == a) ## 1 <= DB < 2.5
    n_t_data <- und_s_1[n_t_rows, ]
    n_t_data <- subset(n_t_data, sup == 1)
    if (nrow(n_t_data) > 0) {
      n_f1 <- predict(gam.f1$gam, newdata = data.frame(ln.h1 = n_t_data$ln.h1, Age = rep(a, nrow(n_t_data))), type = "response", re.form = NA)
      n_f2 <- fitted(gam.f2, newdata = data.frame(ln.h1 = n_t_data$ln.h1, Age = rep(a, nrow(n_t_data))), type = "response", re_formula = NA)[, 1]
      n_f3 <- fitted(gam.f3, newdata = data.frame(ln.h1 = n_t_data$ln.h1, Age = rep(a, nrow(n_t_data))), type = "response", re_formula = NA)[, 1]
      n_f123_t <- sum(n_f1*n_f2*n_f3)*4 # 100m2 to 400m2
      tot_seeds_t$Seeds1[a] <- n_f123_t
    }
  }
  for (a in und_age_range) {
  	n_t_rows <- which(und_s_2.5$Age == a) ## 2.5 <= DB < 5
  	n_t_data <- und_s_2.5[n_t_rows, ]
  	n_t_data <- subset(n_t_data, sup == 1)
  	if (nrow(n_t_data) > 0) {
    	n_f1 <- predict(gam.f1$gam, newdata = data.frame(ln.h1 = n_t_data$ln.h1, Age = rep(a, nrow(n_t_data))), type = "response", re.form = NA)
    	n_f2 <- fitted(gam.f2, newdata = data.frame(ln.h1 = n_t_data$ln.h1, Age = rep(a, nrow(n_t_data))), type = "response", re_formula = NA)[, 1]
    	n_f3 <- fitted(gam.f3, newdata = data.frame(ln.h1 = n_t_data$ln.h1, Age = rep(a, nrow(n_t_data))), type = "response", re_formula = NA)[, 1]
    	n_f123_t <- sum(n_f1*n_f2*n_f3)*2 # 200m2 to 400m2
    	tot_seeds_t$Seeds2.5[a] <- n_f123_t
  	  }
    }
  for (a in und_age_range) {
	n_t_rows <- which(und_s_5$Age == a) ## DB >= 5
	n_t_data <- und_s_5[n_t_rows, ]
	n_t_data <- subset(n_t_data, sup == 1)
	if (nrow(n_t_data) > 0) {
  	n_f1 <- predict(gam.f1$gam, newdata = data.frame(ln.h1 = n_t_data$ln.h1, Age = rep(a, nrow(n_t_data))), type = "response", re.form = NA)
  	n_f2 <- fitted(gam.f2, newdata = data.frame(ln.h1 = n_t_data$ln.h1, Age = rep(a, nrow(n_t_data))), type = "response", re_formula = NA)[, 1]
  	n_f3 <- fitted(gam.f3, newdata = data.frame(ln.h1 = n_t_data$ln.h1, Age = rep(a, nrow(n_t_data))), type = "response", re_formula = NA)[, 1]
  	n_f123_t <- sum(n_f1*n_f2*n_f3) # 400m2 to 400m2
  	tot_seeds_t$Seeds5[a] <- n_f123_t
    }
  }
 tot_seeds_t$Seeds <- rowSums(tot_seeds_t[,c('SeedsUnd', 'Seeds1', 'Seeds2.5', 'Seeds5')], na.rm = TRUE)
 tot_seeds_t <- tot_seeds_t[und_age_range,]
  
 smalltrees_km <- kmeans(smalltrees.h, 6, nstart = 25, iter.max = 1000)
 min_cluster <- which(smalltrees_km$centers == min(smalltrees_km$centers))
 recruits_rows <- which(smalltrees_km$cluster == min_cluster)
 recruits_h <- smalltrees.h[recruits_rows]
 recruits_range <- range(recruits_h)
 recruits_t1 <- c()
 und_age_range_r <- und_age_range+1
 und_ages <- which(und_age_range_r %in% und_age_range)
 und_age_range1 <- und_age_range_r[und_ages]
 understory_raw <- read.csv("./Data/raw_data/raw_understory.csv") # lo necesito porque hay algunas observaciones que son sólo de un año, se desechan para supervivencia y crecimiento, pero sirve para ver el número de reclutas
 for (a in und_age_range) {
   if (a %in% und_age_range1) {
     n_t_rows <- which(understory_raw$Ecenso == a) # Ecenso = Age
     n_t_data <- understory_raw[n_t_rows, ]
     n_t_recruits <- subset(n_t_data, Altura <= recruits_range[2]) # Altura = h1
     n_t_size <- nrow(n_t_recruits)
     scaled_n_t_size <- n_t_size*area_scale # Solo se necesita reescalar por el �rea de la parcela porque en el sotobosque no hay divisi�n por di�metro basal y porque estamos usando las mismas parcelas para calcular f123
     recruits_t1 <- c(recruits_t1, scaled_n_t_size)    
   } else {
    recruits_t1 <- c(recruits_t1, NA)
   }
 }
 ####
 undn_plots <- c()
 undn_n <- c()
 for (p in unique(understory_raw$Plot)) {
    np <- nrow(subset(understory_raw, Plot == p & Supv == 1))
    undn_plots <- c(undn_plots, p)
    undn_n <- c(undn_n, np)
 }
 und_n_dist <- data.frame(plot = undn_plots, n = undn_n)
 sampled_plots_n <- length(unique(understory_raw[which(understory_raw$Ecenso %in% und_age_range),]$Plot)) # 13 PARCELAS DE DONDE SE SACÓ EL NÚMERO DE BRINZALES
 ####
 tot_seeds_t1 <- transform(tot_seeds_t, Recruits_t1 = recruits_t1, Establishment = recruits_t1/Seeds+0.0000001)
 tot_seeds_t1 <- subset(tot_seeds_t1, Establishment <= 1)
 p <- tot_seeds_t1$Establishment
# gam.f4.gam.br <- gam(Establishment ~ s(Age, k = 5), family = betar(link = 'logit'), data = tot_seeds_t1)
 f4.spline <- smooth.spline(tot_seeds_t1$Age, log(p/(1-p))) # f4
 f4.spline.pred <- predict(f4.spline, x = Age.pred)$y
# gam.f4.gam.br.pred <- predict(gam.f4.gam.br, newdata = data.frame(Age = Age.pred), type = 'response')
 p.logistic <- 1/(1+exp(-f4.spline.pred))
 tot_seeds_line <- data.frame(Age = Age.pred, Establishment = p.logistic)
 plot(tot_seeds_t1$Age, tot_seeds_t1$Establishment, xlim = c(1,100), ylim = c(0,1), xlab = 'Succesional age (years)', ylab = 'Establishment probability')
 lines(p.logistic, col = "blue")
# lines(gam.f4.gam.br.pred, col = 'red')
}
# sprouting and recruitment
{
understory_raw <- read.csv("./Data/raw_data/raw_understory.csv") 
understory_raw_corrected <- subset(understory_raw, DB.max < 2)
understory_raw_corrected <- droplevels(understory_raw_corrected)
und_raw_age_range <- sort(unique(understory_raw_corrected$Ecenso))
sprouts_h <- c()
sprouts_a <- c()
sprouts_ai <- c()
sprouts_dbmax <- c()
sprouts_dbmin <- c()
init_ages <- c(1, 9, 23, 62)
understory_raw_corrected$Etiqueta <- as.factor(as.character(understory_raw_corrected$Etiqueta))
for (a in und_raw_age_range) {
  und_a <- droplevels(subset(understory_raw_corrected, Ecenso == a))
  for (i in levels(und_a$Etiqueta)) {
    und_a_p <- subset(und_a, Etiqueta == i)
    if (nrow(und_a_p) > 0)
      if(!any(und_a_p$Age[1] == init_ages) & !is.na(und_a_p$DB3)[1]) {      
        sprouts_h <- c(sprouts_h, und_a_p$Altura[1])
        sprouts_a <- c(sprouts_a, und_a_p$Ecenso[1])
        sprouts_ai <- c(sprouts_ai, und_a_p$Einicial[1])
        sprouts_dbmin <- c(sprouts_dbmin, und_a_p$DB.min[1])
        sprouts_dbmax <- c(sprouts_dbmax, und_a_p$DB.max[1])

      }
  }
}
sprouting <- data.frame(Age = sprouts_a, Age_init = sprouts_ai, h1 = sprouts_h, DBmin = sprouts_dbmin, DBmax = sprouts_dbmax)
sprouts_a <- c()
for (a in sort(unique(sprouting$Age))) {
  sprouting_a <- subset(sprouting, Age == a)
  sprout_n <- nrow(sprouting_a)*area_scale
  sprouts_a <- c(sprouts_a, sprout_n)
}
nsprouts <- data.frame(Age = sort(unique(sprouting$Age)), n = sprouts_a)
##
tot.p.ages <- c()
for (p in unique(understory_raw$Plot)) {
  und.p <- subset(understory_raw, Plot == p & Supv != 0)
  p.ages <- unique(und.p$Ecenso)
  tot.p.ages <- c(tot.p.ages, p.ages)
}
plotsperage <- as.data.frame(table(tot.p.ages))
spr.ppa <- plotsperage[which(plotsperage$tot.p.ages %in% nsprouts$Age),]$Freq
nsprouts <- transform(nsprouts, ppa = spr.ppa)
nsprouts <- transform(nsprouts, nscaled = n/ppa)
f6.spline <- smooth.spline(nsprouts$Age, log(nsprouts$nscaled), spar = 0) # f6
f6.pred <- exp(predict(f6.spline, x = Age.pred)$y)
f6.line.df <- data.frame(Age = Age.pred, Nsprouts = f6.pred)
# density f5 & f7
den_recruits <- approxfun(density(x = log(recruits_h), n = m, na.rm = TRUE, from = e.pred[1], to = e.pred[m+1]))
den_recruits_xpred <- den_recruits(x.pred)
f5.den <- den_recruits_xpred/sum(den_recruits_xpred) # f5
#plot(exp(x.pred), density_recruits, col = 'blue', type = 'l')
den_sprouts <- approxfun(density(x = log(sprouting$h1), n = m, na.rm = TRUE, from = e.pred[1], to = e.pred[m+1]))
den_sprouts_xpred <- den_sprouts(x.pred)
f7.den <- den_sprouts_xpred/sum(den_sprouts_xpred) # f7
#lines(exp(x.pred), density_sprouts, col = 'red', type = 'l')
f57.df <- data.frame(height = c(recruits_h, sprouting$h1), group = as.factor(c(rep("recruits", length(recruits_h)), rep("sprouts", length(sprouting$h1)))))
}
# discretizing vital rates functions
{
   # g array
   GGG <- matrix(NA, ncol = Age.mature, nrow = M)
   for (a in 1:Age.mature)
    GGG[, a] <- predict(gam.g$gam, newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), type = "response", re.form = NA)
   MUG <- matrix(NA, ncol = Age.mature, nrow = m)
   for (a in 1:Age.mature) {
    init.g <- 1
    for (i in 1:m) {
       MUG[i, a] <- mean(GGG[init.g:(init.g+sM-1), a])
     init.g <- init.g + sM
    }
   }
   sdg <- sigma(gam.g$mer)
   G <- array(NA, dim = c(Age.mature, m, m))
   for (a in 1:Age.mature) {
    for (i in 1:m) {
     G[a, 1, i] <- pnorm(e.pred[2], MUG[i, a], sdg)
     for (j in 2:(m-1))
      G[a, j, i] <- pnorm(e.pred[j+1], MUG[i, a], sdg) - pnorm(e.pred[j], MUG[i, a], sdg)
     G[a, m, i] <- 1-pnorm(e.pred[m], MUG[i, a], sdg) # avoiding probability eviction
    }
   }
   # s matrix
   SSS <- matrix(NA, ncol = Age.mature, nrow = M)
   for (a in 1:Age.mature)
    SSS[ ,a] <- predict(gam.s$gam, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), re.form = NA)
   MUS <- matrix(NA, ncol = Age.mature, nrow = m)
   for (a in 1:(Age.mature)) {
    init.s <- 1
    for (i in 1:m) {
       MUS[i, a] <- mean(SSS[init.s:(init.s+sM-1), a])
       init.s <- init.s + sM
    }
   }
   S <- MUS
   # f1 matrix
   FF1 <- matrix (NA, nrow = M, ncol = Age.mature)
   for (a in 1:Age.mature)
    FF1[ ,a] <- predict(gam.f1$gam, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), re.form = NA)
   MUF1 <- matrix(NA, ncol = Age.mature, nrow = m)
   for (a in 1:Age.mature) {
      init.f1 <- 1
      for (i in 1:m) {
       MUF1[i, a] <- mean(FF1[init.f1:(init.f1+sM-1), a])
       init.f1 <- init.f1 + sM
    }
   }
   F1 <- MUF1
   # f2 matrix
   FF2 <- matrix (NA, nrow = M, ncol = Age.mature)
   for (a in 1:Age.mature)
    FF2[ ,a] <- fitted(gam.f2, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), re_formula = NA)[, 1]
   MUF2 <- matrix(NA, ncol = Age.mature, nrow = m)
   for (a in 1:Age.mature) {
      init.f2 <- 1
      for (i in 1:m) {
       MUF2[i, a] <- mean(FF2[init.f2:(init.f2+sM-1), a])
       init.f2 <- init.f2 + sM
    }
   }
   F2 <- MUF2
   # f3 matrix
   FF3 <- matrix (NA, nrow = M, ncol = Age.mature)
   for (a in 1:Age.mature)
    FF3[ ,a] <- fitted(gam.f3, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), re_formula = NA)[, 1]
   MUF3 <- matrix(NA, ncol = Age.mature, nrow = m)
   for (a in 1:Age.mature) {
      init.f3 <- 1
      for (i in 1:m) {
       MUF3[i, a] <- mean(FF3[init.f3:(init.f3+sM-1), a])
       init.f3 <- init.f3 + sM
    }
   }
   F3 <- MUF3
   # f4 vector
   F4 <- p.logistic
   # f5 vector
   F5 <- f5.den
}
# creating kernel
{
  k.i.j.a <- s.i.j.a <- g.i.j.a <- p.i.j.a <- f.i.j.a <- F5.i.j.a <- array(NA, dim = c(Age.mature, m, m))
  s.i.a <- S
  f.i.a <- matrix(NA, nrow = m, ncol = Age.mature)
  for (a in 1:Age.mature) {
    for (i in 1:m) {
     f.i.a[i, a] <- F1[i, a]*F2[i, a]*F3[i, a]*F4[a]
     for (j in 1:m) {
      p.i.j.a[a, j, i] <- s.i.a[i, a]*G[a, j, i]
      f.i.j.a[a, j, i] <- f.i.a[i, a]*F5[j]
      k.i.j.a[a, j, i] <- p.i.j.a[a, j, i]+f.i.j.a[a, j, i]
     }
    }
   }
}
# population growth rate
{
  understory_raw <- read.csv("./Data/raw_data/raw_understory.csv")
  Age1_n_rows <- which(understory_raw$Ecenso == 1)
  Age1_data <- understory_raw[Age1_n_rows,]
  Age1_data <- subset(Age1_data, Altura <= recruits_range[2])
  Age1_n <- nrow(Age1_data)*area_scale/length(unique(Age1_data$Plot)) # average number of recluits per plot at age 1
  str.size <- matrix(NA, m, Age.mature+1) ## para almancenar estructura de tamaños
  lam.list <- c() ## para almacenar lambdas
  init.n.a.v <- Age1_n*F5+f6.pred[1]*f7.den ## se agregan reclutas aquí
  str.size[ ,1] <- init.n.a.v/sum(init.n.a.v)
  n.list <- c(sum(init.n.a.v))
  n.pre <- c() ## número de arbolitos preexistentes
  n.spr <- c() ## número de arbolitos nuevos de rebrote
  for (a in 1:Age.mature) {
    n.a.v <- k.i.j.a[a,,]%*%init.n.a.v
    n.s.a <- f6.pred[a+1]*f7.den ## se agregan reclutas aquí
    lam.a <- (sum(n.a.v)+sum(n.s.a))/sum(init.n.a.v)
    str.size[, a+1] <- (n.a.v+n.s.a)/sum(n.a.v+n.s.a)
    n.pre <- c(n.pre, sum(n.a.v))
    n.spr <- c(n.spr, sum(n.s.a))
    n.list <- c(n.list, (sum(n.a.v)+sum(n.s.a)))
    lam.list <- c(lam.list, lam.a)
    init.n.a.v <- n.a.v+n.s.a
    }
  #tot.lam.pred <- lam.list.NM <- lam.list # lam.list is the transitory lambda vector
  pop_growth_rate.df <- data.frame(N = n.list, Age = c(0, Age.pred), lambda = c(NA, lam.list))
  pop_growth_rate.df <- subset(pop_growth_rate.df, N >= 0.9)
  max.age <- 38
  #pdf('lambda.pdf')
  plot(Age.pred, log(lam.list), type = 'l')
  abline(0,0)
}
# Size structure change
{
	str.size.ggdf <- as.data.frame(matrix(ncol = 3, nrow = m*Age.mature))
	names(str.size.ggdf) <- c('Age', 'Size', 'Density')
	for (a in 1:Age.mature) {
        str.size.ggdf$Age[((a-1)*Age.mature+1):(Age.mature*a)] <- rep(a, m)
        str.size.ggdf$Size[((a-1)*Age.mature+1):(Age.mature*a)] <- h.pred
        str.size.ggdf$Density[((a-1)*Age.mature+1):(Age.mature*a)] <- str.size[,a]
	}
	str.size.df <- subset(str.size.ggdf, Age <= max.age)
	str.size.plot <- ggplot(str.size.df, aes(x = Age, y = Size)) +
	                  theme_minimal()+
	                  coord_cartesian(xlim = c(1, max.age), ylim = c(min(h.pred), max(h.pred))) +
	                  geom_raster(aes(fill = Density)) +
	                  scale_fill_gradientn(colours = hcl.colors(12, "YlOrRd"))+
	                  geom_contour(aes(Age, Size, z = Density), color = 'black') +
	                  geom_label_contour(aes(Age, Size, z = Density)) +
	                  labs(x = 'Successional age (years)', y = 'Height (m)')
	  
	str.size.plot
}
# lambda by plot
{
  canopy_data <- read.csv('./Data/cleaned_data/raw_canopy_with_age.csv')
  canopy_data <- transform(canopy_data, h1 = Height, sup = as.factor(as.character(Sup)))
  canopy_data <- subset(canopy_data, !is.na(h1) & sup == 1 & PlotCode != '') 
  understory_data <- read.csv("./Data/raw_data/raw_understory.csv")
  understory_data <- transform(understory_data, h1 = Altura, sup = as.factor(as.character(Supv)), PlotCode = rep(NA, nrow(understory_data)), Age = Ecenso)
  canopy_data <- droplevels(canopy_data)
  understory_data <- droplevels(understory_data)  
  understory_data <- subset(understory_data, !is.na(h1) & sup == 1)
  # unifying plot coding between canopy and understory data
  understory_data$PlotCode <- as.character(understory_data$Plot)
  for (c in 1:length(understory_data$PlotCode)) {
    if (nchar(understory_data$PlotCode[c]) == 1) {
      understory_data$PlotCode[c] <- paste0("SC0", understory_data$PlotCode[c])
    } else if (nchar(understory_data$PlotCode[c]) == 2) {
      understory_data$PlotCode[c] <- paste0("SC", understory_data$PlotCode[c])
    }
  }
  lam.und <- lam.can <- as.data.frame(matrix(NA, 0, 4))
  names(lam.und) <- names(lam.can) <- c('PlotCode', 'Age', 'n', 'lambda')
  canopy_data$PlotCode <- as.factor(canopy_data$PlotCode)
  understory_data$PlotCode <- as.factor(understory_data$PlotCode)
  for (p in levels(canopy_data$PlotCode)) { # getting observed CANOPY N and lambda by plot by year
    can.p <- subset(canopy_data, PlotCode == p)
    can.p.a.n.list <- c()
    for (a in sort(unique(can.p$Age))) {
      can.p.a <- subset(can.p, Age == a)
      can.p.can <- length(which(can.p.a$maxDBH < 1))*area_scale
	  can.p.1 <- length(which(can.p.a$maxDBH >= 1 & can.p.a$maxDBH < 2.5))*4
	  can.p.2.5 <- length(which(can.p.a$maxDBH >= 2.5 & can.p.a$maxDBH < 5))*2
	  can.p.5 <- length(which(can.p.a$maxDBH >= 5))
	  can.p.a.n <- can.p.can + can.p.1 + can.p.2.5 + can.p.5 
      can.p.a.n.list <- c(can.p.a.n.list, can.p.a.n)
      }
    can.p.a.l.list <- c(NA, can.p.a.n.list[2:length(can.p.a.n.list)]/can.p.a.n.list[1:(length(can.p.a.n.list)-1)])
    can.p.a.df <- data.frame(PlotCode = rep(can.p$PlotCode[1], length(can.p.a.l.list)), Age = sort(unique(can.p$Age)), n = can.p.a.n.list, lambda = can.p.a.l.list)
    lam.can <- rbind(lam.can, can.p.a.df)
  }  
  #
  for (p in levels(understory_data$PlotCode)) { # getting observed UNDERSTORY N and lambda by plot by year
    und.p <- subset(understory_data, PlotCode == p)
    und.p.a.n.list <- c()
    for (a in sort(unique(und.p$Age))) {
      und.p.a <- subset(und.p, Age == a)
      und.p.a.n <- nrow(und.p.a)*area_scale
      und.p.a.n.list <- c(und.p.a.n.list, und.p.a.n)
      }
    und.p.a.l.list <- c(NA, und.p.a.n.list[2:length(und.p.a.n.list)]/und.p.a.n.list[1:(length(und.p.a.n.list)-1)])
    und.p.a.df <- data.frame(PlotCode = rep(und.p$PlotCode[1], length(und.p.a.l.list)), Age = sort(unique(und.p$Age)), n = und.p.a.n.list, lambda = und.p.a.l.list)
    lam.und <- rbind(lam.und, und.p.a.df)
  }
  tot_lam_rows <- which(levels(lam.can$PlotCode) %in% levels(lam.und$PlotCode)) # which plots have data both from understory and from canopy (lam.tot)
  tot_lam_plt <- levels(lam.can$PlotCode)[tot_lam_rows]
  lam.tot <- as.data.frame(matrix(NA, 0, 4))
  names(lam.tot) <- c('PlotCode', 'Age', 'n', 'lambda') 
  for (p in tot_lam_plt) { # getting observed both CANOPY and UNDERSTORY N and lambda by plot by year 
    tot.p.a.n.list <- c()
    tot.c.p <- subset(lam.can, PlotCode == p)
    tot.u.p <- subset(lam.und, PlotCode == p)
    tot.p.rows <- which(tot.c.p$Age %in% tot.u.p$Age)
    tot.p.a <- tot.c.p$Age[tot.p.rows]
    for (a in tot.p.a) {
      n <- tot.c.p$n[which(tot.c.p$Age == a)] + tot.u.p$n[which(tot.u.p$Age == a)]
      tot.p.a.n.list <- c(tot.p.a.n.list, n)
    }
    tot.p.a.l.list <- c(NA, tot.p.a.n.list[2:length(tot.p.a.n.list)]/tot.p.a.n.list[1:(length(tot.p.a.n.list)-1)])
    tot.p.a.df <- data.frame(PlotCode = rep(p, length(tot.p.a.n.list)), Age = tot.p.a, n = tot.p.a.n.list, lambda = tot.p.a.l.list)
    lam.tot <- rbind(lam.tot, tot.p.a.df)
  }
  # mantener lambda del dosel solo cuando no haya informaci�n completa
  only.can.lam.can <- data.frame(matrix(ncol = 4, nrow = 0))
  x <- c('PlotCode', 'Age', 'n', 'lambda')
  colnames(only.can.lam.can) <- x
  ages.tot <- unique(lam.tot$Age)
  lam.can$PlotCode <- as.factor(lam.can$PlotCode)
  lam.tot$PlotCode <- as.factor(lam.tot$PlotCode)
  for (p in unique(lam.can$PlotCode)) {
    lam.can.p <- subset(lam.can, PlotCode == p)
    lam.tot.p <- subset(lam.tot, PlotCode == p)
    only.can.ages <- which(!(unique(lam.can.p$Age) %in% unique(lam.tot.p$Age)))
    lam.can.p <- lam.can.p[only.can.ages,]
    only.can.lam.can <- rbind(only.can.lam.can, lam.can.p)
    }
  #
  Age.list <- c(pop_growth_rate.df$Age,lam.und$Age, lam.can$Age)
  Lam.list <- c(pop_growth_rate.df$lambda, lam.und$lambda, lam.can$lambda)
  Src.list <- c(rep('Model', length(pop_growth_rate.df$Age)), rep('Understory', length(lam.und$Age)), rep('Canopy', length(lam.can$Age)))
  lam.plot.df <- data.frame(Age = Age.list, lambda = Lam.list, source = Src.list)
  lam.plot.df <- subset(lam.plot.df, source != "Model") ##
  pop_growth_rate.df <- data.frame(N = n.list, Age = c(0, Age.pred), lambda = c(NA, lam.list)) ## para comparar todo
  # set maximum age at N < 1
  pop_growth_rate.df.g1 <- subset(pop_growth_rate.df, Age <= max.age)
  pop_growth_rate.df.l1 <- subset(pop_growth_rate.df, Age > max.age)
  #
  lam.p.plot <- ggplot(data = na.omit(pop_growth_rate.df.g1))+##, #aes(x = Age, y = lambda, group = source)) +
                theme_classic() +
                coord_cartesian(xlim = c(0,m)) +
                scale_y_continuous(trans = "log10") +
                labs(x = 'Succesional age (years)', y = 'Population growth rate') +
                ##geom_line(data = na.omit(only.can.lam.can), aes(x = Age, y = lambda, colour = PlotCode), linetype = 'dashed') +
                geom_line(data = na.omit(only.can.lam.can), aes(x = Age, y = lambda, group = PlotCode), color = 'blue', linetype = 'dotted') +
                ##geom_point(data = na.omit(only.can.lam.can), aes(x = Age, y = lambda, colour = PlotCode), shape = 21) +
                ##geom_point(data = na.omit(only.can.lam.can), aes(x = Age, y = lambda), shape = 21, color = "blue") +
                ##scale_color_manual(name = "Source", values = c("Model" = "red", "Observed\n(understory)" = "forestgreen", "Observed\n(canopy)" = "royalblue1"))
                ##geom_line(data = na.omit(lam.tot), aes(x = Age, y = lambda, colour = PlotCode), linetype = 'solid', size = 1) +
                geom_line(data = na.omit(lam.tot), aes(x = Age, y = lambda, group = PlotCode), color = 'blue', linetype = 'solid') +
                ##geom_point(data = na.omit(lam.tot), aes(x = Age, y = lambda, colour = PlotCode), shape = 16) +
                ##geom_point(data = na.omit(lam.tot), aes(x = Age, y = lambda), shape = 16, color = "blue") +
                geom_point(aes(x=100, y=1), color = 'blue', shape = 21) +
                geom_line(data = na.omit(pop_growth_rate.df.l1), aes(x = Age, y = lambda), size = 1, linetype = 'dotted') +
                geom_abline(intercept = 0, slope = 0) +
                geom_line(aes(x = Age, y = lambda), size = 1.2)

  lam.p.plot
  #
  N.p.plot <- ggplot(data = pop_growth_rate.df.g1)+##, #aes(x = Age, y = n, group = source)) +
                theme_classic() +
                coord_cartesian(xlim = c(0,m)) +
                scale_y_continuous(trans = "log10") +
                labs(x = 'Succesional age (years)', y = 'Population size per plot') +
                ##geom_line(data = na.omit(only.can.lam.can), aes(x = Age, y = n, colour = PlotCode), linetype = 'dashed') +
                geom_line(data = only.can.lam.can, aes(x = Age, y = n, group = PlotCode), color = 'blue', linetype = 'dotted') +
                ##geom_point(data = na.omit(only.can.lam.can), aes(x = Age, y = n, colour = PlotCode), shape = 21) +
                ##geom_point(data = na.omit(only.can.lam.can), aes(x = Age, y = n), shape = 21, color = "blue") +
                ##scale_color_manual(name = "Source", values = c("Model" = "red", "Observed\n(understory)" = "forestgreen", "Observed\n(canopy)" = "royalblue1"))
                ##geom_line(data = na.omit(lam.tot), aes(x = Age, y = n, colour = PlotCode), linetype = 'solid', size = 1) +
                geom_line(data = lam.tot, aes(x = Age, y = n, group = PlotCode), color = 'blue', linetype = 'solid') +
                ##geom_point(data = na.omit(lam.tot), aes(x = Age, y = n, colour = PlotCode), shape = 16) +
                ##geom_point(data = na.omit(lam.tot), aes(x = Age, y = n), shape = 16, color = "blue") +
                ##geom_line(data = na.omit(pop_growth_rate.df.l1), aes(x = Age, y = N), size = 1, linetype = 'dotted') +
                geom_line(aes(x = Age, y = N), size = 1.2) +
                geom_hline(yintercept = 1)

  N.p.plot
  # Multiplot
  multiplot.pop <- ggarrange(N.p.plot, lam.p.plot, str.size.plot,
                    ncol = 3, nrow = 1, labels = "AUTO")
  multiplot.pop
  ggsave(multiplot.pop, file = 'plot_pop_panel_classic.pdf', height=4, width=12, units="in")

} 
# vital rates plots
#load('./Data/r_data/plots.RData')
{
 # definitions
 {
    bin.size <- 0.1
    n.sim <- 1000
    #h2.pred.f <- seq(min(f1.rep.data$h2, na.rm = TRUE), ceiling(max(f1.rep.data$h2, na.rm = TRUE)), length.out = m+1)
    h1.pred <- seq(min(survival$h1, na.rm = TRUE), max(survival$h1, na.rm = TRUE), length.out = m)
    h2.pred <- seq(min(survival$h2, na.rm = TRUE), max(survival$h2, na.rm = TRUE), length.out = m+1)
    mids.h2 <- hist(survival$h2, breaks = h2.pred, plot = FALSE)$mids
    ln.mids.h2 <- log(mids.h2)
    lh1.pred <- log(h1.pred)
    data.s.j <- subset(survival, !is.na(h1) & !is.na(sup))
    data.s.j <- droplevels(data.s.j)
    data.g.j <- subset(data.s.j, !is.na(h2))
    data.g.j <- droplevels(data.g.j)
    h1.pred.1 <- seq(min(survival$h1, na.rm = TRUE), ceiling(max(survival$h1, na.rm = TRUE)), length.out = m+1)
    #h1.pred.1 <- seq(min(survival$h1, na.rm = TRUE), max(survival$h1, na.rm = TRUE), length.out = m+1)
    mids <- hist(survival$h1, breaks = h1.pred.1, plot = FALSE)$mids
    ln.mids <- log(mids)
 }
 # Survival
 {
   data.s.j.1 <- subset(data.s.j, sup == 1)
   count.s.j.1 <- hist(data.s.j.1$h1, breaks = h1.pred, plot = FALSE)$counts
   count.s.j <- hist(data.s.j$h1, breaks = h1.pred, plot = FALSE)$counts
   prop.s.j <- count.s.j.1/count.s.j
   plot.data.s.j <- as.data.frame(matrix(ncol = 4, nrow = 4))
   names(plot.data.s.j) <- c("age", "h1", "prop", "n")
   for (a in Age.pred) {
    data.a <- subset(data.s.j, Age == a)
    data.a.1 <- subset(data.a, sup == 1)
    count.a.1 <- hist(data.a.1$h1, breaks = h1.pred.1, plot = FALSE)$counts
    count.a <- hist(data.a$h1, breaks = h1.pred.1, plot = FALSE)$counts
    plot.data.s.j[(a-1)*m+1:m*a, ] <- cbind(rep(a, m), mids, count.a.1/count.a, count.a)
   }
   plot.data.s.j <- subset(plot.data.s.j, !is.na(prop))
   pred.s.j <- as.data.frame(matrix(ncol = 3, nrow = length(Age.pred)*length(mids)))
   names(pred.s.j) <- c("age.pred", "h1.pred", "s.pred")
   for (a in Age.pred) {
    s.pred <- as.vector(predict(gam.s$gam, newdata = data.frame(ln.h1 = ln.mids, Age = rep(a, length(mids))), type = "response", re.form = NA))
    pred.s.j[(a-1)*length(mids)+1:length(mids), ] <- cbind(rep(a, length(mids)), mids, s.pred)
   }
   s.plot <- ggplot(pred.s.j) +
     theme_minimal() +
     #theme(legend.position="none") +
     #theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
     xlab("Successional age (years)") +
     ylab("Height (m)") +
     scale_x_continuous(expand = c(0,0)) +
     scale_y_continuous(expand = c(0,0)) +
     scale_fill_gradientn(colors = hcl.colors(12, "YlOrRd")) +
     #scale_alpha(guide = "none") +
     geom_raster(aes(age.pred, h1.pred, fill = s.pred)) +
     scale_colour_gradient(low = "#d90000", high = "#fcdf03") +
     geom_point(data = plot.data.s.j, aes(age, h1, fill = prop, size = n), color = 'black', shape = 21) +
     scale_size_area(breaks = c(1,5,10,20)) +
     geom_contour(aes(age.pred, h1.pred, z = s.pred), color = 'black') +
     geom_label_contour(aes(age.pred, h1.pred, z = s.pred)) +
     guides(color = guide_colourbar(title = "Survival\nprobability", order = 1), size = guide_legend(title = "Sample\nsize", order = 2))
     #dev.new(width = 10, height = 8)
     s.plot
     #ggsave("survival.pdf", s.plot, device = "pdf", width = 9, height = 6, units = "in", dpi = 180*2)
 }
 # growth
 { 
   # create matrix with observed data
   #h1.pred.g <- seq(min(growth$h1, na.rm = TRUE), ceiling(max(growth$h1, na.rm = TRUE)), length.out = m)
   #h1.pred.1.g <- seq(min(growth$h1, na.rm = TRUE), ceiling(max(growth$h1, na.rm = TRUE)), length.out = m+1)
   h1.pred.g <- seq(min(growth$h1, na.rm = TRUE), max(growth$h1, na.rm = TRUE), length.out = m)
   h1.pred.1.g <- seq(min(growth$h1, na.rm = TRUE), max(growth$h1, na.rm = TRUE), length.out = m+1)
   data.g.mean <- subset(growth, !is.na(h1) & !is.na(h2))
   count.g.1 <- hist(data.g.mean$h1, breaks = h1.pred.g, plot = FALSE)$counts
   data.g.matrix <- matrix(nrow = m, ncol = m, NA)
   for (j in 1:m) {
     g.mean.list <- c()
     for (a in Age.pred) {
       # for first height
       if (j == 1) {
        data.a.j <- subset(data.g.mean, Age == a & h1 <= h1.pred.g[j])
        del.h.list <- c()
        for (i in 1:nrow(data.a.j))
          del.h.list[i] <- data.a.j$h2[i]-data.a.j$h1[i]
        g.mean.a.j <- mean(del.h.list)
        data.g.matrix[1, a] <- g.mean.a.j
       }
       # for the rest
       else {
        data.a.j <- subset(data.g.mean, Age == a &  h1 <= h1.pred.g[j] & h1 > h1.pred.g[j-1])
        del.h.list <- c()
        for (i in 1:nrow(data.a.j))
          del.h.list[i] <- data.a.j$h2[i]-data.a.j$h1[i]
        g.mean.a.j <- mean(del.h.list)
        data.g.matrix[j, a] <- g.mean.a.j
       }
      }
    }
   # create data frame with observed data
   data.g.df <- as.data.frame(matrix(ncol = 4, nrow = 4))
   names(data.g.df) <- c("age", "h1", "gmean", "n")
   for (a in Age.pred) {
    data.a <- subset(data.g.mean, Age == Age.pred[a])
    count.a <- hist(data.a$h1, breaks = h1.pred.1.g, plot = FALSE)$counts
    data.g.df[(a-1)*length(mids)+1:length(mids), ] <- cbind(rep(Age.pred[a], length(mids)), mids, data.g.matrix[ ,a], count.a)
   }
   # create data frame for plotting
   pred.g.mean <- as.data.frame(matrix(ncol = 3, nrow = length(Age.pred)*length(mids)))
   names(pred.g.mean) <- c("age.pred", "h1.pred", "g.mean")
   for (a in 1:length(Age.pred)) {
    g.m <- predict(gam.g$gam, newdata = data.frame(ln.h1 = x.pred, Age = rep(a, m)), type = "response", re.form = NA)
          delta.g <- exp(g.m)-exp(x.pred)
    pred.g.mean[(a-1)*length(mids)+1:length(mids), ] <- cbind(rep(Age.pred[a], length(mids)), mids, delta.g)
   }
   g.cont.seq <- c(0.40,0.30, 0.25, 0, -0.5, -1, -1.5, -2)
   # plot
   g.mean.plot <- ggplot(pred.g.mean) +
   theme_minimal() +
   #theme(legend.position="none") +
   #theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
   xlab("Successional age (years)")+
   ylab("Height (m)") +
   scale_x_continuous(expand = c(0,0)) +
   scale_y_continuous(expand = c(0,0)) +
   scale_fill_gradientn(colours = hcl.colors(12, "YlOrRd")) +
   scale_alpha(guide = "none") +
   geom_raster(aes(age.pred, h1.pred, fill = g.mean)) +
   scale_colour_gradient(low = "#d90000", high = "#fcdf03") +
   geom_point(data = na.omit(data.g.df), aes(age, h1, fill = gmean, size = n), shape = 21, color = 'black') +
   scale_size_area(breaks = c(1,5,10,20)) +
   geom_contour(aes(age.pred, h1.pred, z = g.mean), breaks = g.cont.seq, color = 'black', binwidth = 0.1) +
   geom_label_contour(aes(age.pred, h1.pred, z = g.mean), breaks = g.cont.seq) +
   guides(colour = guide_colourbar(title = "Average\ngrowth (m)", order = 1), size = guide_legend(title = "Sample\nsize", order = 2))
   #dev.new(width = 10, height = 8)
   g.mean.plot
   #ggsave("mean-growth.pdf", g.mean.plot, device = "pdf", width = 9, height = 6, units = "in", dpi = 180*2)
 }
 # Fecundity 1
 {
   data.rep <- subset(reproduction, !is.na(h1) & !is.na(Sup))
   data.rep <- transform(data.rep, id = paste0('p', PlotName, 'i', Tag))
   data.rep$id <- as.factor(data.rep$id)
   seeds <- transform(seeds, id = paste0('p', Plot, 'i', Tag))
   seeds$id <- as.factor(seeds$id)
   mean.nseed <- c()
   for (i in levels(seeds$id)) {
	   seeds.i <- subset(seeds, id == i)
	   seeds.i.mean <- mean(seeds.i$N.seed)
	   mean.nseed <- c(mean.nseed, seeds.i.mean)
   }
   seeds.mean <- data.frame(id = levels(seeds$id), m.nseed = mean.nseed)
   data.rep <- droplevels(data.rep)
   data.rep <- transform(data.rep, nseeds = mean(seeds$N.seed), sourc = 'allmean')
   data.rep <- subset(data.rep, Rep == 1 | Rep == 0)
   for (i in levels(seeds$id)) {
	   if (i %in% levels(data.rep$id)) {
		   if (data.rep[which(data.rep$id == i),]$Rep == 1) {
			   data.rep[which(data.rep$id == i),]$nseeds <- seeds.mean[which(seeds.mean$id == i),]$m.nseed
			   data.rep[which(data.rep$id == i),]$sourc <- 'observed'
			}
	   }
   }
   data.rep$sourc <- as.factor(data.rep$sourc)
   data.rep <- transform(data.rep, totalseeds = Rep*FruitsPerBranch*nseeds)
   data.rep1 <- subset(data.rep, Rep == 1)
   count.rep1 <- hist(data.rep1$h1, breaks = h1.pred, plot = FALSE)$counts
   count.rep <- hist(data.rep$h1, breaks = h1.pred, plot = FALSE)$counts
   prop.rep <- count.rep1/count.rep
   plot.data.rep <- as.data.frame(matrix(ncol = 4, nrow = 4))
   names(plot.data.rep) <- c("age", "h1", "prop.f1", "n")
   for (a in Age.pred) {
       data.rep.a <- subset(data.rep, Age == a)
       data.rep.a.1 <- subset(data.rep.a, Rep == 1)
       count.rep.a.1 <- hist(data.rep.a.1$h1, breaks = h1.pred.1, plot = FALSE)$counts
       count.rep.a <- hist(data.rep.a$h1, breaks = h1.pred.1, plot = FALSE)$counts
       plot.data.rep[(a-1)*m+1:m, ] <- cbind(rep(a, m), mids, count.rep.a.1/count.rep.a, count.rep.a)
       }
   f1.rep.data <- subset(plot.data.rep, !is.na(prop.f1))
   h1.pred.f <- seq(min(f1.rep.data$h1, na.rm = TRUE), ceiling(max(f1.rep.data$h1, na.rm = TRUE)), length.out = m)
   pred.f1.j <- as.data.frame(matrix(ncol = 3, nrow = length(Age.pred)*length(mids)))
   names(pred.f1.j) <- c("age.pred", "h1.pred", "f1.pred")
   for (a in Age.pred) {
    pred.f1 <- as.vector(predict(gam.f1$gam, newdata = data.frame(ln.h1 = ln.mids, Age = rep(a, m)), type = "response", re.form = NA))
    pred.f1.j[(a-1)*m+1:m, ] <- cbind(rep(a, m), mids, pred.f1)
    }
   f1.plot <- ggplot(pred.f1.j) +
   theme_minimal() +
   #theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
   xlab("Successional age (years)") +
   ylab("Height (m)") +
   scale_x_continuous(expand = c(0,0)) +
   scale_y_continuous(expand = c(0,0)) +
   scale_fill_gradientn(colours = hcl.colors(12, "YlOrRd")) +
   scale_alpha(guide = "none") +
   geom_raster(aes(age.pred, h1.pred, fill = f1.pred)) +
   scale_colour_gradient(low = "#d90000", high = "#fcdf03") +
   geom_point(data = f1.rep.data, aes(age, h1, fill = prop.f1, size = n), color = 'black', shape = 21) +
   geom_contour(aes(age.pred, h1.pred, z = f1.pred), color = 'black', binwidth = 0.1) +
   geom_label_contour(aes(age.pred, h1.pred, z = f1.pred)) +
   guides(colour = guide_colourbar(title = "Reproduction\nprobability", order = 1), size = guide_legend(title = "Sample\nsize", order = 2))
   #dev.new(width = 10, height = 8)
   f1.plot
   #ggsave("f1.pdf", f1.plot, device = "pdf", width = 9, height = 6, units = "in", dpi = 180*2)
 }
 # Fecundity 2
 {
   plot.data.f2 <- as.data.frame(matrix(ncol = 3, nrow = m*Age.mature))
   names(plot.data.f2) <- c("age.pred.f2", "h1.pred.f2", "frt.n.pred")
   for (a in Age.pred) {
       frt.n.pred <- as.vector(fitted(gam.f2, newdata = data.frame(ln.h1 = x.pred, Age = rep(a, m)), type = "response", re_formula = NA)[,1])
       plot.data.f2[(a-1)*m+1:m, ] <- cbind(rep(a, m), mids, frt.n.pred)
       }
   f2.rep.data <- subset(reproduction, Rep == 1)
   my_varsf2 <- c("Age", "h1", "TotalFruitset")
   f2.rep.data <- f2.rep.data[my_varsf2]
   f2.plot <- ggplot(plot.data.f2) +
     #theme_minimal() +
     #theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
     xlab("Successional age (years)") +
     ylab("Height (m)") +
     scale_x_continuous(expand = c(0,0)) +
     scale_y_continuous(expand = c(0,0)) +
     scale_fill_gradientn(colours = hcl.colors(12, "YlOrRd")) +
     scale_alpha(guide = "none") +
     geom_raster(aes(age.pred.f2, h1.pred.f2, fill = frt.n.pred)) +
     scale_colour_gradient(low = "#d90000", high = "#fcdf03") +
     geom_point(data = f2.rep.data, aes(Age, h1, fill = TotalFruitset), size = 2, shape = 21, color = 'black') +
     geom_contour(aes(age.pred.f2, h1.pred.f2, z = frt.n.pred), color = 'black') +
     geom_label_contour(aes(age.pred.f2, h1.pred.f2, z = frt.n.pred)) +
     guides(colour = guide_colourbar(title = "Number of\nfruits per tree"))
     #dev.new(width = 10, height = 8)
     f2.plot
     #ggsave("f2.pdf", f2.plot, device = "pdf", width = 9, height = 6, units = "in", dpi = 180*2)
 }
 # Fecundity 3
 {
   plot.data.f3 <- as.data.frame(matrix(ncol = 3, nrow = m*Age.mature))
   names(plot.data.f3) <- c("Age.pred", "h1.pred.f3", "sd.n.pred")
   ln.mids.f3 <- seq(min(seeds$ln.h1), max(seeds$ln.h1), length.out = m)
   age.mids.f3 <- seq(min(seeds$Age), max(seeds$Age))
   for (a in Age.pred) {
       sd.n.pred <- as.vector(fitted(gam.f3, newdata = data.frame(ln.h1 = x.pred, Age = rep(a, m)), type = "response", re_formula = NA)[,1])
       plot.data.f3[(a-1)*m+1:m, ] <- cbind(rep(a, m), mids, sd.n.pred)
       }
   f3.rep.data <- seeds
   my_varsf3 <- c("Age", "h1", "N.seed")
   f3.rep.data <- f3.rep.data[my_varsf3]
   f3.plot <- ggplot(plot.data.f3) +
   theme_minimal() +
   #theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
   xlab("Successional age (years)") +
   ylab("Height (m)") +
   scale_x_continuous(expand = c(0,0)) +
   scale_y_continuous(expand = c(0,0)) +
   scale_fill_gradientn(colours = hcl.colors(12, "YlOrRd")) +
   scale_alpha(guide = "none") +
   geom_raster(aes(Age.pred, h1.pred.f3, fill = sd.n.pred)) +
   scale_colour_gradient(low = "#d90000", high = "#fcdf03") +
   geom_point(data = f3.rep.data, aes(Age, h1, fill = N.seed), size = 2, color = 'black', shape = 21) +
   geom_contour(aes(Age.pred, h1.pred.f3, z = sd.n.pred), color = 'black') +
   geom_label_contour(aes(Age.pred, h1.pred.f3, z = sd.n.pred)) +
   guides(colour = guide_colourbar(title = "Number of\nseeds per fruit"))
   #dev.new(width = 10, height = 8)
   f3.plot
   #ggsave("f3.pdf", g.f3.j, device = "pdf", width = 9, height = 6, units = "in", dpi = 180*2)
 }
 # Fecundity 1*2*3
 {
    f1.rep.data <- transform(f1.rep.data, exito = prop.f1*n/2)
    plot.data.f123 <- as.data.frame(matrix(ncol = 3, nrow = m*Age.mature))
    names(plot.data.f123) <- c("Age.pred", "h1.pred", "f123.pred")
    for (a in Age.pred) {
        f1.pred <- as.vector(predict(gam.f1$gam, newdata = data.frame(ln.h1 = log(h1.pred.f), Age = rep(a, m)), type = "response", re_formula = NA))
        f2.pred <- as.vector(fitted(gam.f2, newdata = data.frame(ln.h1 = log(h1.pred.f), Age = rep(a, m)), type = "response", re_formula = NA)[,1])
        f3.pred <- as.vector(fitted(gam.f3, newdata = data.frame(ln.h1 = log(h1.pred.f), Age = rep(a, m)), type = "response", re_formula = NA)[,1])
        f123 <- f1.pred*f2.pred*f3.pred
        plot.data.f123[(a-1)*m+1:m, ] <- cbind(rep(a, m), h1.pred.f, f123)
    }
    f123.cont.seq <- c(10,20,30,40,50,100,150,200,250,300,350,400)
    f123.plot <- ggplot(plot.data.f123) + ## 1
        theme_minimal() +
        #theme(legend.position="none") +
        #theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
        xlab("Successional age (years)") +
        ylab("Height (m)") +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        scale_fill_gradientn(colours = hcl.colors(12, "YlOrRd")) +
        scale_alpha(guide = "none") +
        geom_raster(aes(Age.pred, h1.pred, fill = f123.pred)) +
        scale_colour_gradient(low = "#d90000", high = "#fcdf03") +
        geom_point(data = data.rep, aes(Age, h1, fill = totalseeds), color = 'black', shape = 21, size = 3) +
        geom_contour(aes(Age.pred, h1.pred, z = f123.pred), breaks = f123.cont.seq, color = 'black') +
        geom_label_contour(aes(Age.pred, h1.pred, z = f123.pred), breaks = f123.cont.seq) +
        scale_fill_gradient(low = "#d90000", high = "#fcdf03") 
        #new_scale_fill() +
        #geom_point(data = f2.rep.data, aes(Age, h1, fill = TotalFruitset), size = 2, shape = 22, color = 'black') +
        #scale_fill_gradient(low = "#d90000", high = "#fcdf03") 
        #new_scale_fill() +
        #geom_point(data = f3.rep.data, aes(Age, h1, fill = N.seed), size = 2, color = 'black', shape = 23) +
        #scale_fill_gradient(low = "#d90000", high = "#fcdf03")
        #dev.new(width = 10, height = 8)
    f123.plot
    #ggsave("f123-1.pdf", f123.plot.1, device = "pdf", width = 9, height = 6, units = "in", dpi = 180*2)
 }
 # Fecundity 4
 {
  f4.plot <- ggplot(tot_seeds_line, aes(x = Age, y = Establishment)) +
            theme_minimal() +
            #theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
            xlab('Succesional age (years)') +
            ylab('Establishment probability') +
            geom_line(size = 1, color = 'red', alpha = 0.9) +
            geom_point(data = tot_seeds_t1, aes(x = Age, y = Establishment), shape = 21)
            #dev.new(width = 10, height = 8)
   f4.plot
   #ggsave("f4.pdf", f4.plot, device = "pdf", width = 9, height = 6, units = "in", dpi = 180*2)
 }
 # Fecundity 6
 {
    f6.plot <- ggplot(f6.line.df, aes(x = Age, y = Nsprouts)) +
        theme_minimal() +
        ##scale_y_continuous(trans = "log10") +
        #theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
        xlab("Succesional age (years)") +
        ylab("Number of sprouts") +
        geom_line(size = 1, color = 'red', alpha = 0.9) +
        geom_point(data = nsprouts, aes(x = Age, y = nscaled), shape = 21)
        #dev.new(width = 10, height = 8)
    f6.plot
   #ggsave("f6.pdf", f6.plot, device = "pdf", width = 9, height = 6, units = "in", dpi = 180*2)
 }
 # Density (F5 & F7)
 {
  f57.plot <- ggplot(f57.df, aes(x=height, color=group, fill=group)) + 
  theme_minimal() +
  #theme(legend.position="none") +
  xlab("Height (m)") +
  ylab("Density") +
  geom_histogram(aes(y=..density..), alpha=0.5, 
                position="identity", binwidth = 0.06) +
  geom_density(alpha=.2)
  #dev.new(width = 10, height = 8)
  f57.plot
  #f57.gden <- ggdensity(f57.df, x = "height",
  #add = "mean", rug = TRUE, color = "group", fill = "group")
  #f57.gden
  
  ggplot(f57.df, aes(x=height, fill=group)) + geom_density(alpha=.3)
 }
  # Multiplot
  multiplot <- ggarrange(s.plot, g.mean.plot, f123.plot, f57.plot, f4.plot, f6.plot,
                    ncol = 2, nrow = 3, labels = "AUTO")
  multiplot
  ggsave(multiplot, file = 'plot_panel.pdf', height=10, width=8, units="in")
  ##
}


































{
  # building gam for understory n
  understory_n.gam_model.1 <- gam(round(stb.n$SNB)~s(log(x+1), k = 3), family = nb())
  understory_n.gam_model.2 <- gam(round(stb.n$SNB)~s(log(x+1), k = 3), family = poisson)
  understory_n.gam_model.3 <- gam(round(stb.n$SNB)~s(log(x+1), k = 5), family = nb())
  understory_n.gam_model.4 <- gam(round(stb.n$SNB)~s(log(x+1), k = 7), family = nb())
  understory_n.gam_model.5 <- gam(round(stb.n$SNB)~s(log(x+1), k = 10), family = nb())
  #> AICc(understory_n.gam_model.1)
  # [1] 443.6209
  #> AICc(understory_n.gam_model.2)
  #[1] 1528.83
  #> AICc(understory_n.gam_model.3)
  #[1] 436.3641
  #> AICc(understory_n.gam_model.4)
  #[1] 438.608
  #> AICc(understory_n.gam_model.5)
  #[1] 439.3461
  #
  plot(0:100, predict(understory_n.gam_model.1, newdata= data.frame(x=0:100), type = 'response'), type = 'l', col = 'black')
  lines(0:100, predict(understory_n.gam_model.2, newdata= data.frame(x=0:100), type = 'response'), col = 'blue')
  lines(0:100, predict(understory_n.gam_model.3, newdata= data.frame(x=0:100), type = 'response'), col = 'red')
  lines(0:100, predict(understory_n.gam_model.4, newdata= data.frame(x=0:100), type = 'response'), col = 'green')
  lines(0:100, predict(understory_n.gam_model.5, newdata= data.frame(x=0:100), type = 'response'), col = 'pink')
  points(x, stb.n$SNB)
  understory_n <- predict(understory_n.gam_model.3, newdata = data.frame(x = 0:100), type = "response")
}
# total predicted lambda
{
  lam.list <- c()
  n.0 <- understory_n[1]
  init.n.a.v <- n.0*F5
  size.v.a.NM <- list(init.n.a.v) # create list for size structure change, NM = no migration
  n.list <- c(sum(init.n.a.v))
  for (a in 1:Age.mature) {
       		n.a.v <- k.i.j.a[a,,]%*%init.n.a.v
       		lam.a <- sum(n.a.v)/sum(init.n.a.v)
       		size.v.a.NM[[a+1]] <- n.a.v
        	n.list <- c(n.list, sum(n.a.v))
       		lam.list <- c(lam.list, lam.a)
       		init.n.a.v <- n.a.v
  	}
  tot.lam.pred <- lam.list.NM <- lam.list # lam.list is the transitory lambda vector
  prd.N.total <- data.frame(N = n.list, Age = c(0,Age.pred))
  #pdf('lambda.pdf')
  plot(Age.pred, tot.lam.pred, type = 'l')
  #dev.off()
}
# scaling observed lambda
{
  canopy <- read.csv("./Data/Additional/canopy_cleaned.csv")
  canopy$plot <- as.factor(as.character(canopy$plot))
  canopy$Census <- as.factor(as.character(canopy$Census))
  total_n <- c()
  canopy_n <- c()
  total_l <- c()
  # getting n
  for (a in 0:100) {
    canopy.a <- droplevels(subset(canopy, Age == a & sup == 1))
    if (a == 100) {
      can.n <- c()
      for (i in levels(canopy.a$Census))
        canopy.a.i <- subset(canopy.a, Census == i)
        can.n <- c(nrow(canopy.a.i), can.n)
      canopy_n[a+1] <- mean(can.n)
    } else {
      if (nlevels(canopy.a$plot) > 1) {
        canopy_n[a+1] <- sum(nrow(canopy.a)) / nlevels(canopy.a$plot)
      } else {
        canopy_n[a+1] <- sum(nrow(canopy.a))
      }
    }
    total_n[a+1] <- understory_n[a+1] + canopy_n[a+1] 
  }
  # plotting n
  #pdf('n_change_0.pdf')
  plot(0:100, total_n, type = 'b', col = 'blue', xlab = 'Succesional age (years)', ylab = 'N')
  lines(0:100, canopy_n, type = 'b', col = 'red')
  lines(0:100, understory_n, type = 'b',col = 'green')
  lines(0:100, rep(0, 101), col = 'black')
  lines(0:100, n.list, col = 'goldenrod3')
  #dev.off()
}
# total observed lambda
{
  ob.lam.list <- c()				# create observed lambda vector
  N.0 <- subset(s1.db, Age == 0)		# get row number of trees in first year
  N.0 <- subset(N.0, sup == 1)			# filtering to obtain only alive trees
  N.init <- nrow(N.0)				# get number of alive trees in first year
  N.list <- c(N.init)				# create a new list of observed population sizes and setting first value
  for (i in 1:(Age.mature-1)) {
    ns.i <- subset(s1.db, Age == i)	# create a set of trees in year a, starting from a = 1
    ns.i <- subset(ns.i, sup == 1)	# filtering to get only alive trees
    N.ns.i <- nrow(ns.i)			# counting trees
    if (N.ns.i != 0 & N.init != 0) {	# if it is not an empty vector:
     lambda <- N.ns.i/N.init	# get lambda
    } else {			# if it is:
     lambda <- NA		# set lambda to NA (avoiding to get 0 instead of NA)
    }
      N.init <- N.ns.i			# set a new initial population for next lambda
      N.list <- c(N.list, N.init)		# add new population size to vector
      ob.lam.list <- c(ob.lam.list, lambda)# aff new lambda to observed lambda vector
  }
  ns.mature <- subset(s1.db, Age == Age.mature) 	# get mature forest subset (avoid population overestimating due to individuals repeat at Age.mature)
  ns.mature <- subset(ns.mature, sup == 1)		# get alive trees only
  ns.mature <- transform(ns.mature, Census = as.factor(as.character(Census))) # transform census year as factor
  ns.mature <- droplevels(ns.mature) # drop useless levels
  N.mat.list <- c()  # create new vector of population sizes
  for (i in levels(ns.mature$Census)) {		# for each year in mature forest:
    ns.mat.yr <- subset(ns.mature, Census == i)	# get a subset of population in census year = i
    ns.mat.N <- nrow(ns.mat.yr)			# count trees
    N.mat.list <- c(N.mat.list, ns.mat.N)		# add population size to vector
  }
  lam.m.list <- c()					# create a new vector of lambdas of mature forest
  init.n.mat <- N.mat.list[1]			# set first population size
  for (j in 2:length(N.mat.list)) {			# for every j in population size list:
    lam.j <- N.mat.list[j]/init.n.mat		# get lambda
    lam.m.list <- c(lam.m.list,lam.j)		# add lambda to list
    init.n.mat <- N.mat.list[j]			# reset initial population size
  }
  N.list[100] <- round(mean(N.mat.list))
  lam.mature <- geoMean(lam.m.list)			# apply geometric mean to get mature forest lambda
  ob.lam.list <- c(ob.lam.list, lam.mature)		# add observed lambda in mature forest to the whole succession list
  ob.lam.df <- as.data.frame(matrix(ncol = 2, nrow = Age.mature)) # create database for observed lambdas
  names(ob.lam.df) <- c("Age", "ob.lambda")		# naming database
  ob.lam.df$Age <- Age.pred			# set age range
  ob.lam.df$ob.lambda[Age.pred] <- ob.lam.list # set observed lambda value
  lambda.df <- as.data.frame(list(lambda = lam.list, Age = 1:Age.mature)) # create lambda data frame for plottin
}
# predicted lambda by plot
{
  load("./Data/kernel-plot-list.RData")					# read list of kernels calculated by plot: poner link del código?
  s1.db <- transform(s1.db, PLOT = as.factor(as.character(plot)))		# transform plot variable to factor
  nplot <- nlevels(s1.db$PLOT)						# get number of plots
  est.lam.df <- as.data.frame(matrix(nrow = Age.mature, ncol = nplot, NA))	# create new database for estimated lambda values
  names(est.lam.df) <- 1:nplot						# plot names
  for (i in 1:nplot) {								# for each plot:
    Mim.s.p <- subset(s1.db, PLOT == i & !is.na(ln.h2))			# get a subset of plot i omitting NA's
    if (Mim.s.p$PLOT[1] != 5) {						# if plot number is 5 (mature forest):
      if (all(!is.na(Mim.s.p$Age))) {					# if it is not an empty list:
        lam.list <- c()							# create new empty list for lambda values
        min.a <- min(Mim.s.p$Age, na.rm = TRUE)				# set maximum height for that plot
        max.a <- max(Mim.s.p$Age, na.rm = TRUE)				# set minimum height for that plot
        n.1.i <- which(Mim.s.p$Age == min.a)				# get row number of trees
        n.1.h.i <- log(Mim.s.p[n.1.i,]$h2)[order(log(Mim.s.p[n.1.i,]$h2))] # get height of trees
        init.n.a.v.i <- n.1.v.i <- hist(n.1.h.i, breaks = e.pred, plot = FALSE)$counts # get size structure vector and set initial vector
        for (a in ((min.a+1):max.a)) {					# for every year:
          n.a.v.i <- k.p.list[[as.integer(i)]][(a-min.a),,]%*%init.n.a.v.i # multiply size vector by plot kernel
          lam.a.i <- sum(n.a.v.i)/sum(init.n.a.v.i)			# get lambda
          init.n.a.v.i <- n.a.v.i					# reset initial vector
          lam.list <- c(lam.list, lam.a.i)				# add lambda to list
        }
        est.lam.df[(min.a+1):max.a, as.numeric(i)] <- lam.list		# fill database
      }
    }
  }
  est.lam.list <- list(NA)							# create new list
  for (i in 1:nplot) {							# for each plot
    if (i != 5) {								# if nor plot 5 (not mature forest):
      plot.l.i <- as.data.frame(matrix(nrow = length(which(!is.na(est.lam.df[,i]) == TRUE)), ncol = 2)) # create new database for that plot
      names(plot.l.i) <- c("Age", "lambda")					# name columns
      plot.l.i$Age <- which(!is.na(est.lam.df[,i]))				# set ages
      plot.l.i$lambda <- est.lam.df[,i][which(!is.na(est.lam.df[,i]) == TRUE)]# set lambda values
      est.lam.list[[i]] <- plot.l.i 						# set plot number
    } else {								# if it is plot 5:
      plot.l.5 <- as.data.frame(matrix(nrow = 1, ncol = 2))			# create new database
      names(plot.l.5) <- c("Age", "lambda")					# name columns
      plot.l.5$Age[1] <- Age.mature						# set age
      plot.l.5$lambda[1] <- lam.mature					# set lambda value
      est.lam.list[[5]] <- plot.l.5						# set plot
    }
  }
}
# observed lambda by plot
{
  nplot <- nlevels(s1.db$PLOT)						# get number of plots
  s1.db <- transform(s1.db, PLOT = as.factor(as.character(plot)))	# transoform plots as factor
  plot.lam.list <- as.data.frame(matrix(nrow = Age.mature, ncol = length(levels(s1.db$PLOT)), NA)) # create new database for observed lambda list by plot
  names(plot.lam.list) <- 1:nplot					# name plot columns
  plot.N.list <- as.data.frame(matrix(nrow = Age.mature, ncol = length(levels(s1.db$PLOT)), NA)) # create new database for observed lambda list by plot
  names(plot.N.list) <- 1:nplot					# name plot columns
  for (i in levels(s1.db$PLOT)) {					# for every plot:
    Mim.s.p <- subset(s1.db, PLOT == i)				# get a subset of the whole database
    if (Mim.s.p$PLOT[1] != 5) {					# if plot number is not 5: !!! (tengo que ver por qué, está vacío)
      if (all(!is.na(Mim.s.p$Age))) {				# if subset has all ages registered:
        lam.list <- c()						# create a new lambda list
        n.list.p <- c()
        min.a <- min(Mim.s.p$Age, na.rm = TRUE)			# get maximum age
        max.a <- max(Mim.s.p$Age, na.rm = TRUE)			# get minimum age
        n.1.i <- which(Mim.s.p$Age == min.a)			# get rows containing trees in year a
        n.1.h.i <- log(Mim.s.p[n.1.i,]$h2)[order(log(Mim.s.p[n.1.i,]$h2))] # get tree hight and order
        init.n.a.v.i <- n.1.v.i <- hist(n.1.h.i, breaks = e.pred, plot = FALSE)$counts # count trees in each size class and set a population structure vector
        for (a in ((min.a+1):max.a)) {				# for each age in plot subset:
          n.a.db <- droplevels(subset(Mim.s.p, Age == a-1))	# drop levels
          n.a.v.i <- nlevels(n.a.db$id)				# count trees
          n.list.p <- c(n.list.p, n.a.v.i)
          lam.a.i <- n.a.v.i/sum(init.n.a.v.i)			# get lambda
          init.n.a.v.i <- n.a.v.i				# reset initial n
          lam.list <- c(lam.list, lam.a.i)			# add lambda to list
        }
        plot.lam.list[(min.a+1):max.a, as.numeric(i)] <- lam.list	# set lambda values and plot
        plot.N.list[(min.a+1):max.a, as.numeric(i)] <- n.list.p	# set lambda values and plot
      }
    }
  }
  nplot <- nlevels(s1.db$PLOT)
  obs.lam.plot <- list(NA)						# create lambda list
  obs.N.plot <- list(NA)
  for (i in 1:nplot) {							# for each plot:
    if (i != 5) {							# if it is not 5:
      plot.l.i <- as.data.frame(matrix(nrow = length(which(!is.na(plot.lam.list[,i]) == TRUE)), ncol = 2)) # create new database: cols = 2 & rows = number of not NA in list
      names(plot.l.i) <- c("Age", "lambda")			# set names
      plot.l.i$Age <- which(!is.na(plot.lam.list[,i]))	# set valid ages
      plot.l.i$lambda <- plot.lam.list[,i][which(!is.na(plot.lam.list[,i]) == TRUE)] # set lambda values
      obs.lam.plot[[i]] <- plot.l.i 				# add lambda matrix to list
      plot.n.i <- as.data.frame(matrix(nrow = length(which(!is.na(plot.N.list[,i]) == TRUE)), ncol = 2)) # create new database: cols = 2 & rows = number of not NA in list
      names(plot.n.i) <- c("Age", "N")			# set names
      plot.n.i$Age <- which(!is.na(plot.N.list[,i]))	# set valid ages
      plot.n.i$N <- plot.N.list[,i][which(!is.na(plot.N.list[,i]) == TRUE)] # set lambda values
      obs.N.plot[[i]] <- plot.n.i
    }
    else if (i == 5) {						# if plot is 5:
      plot.l.5 <- as.data.frame(matrix(nrow = 1, ncol = 2))	# create new database for mature forest
      names(plot.l.5) <- c("Age", "lambda")			# name it
      plot.l.5$Age[1] <- Age.mature				# set age as mature
      plot.l.5$lambda[1] <- lam.mature			# set lambda as calculated for mature forest above
      obs.lam.plot[[5]] <- plot.l.5				# set plot number
      plot.n.5 <- as.data.frame(matrix(nrow = 1, ncol = 2))	# create new database for mature forest
      names(plot.n.5) <- c("Age", "N")			# name it
      plot.n.5$Age[1] <- Age.mature				# set age as mature
      plot.n.5$N[1] <- N.list[100] 				# for mature forest above
      obs.N.plot[[5]] <- plot.n.5				# set plot number
    }
  }
}
# migration estimating function (for population size)
{
  obs.N.total <- data.frame(Age.pred,N.list)
  names(obs.N.total) <- c("Age", "N")
  distance.N <- function(parm) {
    dist.N <- c()
    n.0 <- which(s1.db$Age == 0)	#
    n.0.h <- log(s1.db[n.0,]$h2)	#
    n.0.v <- hist(n.0.h, breaks = e.pred, plot = FALSE)$counts
    b0 <- parm
    c <- b0
    init.n.a.v <- n.0.v+c*F5
    pred.N.total <- c() # size structure vector by year, WM = with migration
    for (a in 1:Age.mature) {
      n.a.v <- k.i.j.a[a,,]%*%init.n.a.v
      n.a.v <- n.a.v+c*F5
      pred.N.total <- c(pred.N.total, sum(n.a.v))
      init.n.a.v <- n.a.v			# add size structure
    }
    for (p in 1:nplot) {
      dist.N.p <- 0
      N.p <- obs.N.plot[[p]]
      max.a <- max(N.p$Age)
      for (j in 1:nrow(N.p))
      dist.N.p <- dist.N.p + (N.p$N[j]-pred.N.total[max.a-nrow(N.p)+j])^2
      dist.N <- c(dist.N, dist.N.p)
    }
    dist.N <- sqrt(sum(dist.N))
    cat(paste0("b0 = ", b0, ", dist = ", dist.N, "\n"))
    return(dist.N)
  }
}
# migration estimating function (for lambda)
{
  distance.L <- function(c) {
    dist.L <- c()
    n.0 <- which(s1.db$Age == 0)	# get row numbers with age = 0
    n.0.h <- log(s1.db[n.0,]$h2)	#
    n.0.v <- hist(n.0.h, breaks = e.pred, plot = FALSE)$counts
    init.n.a.v <- n.0.v+c*F5
    pred.L.total <- c() # size structure vector by year, WM = with migration
    for (a in 1:Age.mature) {
      n.a.v <- k.i.j.a[a,,]%*%init.n.a.v
      n.a.v <- n.a.v+c*F5
      lam.p <- sum(n.a.v)/sum(init.n.a.v)
      pred.L.total <- c(pred.L.total, lam.p)
      init.n.a.v <- n.a.v			# add size structure
    }
    for (p in 1:nplot) {
      dist.L.p <- 0
      L.p <- obs.lam.plot[[p]]
      max.a <- max(L.p$Age)
      for (j in 1:nrow(L.p))
      dist.L.p <- dist.L.p + (L.p$lambda[j]-pred.L.total[max.a-nrow(L.p)+j])^2
      dist.L.p <- sqrt(dist.L.p)
      dist.L <- c(dist.L, dist.L.p)
    }
    dist.L <- mean(dist.L)
    cat(paste0("c = ", c, ", dist = ", dist.L, "\n"))
    return(dist.L)
  }
}
# trying different estimators (for population size)
{
  opt.brent <- optim(0, distance.N, method = "Brent", lower = 0, upper = 100)
  c <- opt.brent$par # c = 9.35491439069071, dist = 81.8803591350505
  opt.lbf <- optim(c(5, -2, -2), distance.N, method = "L-BFGS-B", lower = c(0, -15, -15), upper = c(50, 5, 5))
  par_c <- opt.lbf$par # c = 9.35391449733551, dist = 81.8803596331657
  opt.bfgs <- optim(0, distance.N, method = "BFGS")
  opt.bfgs$par # c = 9.35491911006061, dist = 81.8803591350612
}
# trying different estimators (for lambda)
{
  opt.brent <- optim(0, distance.L, method = "Brent", lower = 0, upper = 100000000)
  opt.brent$par # c = 2.47619435505125, dist = 0.530052570070713
  opt.lbf <- optim(0, distance.L, method = "L-BFGS-B")
  opt.lbf$par # c = 2.47519458486379, dist = 0.530052570638488
  opt.bfgs <- optim(0, distance.L, method = "BFGS")
  opt.bfgs$par # c = 15.7093246924741, dist = 0.535953969391775
}
# migration
{
  # total, no migration
  {
     lambda.df.NM <- as.data.frame(list(lambda = lam.list.NM, Age = 1:Age.mature))
     lambda.edad <- qplot(x = Age, y = log(lambda), data = lambda.df.NM, color = "red", geom = "line", xlab = "Succesional age (years)", ylab = expression(italic(r))) +
      theme_minimal() +
      theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15), legend.text = element_text(size = 12), legend.title = element_text(size = 15), legend.position = "none") +
      geom_point(data = ob.lam.df, mapping = aes(x = Age, y = log(ob.lambda)))
      lambda.edad
      #ggsave("lambda.a.pdf", lambda.edad, device = "pdf", width = 9, height = 6, units = "in", dpi = 180*2)
  }
  # by plot, no migration
  {
    lambda.df.2.NM <- transform(lambda.df.NM)
    p.lam.p.2 <- as.data.frame(matrix(nrow = 0, ncol = 2))
    names(p.lam.p.2) <- c("Age", "lambda")
    for (i in 1:nplot) {
      p.lam.p.i <-transform(obs.lam.plot[[i]])
      p.lam.p.2  <- rbind(p.lam.p.2, p.lam.p.i)
    }
    p.lam.ob <- ggplot(p.lam.p.2, aes(x = Age, y = lambda, alpha = 1/2), show.legend = FALSE) +
      theme_minimal() +
      geom_point(size = 1) +
      theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 15)) +
      geom_line(data = lambda.df.NM, aes(x = Age, y = lambda), col = "red", size = 2, alpha = 1/3, show.legend = FALSE) +
      geom_line(data = ob.lam.df, aes(x = Age, y = ob.lambda), col = "black", size = 2, alpha = 1/3, show.legend = FALSE) +
      labs(x = expression(paste("Abandonment time ", italic(t), " (years)")), y = expression(lambda))+
      scale_alpha(guide = "none") +
      scale_y_continuous(limits = c(0, 3))
      p.lam.ob
      #ggsave("no-migration.pdf", p.lam.ob, device = "pdf", width = 9, height = 6, units = "in", dpi = 180*2)
  }
  # by plot, with migration
  {
    obs.lam.p <- as.data.frame(matrix(nrow = 0, ncol = 3))
    names(obs.lam.p) <- c("Age", "lambda", "plot")
    for (i in 1:nplot) {
      obs.lam.p.i <-transform(obs.lam.plot[[i]], plot = as.factor(i))
      obs.lam.p  <- rbind(obs.lam.p, obs.lam.p.i)
    }
    # calculate new lambda list considering migration
    lam.list <- c()			# create lambda vector
    n.0 <- which(s1.db$Age == 0)	# get row numbers with age = 0
    n.0.h <- log(s1.db[n.0,]$h2)	#
    n.0.v <- hist(n.0.h, breaks = e.pred, plot = FALSE)$counts
    c <- c_0
    init.n.a.v <- n.0.v+c*F5
    size.v.a.WM <- list(init.n.a.v) # size structure vector by year, WM = with migration
    n.list <- c(sum(init.n.a.v))
    for (a in 1:Age.mature) {
      n.a.v <- k.i.j.a[a,,]%*%init.n.a.v
      lam.a <- sum(n.a.v)/sum(init.n.a.v) # add c individuals
      size.v.a.WM[[a+1]] <- n.a.v
      n.list <- c(n.list, sum(n.a.v))
      lam.list <- c(lam.list, lam.a)
      init.n.a.v <- n.a.v			# add size structure
    }
    tot.lam.pred <- lam.list
    Age = rep(Age.pred,2)
    lambda = c(lambda.df.NM$lambda, lam.list)
    Migration = c(rep("no",100),rep("yes",100))
    lambda.df.2 <- data.frame(Age, lambda, Migration)
    plot.l.graf.c <- ggplot(lambda.df.2, aes(x = Age, y = lambda)) +
      theme_minimal() +
      xlab(expression(paste("Abandonment time ", italic(t), " (years)"))) +
      ylab(expression(lambda)) +
      theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 15)) +
      geom_line(size = 2, alpha = 0.8, aes(color=Migration)) +
      geom_point(data = obs.lam.p, aes(x = Age, y = lambda), alpha = 1/2, size = 1) +
      geom_line(data = ob.lam.df, aes(x = Age, y = ob.lambda), col = "black", size = 2, alpha = 1/3, show.legend = FALSE) +
      scale_alpha(guide = "none") +
      scale_y_continuous(limits = c(0, 3))
      plot.l.graf.c
      #ggsave("with-migration.pdf", plot.l.graf.c, device = "pdf", width = 9, height = 6, units = "in", dpi = 180*2)
  }
}
# sensitivity
{
  # function
  {
    sensitivity <- function(mat, string, mod) {
      if (string == 'S') {
        k.i.j.a <- s.i.j.a <- g.i.j.a <- p.i.j.a <- f.i.j.a <- F5.i.j.a <- array(NA, dim = c(Age.mature, m, m))
        s.i.a <- mat*mod
        f.i.a <- matrix(NA, nrow = m, ncol = Age.mature)
        for (a in 1:Age.mature) {
          for (i in 1:m) {
            f.i.a[i, a] <- F1[i, a]*F2[i, a]*F3[i, a]*F4[a]
            for (j in 1:m) {
              p.i.j.a[a, j, i] <- s.i.a[i, a]*G[a, j, i]
              f.i.j.a[a, j, i] <- f.i.a[i, a]*F5[j]
              k.i.j.a[a, j, i] <- p.i.j.a[a, j, i]+f.i.j.a[a, j, i]
            }
          }
        }
      }
      else if (string ==  'G') {
        k.i.j.a <- s.i.j.a <- g.i.j.a <- p.i.j.a <- f.i.j.a <- F5.i.j.a <- array(NA, dim = c(Age.mature, m, m))
        s.i.a <- S
        f.i.a <- matrix(NA, nrow = m, ncol = Age.mature)
        for (a in 1:Age.mature) {
          for (i in 1:m) {
            f.i.a[i, a] <- F1[i, a]*F2[i, a]*F3[i, a]*F4[a]
            for (j in 1:m) {
              p.i.j.a[a, j, i] <- s.i.a[i, a]*(mat[a, j, i]*mod)
              f.i.j.a[a, j, i] <- f.i.a[i, a]*F5[j]
              k.i.j.a[a, j, i] <- p.i.j.a[a, j, i]+f.i.j.a[a, j, i]
            }
          }
        }
      }
      else if (string == 'F1') {
        k.i.j.a <- s.i.j.a <- g.i.j.a <- p.i.j.a <- f.i.j.a <- F5.i.j.a <- array(NA, dim = c(Age.mature, m, m))
        s.i.a <- S
        f.i.a <- matrix(NA, nrow = m, ncol = Age.mature)
        for (a in 1:Age.mature) {
          for (i in 1:m) {
            f.i.a[i, a] <- (mat[i, a]*mod)*F2[i, a]*F3[i, a]*F4[a]
            for (j in 1:m) {
              p.i.j.a[a, j, i] <- s.i.a[i, a]*G[a, j, i]
              f.i.j.a[a, j, i] <- f.i.a[i, a]*F5[j]
              k.i.j.a[a, j, i] <- p.i.j.a[a, j, i]+f.i.j.a[a, j, i]
            }
          }
        }
      }
      else if (string == 'F2') {
        k.i.j.a <- s.i.j.a <- g.i.j.a <- p.i.j.a <- f.i.j.a <- F5.i.j.a <- array(NA, dim = c(Age.mature, m, m))
        s.i.a <- S
        f.i.a <- matrix(NA, nrow = m, ncol = Age.mature)
        for (a in 1:Age.mature) {
          for (i in 1:m) {
            f.i.a[i, a] <- F1[i, a]*(mat[i, a]*mod)*F3[i, a]*F4[a]
            for (j in 1:m) {
              p.i.j.a[a, j, i] <- s.i.a[i, a]*G[a, j, i]
              f.i.j.a[a, j, i] <- f.i.a[i, a]*F5[j]
              k.i.j.a[a, j, i] <- p.i.j.a[a, j, i]+f.i.j.a[a, j, i]
            }
          }
        }
      }
      else if (string == 'F3') {
        k.i.j.a <- s.i.j.a <- g.i.j.a <- p.i.j.a <- f.i.j.a <- F5.i.j.a <- array(NA, dim = c(Age.mature, m, m))
        s.i.a <- S
        f.i.a <- matrix(NA, nrow = m, ncol = Age.mature)
        for (a in 1:Age.mature) {
          for (i in 1:m) {
            f.i.a[i, a] <- F1[i, a]*F2[i, a]*(mat[i, a]*mod)*F4[a]
            for (j in 1:m) {
              p.i.j.a[a, j, i] <- s.i.a[i, a]*G[a, j, i]
              f.i.j.a[a, j, i] <- f.i.a[i, a]*F5[j]
              k.i.j.a[a, j, i] <- p.i.j.a[a, j, i]+f.i.j.a[a, j, i]
            }
          }
        }
      }
      else if (string == 'F4') {
        k.i.j.a <- s.i.j.a <- g.i.j.a <- p.i.j.a <- f.i.j.a <- F5.i.j.a <- array(NA, dim = c(Age.mature, m, m))
        s.i.a <- S
        f.i.a <- matrix(NA, nrow = m, ncol = Age.mature)
        for (a in 1:Age.mature) {
          for (i in 1:m) {
            f.i.a[i, a] <- F1[i, a]*F2[i, a]*F3[i, a]*(mat[a]*mod)
            for (j in 1:m) {
              p.i.j.a[a, j, i] <- s.i.a[i, a]*G[a, j, i]
              f.i.j.a[a, j, i] <- f.i.a[i, a]*F5[j]
              k.i.j.a[a, j, i] <- p.i.j.a[a, j, i]+f.i.j.a[a, j, i]
            }
          }
        }
      }
      else if (string == 'F5') {
        k.i.j.a <- s.i.j.a <- g.i.j.a <- p.i.j.a <- f.i.j.a <- F5.i.j.a <- array(NA, dim = c(Age.mature, m, m))
        s.i.a <- S
        f.i.a <- matrix(NA, nrow = m, ncol = Age.mature)
        for (a in 1:Age.mature) {
          for (i in 1:m) {
            f.i.a[i, a] <- F1[i, a]*F2[i, a]*F3[i, a]*F4[a]
            for (j in 1:m) {
              p.i.j.a[a, j, i] <- s.i.a[i, a]*G[a, j, i]
              f.i.j.a[a, j, i] <- f.i.a[i, a]*(mat[j]*mod)
              k.i.j.a[a, j, i] <- p.i.j.a[a, j, i]+f.i.j.a[a, j, i]
            }
          }
        }
      }
      else if (string == 'F') {
        k.i.j.a <- s.i.j.a <- g.i.j.a <- p.i.j.a <- f.i.j.a <- F5.i.j.a <- array(NA, dim = c(Age.mature, m, m))
        s.i.a <- S
        f.i.a <- matrix(NA, nrow = m, ncol = Age.mature)
        for (a in 1:Age.mature) {
          for (i in 1:m) {
            for (j in 1:m) {
              p.i.j.a[a, j, i] <- s.i.a[i, a]*G[a, j, i]
              k.i.j.a[a, j, i] <- p.i.j.a[a, j, i]+(mat[a, j, i]*mod)
            }
          }
        }
      }
      else {
        return('wrong input: array')
        break
      }
      # lambda
      {
        lam.list <- c()
        n.0 <- which(s1.db$Age == 0)
        n.0.h <- log(s1.db[n.0,]$h2)
        n.0.v <- hist(n.0.h, breaks = e.pred, plot = FALSE)$counts
        init.n.a.v <- n.0.v+c*F5
        n.list <- c(sum(init.n.a.v))
        for (a in 1:Age.mature) {
          n.a.v <- k.i.j.a[a,,]%*%init.n.a.v
          n.a.v <- n.a.v+c*F5
          lam.a <- sum(n.a.v)/sum(init.n.a.v) # add c individuals
          n.list <- c(n.list, sum(n.a.v))
          lam.list <- c(lam.list, lam.a)
          init.n.a.v <- n.a.v			# add size structure
        }
      }
      return(lam.list)
    }
  }
  # sensitivity analysis
  {
    F <- f.i.j.a
    mod <- c(0.999,0.9999,0.99999,0.999999,1,1.000001,1.00001,1.0001,1.001)
    vrname <- list('S', 'G', 'F1', 'F2', 'F3', 'F4', 'F5', 'F')
    vrarray <- list(S, G, F1, F2, F3, F4, F5, F)
    for (i in (1:8)) { # 8 is the length of vrname & vrarray
      sens.i <- matrix(nrow = m, ncol = length(mod))
      for (j in 1:9) { # 9 is the length of mod
        sens.i[,j] <- sensitivity(vrarray[[i]], vrname[[i]], mod[j])
        assign(paste0('sens', vrname[i]), sens.i)
      }
    }
    sens.plot <- matrix(ncol = 3, nrow = length(mod)*length(vrname), NA)
    sens.l <- list(sensS, sensG, sensF1, sensF2, sensF3, sensF4, sensF5, sensF)
    for (i in 1:length(sens.l)) {
      for (j in 1:ncol(sens.l[[i]])) {
        sens.plot[((i-1)*length(mod)+1):(i*length(mod)),1] <- mod
        sens.plot[((i-1)*length(mod)+1):(i*length(mod)),2][j] <- mean(sens.l[[i]][,5])-mean(sens.l[[i]][,j]) # column 5 correspond to no modification
        #sens.plot[((i-1)*length(mod)+1):(i*length(mod)),3] <- rep(vrname[[i]], length(mod))
      }
    }
    sens.s.p <- as.data.frame(sens.plot[1:9,])
    colnames(sens.s.p) = c('mod','delta.lam','vr')
    plot(sens.s.p$mod, sens.s.p$delta.lam, type = 'l', col = 'black')
    sens.s.g <- as.data.frame(sens.plot[10:18,])
    colnames(sens.s.g) = c('mod','delta.lam','vr')
    lines(sens.s.g$mod, sens.s.g$delta.lam, col = 'red')
    sens.s.f <- as.data.frame(sens.plot[64:72,])
    colnames(sens.s.f) = c('mod','delta.lam','vr')
    lines(sens.s.f$mod, sens.s.f$delta.lam, col = 'blue')
  }
}
# size structure & population size change over time
{
  pop.size.NM <- c()
  pop.size.WM <- c()
  size.v.a.NM.s <- list(NA) # size structure vector by age, NO MIGRATION, standardized
  size.v.a.WM.s <- list(NA) # size structure vector by age, WITH MIGRATION, standardized
  for (a in 1:(Age.mature+1)) {
    pop.size.NM <- c(pop.size.NM, sum(size.v.a.NM[[a]]))			# get population size
    size.v.a.NM.s[[a]] <- size.v.a.NM[[a]]/sum(size.v.a.NM[[a]])		# standardized sized vector, NO MIGRATION
    pop.size.WM <- c(pop.size.WM, sum(size.v.a.WM[[a]]))
    size.v.a.WM.s[[a]] <- size.v.a.WM[[a]]/sum(size.v.a.WM[[a]])
  }
  obs.N.total <- data.frame(N = c(NA,N.list), Age = c(0,Age.pred))
  ps.df <- data.frame(Age = rep(c(0,Age.pred), 2), N = c(pop.size.NM, pop.size.WM), migration = c(rep("no",Age.mature+1), rep("yes", Age.mature+1)))
  ps.plot <- ggplot(ps.df, aes(x = Age, y = N)) +
    theme_minimal() +
    xlab("Succesional age (years)") +
    ylab("Projected population size") +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15), legend.text = element_text(size = 12), legend.title = element_text(size = 15), legend.position="right") +
    geom_line(aes(color = migration), size = 1) +
    geom_point(data = obs.N.total, aes(x = Age, y = N))
    ps.plot
    #ggsave("porjected-N.pdf", ps.plot, device = "pdf", width = 9, height = 6, units = "in", dpi = 180*2)
  size.str.mat.NM.s <- matrix(unlist(size.v.a.NM.s), ncol = 100, byrow = TRUE) #
  size.str.mat.WM.s <- matrix(unlist(size.v.a.WM.s), ncol = 100, byrow = TRUE)
  zlim <- max(max(size.str.mat.NM.s), max(size.str.mat.WM.s))
  #pdf(file="size.str.NM.pdf",width=8,height=8)
  hist3D(y = exp(x.pred), x = c(0,Age.pred), z = size.str.mat.NM.s, col = "grey", border = "black", xlab = "Successional age (years)", ylab = "Size (m)", zlab = "Density", zlim = c(0, zlim)
    #,theta = -90
    ,ticktype = "detailed"
    )
  #dev.off()
  #pdf(file="size.str.WM.pdf",width=8,height=8)
  hist3D(y = exp(x.pred), x = c(0,Age.pred), z = size.str.mat.WM.s, col = "grey", border = "black", xlab = "Successional age (years)", ylab = "Size (m)", zlab = "Density", zlim = c(0, zlim)
    #,theta = -90
    ,ticktype = "detailed"
    )
  #dev.off()
}
# plotly saving images
{
  widget_file_size <- function(p) {
  d <- tempdir()
  withr::with_dir(d, htmlwidgets::saveWidget(p, "index.html"))
  f <- file.path(d, "index.html")
  mb <- round(file.info(f)$size / 1e6, 3)
  message("File is: ", mb," MB")
}
}

plotsrdt <- load("./Data/r_data/plots.RData")
s.plot <- get(plotsrdt[1])
g.mean.plot <- plotsrdt[2]
f1.plot <- plotsrdt[3]
f2.plot <- plotsrdt[4]
f3.plot <- plotsrdt[5]
f4.plot <- plotsrdt[6]
f6.plot <- plotsrdt[7]



 # Multiplot function
 {
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                      ncol = cols, nrow = ceiling(numPlots/cols))
    }

   if (numPlots==1) {
      print(plots[[1]])

    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }
 }
 multiplot(s.plot, f1.plot, f2.plot, f3.plot, f4.plot, g.mean.plot, f6.plot, cols = 2)
