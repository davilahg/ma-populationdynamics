# Population dynamics of Mimosa acantholoba
#
setwd("/Users/Edgar/Google\ Drive/Edgar/Trabajo/FC/Dirección\ de\ tesis/L\ 2017\ Gerardo\ Dávila/")
# load requiered libraries
library(AICcmodavg)
library(betareg)
library(EnvStats)
library(gamm4)
library(ggplot2)
library(lme4)
library(mgcv)
library(reshape2)
library(plotly)
library(gridExtra)
library(plot3D)
library(RCurl)
library(fields)
# read data
{
  s.db <- read.csv("./Data/survival.csv") # survival with estimated h1
  s1.db <- read.csv("./Data/survival-1.csv") # survival without estimated h1
  g.db <- read.csv("./Data/growth.csv")
  f1.db <- read.csv("./Data/reproduction-probability.csv")
  f2.db <- read.csv("./Data/reproduction-fruits.csv")
  f3.db <- read.csv("./Data/reproduction-seeds.csv")
  f4.db <- read.csv("./Data/establishment.csv")
  f5.db <- read.csv("./Data/understory.csv")
}
# definitions
{
 Age.mature <- 100
 m <- 100 # small-matrix dimension
 M <- 200000 # big-matrix dimension
 min.lh1 <- min(s.db$ln.h1, na.rm = TRUE)-0.0001
 min.lh2 <- min(s.db$ln.h1, na.rm = TRUE)-0.0001
 max.lh1 <- max(s.db$ln.h1, na.rm = TRUE)+0.0001
 max.lh2 <- max(s.db$ln.h2, na.rm = TRUE)+0.0001
 E.pred <- seq(min(min.lh1, min.lh2), max(max.lh1, max.lh2), length.out = M+1)
 h <- E.pred[2]-E.pred[1]
 X.pred <- (E.pred[2:(M+1)]+E.pred[1:M])/2
 t1.pred <- rep(0, length.out = M)
 sM <- M/m
 e.pred <- seq(min(X.pred), max(X.pred), length.out = m+1)
 x.pred <- (e.pred[2:(m+1)]+e.pred[1:m])/2
 f4.pred <- Age.pred <- 1:Age.mature
 nplot <- nlevels(s.db$PLOT)
 ex.pr <- exp(x.pred)
}
# modelling vital rates
{
 glmm.s <- glmer(sup~ln.h1*Age+(0+ln.h1| Census)+(0+Age+ln.h1|plot), s.db, binomial) # survival: s
 lmm.g <- lmer(ln.h2~ln.h1*Age+(1+ln.h1| Census)+(1+ln.h1|plot), g.db) # growth: g
 glmm.f1 <- glmer(Rep~Age*ln.h1 + (1| plot), f1.db, binomial) # reproduction probability: f1
 gam.f2 <- gamm4(TotFrut~t2(ln.h1, Age, k = 4), random = ~ (0 + Age + ln.h1|plot), family = negbin(1.179556), data = f2.db) # fruit number per individual: f2
 f3.db$Plot <- as.factor(as.character(f3.db$Plot))
 glmm.f3 <- glmer(N.seed~ln.h1+Age + (1| Plot), f3.db, poisson) # seed number per fruit: f3
 den.f5 <- density(x = f5.db$h1, n = m, na.rm = TRUE, from = min(f5.db$h1, na.rm = TRUE), to = max(exp(max.lh1), exp(max.lh2)))
 den.f5$y[74:100] <- 0 # no probability of individuals higher than 3.2 m (maximum observed height of understory)
 den.f5$y <- den.f5$y/sum(den.f5$y)
 # gam.f4
 {
	und.ages <- f4.db$Age <- as.factor(as.character(f4.db$Age)) # get ages with register in understory data
	stb.n <- as.data.frame(table(f4.db$Age)) # get number of recruits per year
	names(stb.n) <- c("Age", "SNB") # SNB = recruits
	und.ages <- as.numeric(levels(und.ages))
	new.stb.n <- data.frame(Age = und.ages[order(und.ages)], SNB = rep(NA, length(und.ages))) # create new dataframe
	for (i in 1:nrow(stb.n)) { # order dataframe
		if (stb.n$Age[i] == new.stb.n$Age[which(new.stb.n$Age == stb.n$Age[i])]) {
			new.stb.n$SNB[which(new.stb.n$Age == stb.n$Age[i])] <- stb.n$SNB[i]
		}
	}
	new.stb.n$SNB[which(is.na(new.stb.n$SNB) == TRUE)] <- 0 # setting 0 to ages with register but no recruits
	stb.n <- new.stb.n
	stb.n <- transform(stb.n, Age = as.numeric(as.character(Age)))
	n.ages <- stb.n$Age
	f4.a <- rep(NA, length(n.ages)) # create new dataframe for establishment probability
	for (a in n.ages) {
		stb.a <- subset(s.db, Age == a) # create new data for each year
		n.str.a <- as.data.frame(cbind(x.pred, hist(log(as.numeric(as.character(stb.a$h1))), breaks = e.pred, plot = FALSE)$counts))
		names(n.str.a) <- c("ln.h1", "n")
		n.vec.a <- n.str.a$ln.h1
		pred.f1.a <- predict(glmm.f1, newdata = data.frame(ln.h1 = x.pred, Age = rep(a-1, length(x.pred))), type = "response", re.form = NA)
		pred.f2.a <- predict(gam.f2$gam, newdata = data.frame(ln.h1 = x.pred, Age = rep(a-1, length(x.pred))), type = "response", re.form = NA)
		pred.f3.a <- predict(glmm.f3, newdata = data.frame(ln.h1 = x.pred, Age = rep(a-1, length(x.pred))), type = "response", re.form = NA)
		pred.f.a <- pred.f1.a*pred.f2.a*pred.f3.a # predicting and multiplying reproduction functions
		new.n.a <- sum(pred.f.a) # seed number produced in t-1
		stb.n.a <- stb.n$SNB[which(stb.n$Age == a)] # observed recruits number
		f4.a[which(stb.n$Age == a)] <- stb.n.a/new.n.a # establishment probability from t-1 to t
		}
	f4.a[which(f4.a == 0)] <- 0.00000001 # change 0 to perform beta regression
	x <- n.ages
	gam.f4 <- gam(f4.a~s(x, k = 3), family = gaussian(log))
 }
}
# discretizing vital rates functions
{
 # g array
 GGG <- matrix(NA, ncol = Age.mature, nrow = M)
 for (a in 1:Age.mature)
  GGG[, a] <- predict(lmm.g, newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), type = "response", re.form = NA)
 MUG <- matrix(NA, ncol = Age.mature, nrow = m)
 for (a in 1:Age.mature) {
  init.g <- 1
  for (i in 1:m) {
     MUG[i, a] <- mean(GGG[init.g:(init.g+sM-1), a])
   init.g <- init.g + sM
  }
 }
 sdg <- sigma(lmm.g)
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
  SSS[ ,a] <- predict(glmm.s, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), re.form = NA)
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
  FF1[ ,a] <- predict(glmm.f1, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), re.form = NA)
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
  FF2[ ,a] <- predict(gam.f2$gam, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), re.form = NA)
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
  FF3[ ,a] <- predict(glmm.f3, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), re.form = NA)
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
 F4 <- predict(gam.f4, newdata = data.frame(x = f4.pred), type = "response")
 # f5 vector
 F5 <- den.f5$y
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
# kernel plots
{
  k.ages <- c(1,10,20,30,40,50,60,80,100)
  zlim.k = c(0, max(k.i.j.a[,,]))
  #pdf(file="kernel.a.pdf",width=8,height=8)
  par(mfrow=c(3,3), tcl=-0.5, family="serif", mai=c(0.3,0.3,0.3,0.3))
  for (a in k.ages) {
  	if ( a != 100) {
  		image(ex.pr, ex.pr, k.i.j.a[a,,], zlim = zlim.k, xlab = "", ylab = "", main = paste0("Successional age =  ", a))
  	} else {
  		image(ex.pr, ex.pr, k.i.j.a[a,,], zlim = zlim.k, xlab = "", ylab = "", main = "Mature forest")
  	 }
  	}
  mtext(substitute(paste("Size ", italic(t))), side=1, outer=T, at=0.5)
  mtext(substitute(paste("Size ", italic(t), " + 1")), side=2, outer=T, at=0.5)
  #image.plot(legend.only=TRUE, zlim= zlim.k, col =  heat.colors(12),horizontal = F)
  #dev.off()
  # growth plots
  zlim.g = c(0, max(G[,,]))
  pdf(file="G.a.pdf",width=8,height=8)
  par(mfrow=c(3,3), tcl=-0.5, family="serif", mai=c(0.3,0.3,0.3,0.3))
  for (a in k.ages) {
  	if (a != 100) {
  		image(ex.pr, ex.pr, G[a,,], zlim = zlim.g, xlab = "", ylab = "", main = paste0("Successional age =  ", a))
  	 } else {
  		image(ex.pr, ex.pr, G[a,,], zlim = zlim.g, xlab = "", ylab = "", main = "Mature forest")
  	 }
  	}
  mtext(substitute(paste("Size ", italic(t))), side=1, outer=T, at=0.5)
  mtext(substitute(paste("Size ", italic(t), " + 1")), side=2, outer=T, at=0.5)
  image.plot(legend.only=TRUE, zlim= zlim.g, col =  heat.colors(12),horizontal = F)
  dev.off()
  # p(x,y,t) = g(x,y,t)*s(x,t)
  zlim.p = c(0, max(p.i.j.a[,,]))
  #pdf(file="P.a.pdf",width=8,height=8)
  par(mfrow=c(3,3), tcl=-0.5, family="serif", mai=c(0.3,0.3,0.3,0.3))
  for (a in k.ages) {
  	if (a != 100) {
  		image(ex.pr, ex.pr, p.i.j.a[a,,], zlim = zlim.p, xlab = "", ylab = "", main = paste0("Successional age =  ", a))
  	} else {
  		image(ex.pr, ex.pr, p.i.j.a[a,,], zlim = zlim.p, xlab = "", ylab = "", main = "Mature forest")
  	 }
  	}
  mtext(substitute(paste("Size ", italic(t))), side=1, outer=T, at=0.5)
  mtext(substitute(paste("Size ", italic(t), " + 1")), side=2, outer=T, at=0.5)
  #image.plot(legend.only=TRUE, zlim= zlim.p, col =  heat.colors(12),horizontal = F)
  #dev.off()
  # fecundity
  zlim.f = c(0, max(f.i.a[,]))
  #pdf(file="F.a.pdf",width=4,height=4)
  image(ex.pr, 1:Age.mature, f.i.a, zlim = zlim.p, xlab = "", ylab = "", main = "")
  #image.plot(legend.only=TRUE, zlim= zlim.f, col =  heat.colors(12),horizontal = F)
  #dev.off()
}
# total predicted lambda
{
	lam.list <- c()
	n.0 <- which(s1.db$Age == 0) # get trees in first year (row number)
	n.0.h <- log(s1.db[n.0,]$h2) # get hight ... since this dataframe has a new estimated h1, h2 is the observed first height for the first year
	n.0.v <- hist(n.0.h, breaks = e.pred, plot = FALSE)$counts # count number of trees in each size class
	init.n.a.v <- n.0.v # rename vector
	size.v.a.NM <- list(init.n.a.v) # create list for size structure change, NM = no migration
	n.list <- c(sum(init.n.a.v)) # create population size vector & setting first value
	for (a in 1:Age.mature) {
       		n.a.v <- k.i.j.a[a,,]%*%init.n.a.v
       		lam.a <- sum(n.a.v)/sum(init.n.a.v) # add c individuals
       		size.v.a.NM[[a+1]] <- n.a.v
        	n.list <- c(n.list, sum(n.a.v))
       		lam.list <- c(lam.list, lam.a)
       		init.n.a.v <- n.a.v			# add size structure
		}
	tot.lam.pred <- lam.list.NM <- lam.list # lam.list is the transitory lambda vector
	prd.N.total <- data.frame(N = n.list, Age = c(0,Age.pred))
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
     lambda.df <- as.data.frame(list(lambda = lam.list, Age = 1:Age.mature)) # create lambda data frame for plotting
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
  distance.N <- function(c) {
  	  dist.N <- c()
	  n.0 <- which(s1.db$Age == 0)	#
     	  n.0.h <- log(s1.db[n.0,]$h2)	#
     	  n.0.v <- hist(n.0.h, breaks = e.pred, plot = FALSE)$counts
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
		dist.N.p <- sqrt(dist.N.p)
		dist.N <- c(dist.N, dist.N.p)
		}
	        dist.N <- sqrt(sum(dist.N))
	        cat(paste0("c = ", c, ", dist = ", dist.N, "\n"))
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
    opt.brent$par # c = 9.35491439069071, dist = 81.8803591350505
    opt.lbf <- optim(0, distance.N, method = "L-BFGS-B")
    c <- opt.lbf$par # c = 9.35391449733551, dist = 81.8803596331657
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
       	init.n.a.v <- n.0.v+c*F5
  	size.v.a.WM <- list(init.n.a.v) # size structure vector by year, WM = with migration
       	n.list <- c(sum(init.n.a.v))
       	for (a in 1:Age.mature) {
         		n.a.v <- k.i.j.a[a,,]%*%init.n.a.v
         		n.a.v <- n.a.v+c*F5
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
# vital rates plots
{
 # definitions
 {
	 bin.size <- 0.1
	 n.sim <- 1000
	 h1.pred <- seq(min(s1.db$h1, na.rm = TRUE), ceiling(max(s1.db$h1, na.rm = TRUE)), length.out = m)
	 h2.pred <- seq(min(s1.db$h2, na.rm = TRUE), ceiling(max(s1.db$h2, na.rm = TRUE)), length.out = m+1)
	 mids.h2 <- hist(s1.db$h2, breaks = h2.pred, plot = FALSE)$mids
	 ln.mids.h2 <- log(mids.h2)
	 lh1.pred <- log(h1.pred)
	 data.s.j <- subset(s1.db, !is.na(h1) & !is.na(sup))
	 data.s.j <- droplevels(data.s.j)
	 data.g.j <- subset(data.s.j, !is.na(h2))
	 data.g.j <- droplevels(data.g.j)
 	 h1.pred.1 <- seq(min(s1.db$h1, na.rm = TRUE), ceiling(max(s1.db$h1, na.rm = TRUE)), length.out = m+1)
	 mids <- hist(s1.db$h1, breaks = h1.pred.1, plot = FALSE)$mids
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
	 	s.pred <- as.vector(predict(glmm.s, newdata = data.frame(ln.h1 = ln.mids, Age = rep(a, length(mids))), type = "response", re.form = NA))
	 	pred.s.j[(a-1)*length(mids)+1:length(mids), ] <- cbind(rep(a, length(mids)), mids, s.pred)
	 }
	 g.s.j <- ggplot(pred.s.j) +
	 theme_minimal() +
	 theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
	 xlab("Successional age (years)") +
	 ylab("Height (m)") +
	 scale_x_continuous(expand = c(0,0)) +
	 scale_y_continuous(expand = c(0,0)) +
	 scale_fill_gradientn(colours = hcl.colors(12, "YlOrRd"), guide = "none") +
	 scale_alpha(guide = "none") +
	 geom_raster(aes(age.pred, h1.pred, fill = s.pred, alpha = 0.35)) +
	 scale_colour_gradient(low = "#d90000", high = "#fcdf03") +
	 geom_point(data = plot.data.s.j, aes(age, h1, color = prop, size = n), alpha = 0.5) +
	 guides(colour = guide_colourbar(title = "Survival\nprobability", order = 1), size = guide_legend(title = "Sample\nsize", order = 2))
	 dev.new(width = 10, height = 8)
	 g.s.j
	 #ggsave("survival.pdf", g.s.j, device = "pdf", width = 9, height = 6, units = "in", dpi = 180*2)
 }
 # growth
 {
	# create matrix with observed data
 	 h1.pred.g <- seq(min(g.db$h1, na.rm = TRUE), ceiling(max(g.db$h1, na.rm = TRUE)), length.out = m)
 	 h1.pred.1.g <- seq(min(g.db$h1, na.rm = TRUE), ceiling(max(g.db$h1, na.rm = TRUE)), length.out = m+1)
	 data.g.mean <- subset(g.db, !is.na(h1) & !is.na(h2))
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
	 	g.m <- predict(lmm.g, newdata = data.frame(ln.h1 = x.pred, Age = rep(a, m)), type = "response", re.form = NA)
      	 	delta.g <- exp(g.m)-exp(x.pred)
	 	pred.g.mean[(a-1)*length(mids)+1:length(mids), ] <- cbind(rep(Age.pred[a], length(mids)), mids, delta.g)
	 }
	 # plot
	 g.mean.plot <- ggplot(pred.g.mean) +
	 theme_minimal() +
	 theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
	 xlab("Successional age (years)")+
	 ylab("Height (m)") +
	 scale_x_continuous(expand = c(0,0)) +
	 scale_y_continuous(expand = c(0,0)) +
	 scale_fill_gradientn(colours = hcl.colors(12, "YlOrRd"), guide = "none") +
	 scale_alpha(guide = "none") +
	 geom_raster(aes(age.pred, h1.pred, fill = g.mean, alpha = 0.35)) +
	 scale_colour_gradient(low = "#d90000", high = "#fcdf03") +
	 geom_point(data = na.omit(data.g.df), aes(age, h1, color = gmean, size = n), alpha = 0.5) +
	 guides(colour = guide_colourbar(title = "Average\ngrowth (m)", order = 1), size = guide_legend(title = "Sample\nsize", order = 2))
	 dev.new(width = 10, height = 8)
	 g.mean.plot
	 #ggsave("mean-growth.pdf", g.mean.plot, device = "pdf", width = 9, height = 6, units = "in", dpi = 180*2)
	 }
 }
 # Fecundity 1
 {
 	 data.rep <- subset(f1.db, !is.na(h1) & !is.na(sup))
	 data.rep <- droplevels(data.rep)
	 data.rep1 <- subset(f1.db, Rep == 1)
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
	 plot.data.rep <- subset(plot.data.rep, !is.na(prop.f1))
	 pred.f1.j <- as.data.frame(matrix(ncol = 3, nrow = length(Age.pred)*length(mids)))
	 names(pred.f1.j) <- c("age.pred", "h1.pred", "f1.pred")
	 for (a in Age.pred) {
	  pred.f1 <- as.vector(predict(glmm.f1, newdata = data.frame(ln.h1 = ln.mids, Age = rep(a, m)), type = "response", re.form = NA))
	  pred.f1.j[(a-1)*m+1:m, ] <- cbind(rep(a, m), mids, pred.f1)
	  }
	 g.f1.j <- ggplot(pred.f1.j) +
	 theme_minimal() +
	 theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
	 xlab("Successional age (years)") +
	 ylab("Height (m)") +
	 scale_x_continuous(expand = c(0,0)) +
	 scale_y_continuous(expand = c(0,0)) +
	 scale_fill_gradientn(colours = hcl.colors(12, "YlOrRd"), guide = "none") +
	 scale_alpha(guide = "none") +
	 geom_raster(aes(age.pred, h1.pred, fill = f1.pred, alpha = 0.35)) +
	 scale_colour_gradient(low = "#d90000", high = "#fcdf03") +
	 geom_point(data = plot.data.rep, aes(age, h1, color = prop.f1, size = n), alpha = 0.5) +
	 guides(colour = guide_colourbar(title = "Reproduction\nprobability", order = 1), size = guide_legend(title = "Sample\nsize", order = 2))
	 dev.new(width = 10, height = 8)
	 g.f1.j
	 #ggsave("f1.pdf", g.f1.j, device = "pdf", width = 9, height = 6, units = "in", dpi = 180*2)
 }

 # Fecundity 2
 {
	 plot.data.f2 <- as.data.frame(matrix(ncol = 3, nrow = m*Age.mature))
	 names(plot.data.f2) <- c("age.pred.f2", "h1.pred.f2", "frt.n.pred")
	 for (a in Age.pred) {
	     frt.n.pred <- as.vector(predict(gam.f2$gam, newdata = data.frame(ln.h1 = x.pred, Age = rep(a, m)), type = "response", re.form = NA))
	     plot.data.f2[(a-1)*m+1:m, ] <- cbind(rep(a, m), mids, frt.n.pred)
	     }
	 real.rep.data <- subset(f2.db, Rep == 1)
	 my_varsf2 <- c("Age", "h1", "TotFrut")
	 real.rep.data <- real.rep.data[my_varsf2]
	 g.f2.j <- ggplot(plot.data.f2) +
		 theme_minimal() +
		 theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
		 xlab("Successional age (years)") +
		 ylab("Height (m)") +
		 scale_x_continuous(expand = c(0,0)) +
		 scale_y_continuous(expand = c(0,0)) +
		 scale_fill_gradientn(colours = hcl.colors(12, "YlOrRd"), guide = "none") +
		 scale_alpha(guide = "none") +
		 geom_raster(aes(age.pred.f2, h1.pred.f2, fill = frt.n.pred, alpha = 0.35)) +
		 scale_colour_gradient(low = "#d90000", high = "#fcdf03") +
		 geom_point(data = real.rep.data, aes(Age, h1, color = TotFrut), size = 3, alpha = 0.5) +
		 dev.new(width = 10, height = 8)
		 g.f2.j
	 #ggsave("f2.pdf", g.f2.j, device = "pdf", width = 9, height = 6, units = "in", dpi = 180*2)
 }
 # Fecundity 3
 {
	 plot.data.f3 <- as.data.frame(matrix(ncol = 3, nrow = m*Age.mature))
	 names(plot.data.f3) <- c("Age.pred", "h1.pred.f3", "sd.n.pred")
	 for (a in 1:length(Age.pred)) {
	     sd.n.pred <- as.vector(predict(glmm.f3, newdata = data.frame(ln.h1 = ln.mids, Age = rep(a, m)), type = "response", re.form = NA))
	     plot.data.f3[(a-1)*m+1:m, ] <- cbind(rep(a, length(mids)), mids, sd.n.pred)
	     }
	 real.rep.data <- f3.db
	 my_varsf3 <- c("Age", "h1", "N.seed")
	 real.rep.data <- real.rep.data[my_varsf3]
	 g.f3.j <- ggplot(plot.data.f3) +
	 theme_minimal() +
	 theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
	 xlab("Successional age (years)") +
	 ylab("Height (m)") +
	 scale_x_continuous(expand = c(0,0)) +
	 scale_y_continuous(expand = c(0,0)) +
	 scale_fill_gradientn(colours = hcl.colors(12, "YlOrRd"), guide = "none") +
	 scale_alpha(guide = "none") +
	 geom_raster(aes(Age.pred, h1.pred.f3, fill = sd.n.pred, alpha = 0.35)) +
	 scale_colour_gradient(low = "#d90000", high = "#fcdf03") +
	 geom_point(data = real.rep.data, aes(Age, h1, color = N.seed), size = 3, alpha = 0.5) +
	 dev.new(width = 10, height = 8)
	 g.f3.j
	 #ggsave("f3.pdf", g.f3.j, device = "pdf", width = 9, height = 6, units = "in", dpi = 180*2)
 }
 # Fecundity 4
 {
	 f4.line <- data.frame(Age = Age.pred, f4 = F4)
	 f4.dots <- data.frame(Age = n.ages, f4 = f4.a)
	 f4.p.a <- ggplot(f4.line, aes(x = Age, y = f4)) +
		theme_minimal() +
	 	theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
		xlab("Succesional age (years)") +
		ylab("Probability of establishment") +
		geom_line(size = 1) +
	 	geom_point(data = f4.dots, aes(x = Age, y = f4))
	 f4.p.a
	 #ggsave("f4.pdf", f4.p.a, device = "pdf", width = 9, height = 6, units = "in", dpi = 180*2)
 }
 # Fecundity 5
 {
  	f5.p.df <- data.frame(f5.db$Age, f5.db$h1)
  	f5.p.gg <- ggplot(f5.p.df, aes(x=f5.db.h1)) +
	    theme_minimal() +
	    geom_histogram(aes(y=..density..), colour="black", fill="white")+
	    geom_density(alpha=.2, fill="#FF6666") +
	    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
	    xlab(expression(paste("Recruits height (m)"))) +
	    ylab(expression("Density"))
	f5.p.gg
    	#ggsave("f5.pdf", f5.p.gg, device = "pdf", width = 9, height = 6, units = "in", dpi = 180*2)
 }
