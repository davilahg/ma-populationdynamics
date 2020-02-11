# Population dynamics of Mimosa acantholoba
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
# read data without downloaded directory
{
  s.db.text <- getURL("https://raw.githubusercontent.com/davilahg/ma-populationdynamics/master/survival.csv")
  s.db <- read.csv(text = s.db.text) # survival with estimated h1
  s1.db.text <- getURL("https://raw.githubusercontent.com/davilahg/ma-populationdynamics/master/survival-1.csv")
  s1.db <- read.csv(text = s1.db.text) # survival without estimated h1
  g.db.text <- getURL("https://raw.githubusercontent.com/davilahg/ma-populationdynamics/master/growth.csv")
  g.db <- read.csv(text = g.db.text)
  f1.db.text <- getURL("https://raw.githubusercontent.com/davilahg/ma-populationdynamics/master/reproduction-probability.csv")
  f1.db <- read.csv(text = f1.db.text)
  f2.db.text <- getURL("https://raw.githubusercontent.com/davilahg/ma-populationdynamics/master/reproduction-fruits.csv")
  f2.db <- read.csv(text = f2.db.text)
  f3.db.text <- getURL("https://raw.githubusercontent.com/davilahg/ma-populationdynamics/master/reproduction-seeds.csv")
  f3.db <- read.csv(text = f3.db.text)
  f4.db.text <- getURL("https://raw.githubusercontent.com/davilahg/ma-populationdynamics/master/establishment.csv")
  f4.db <- read.csv(text = f4.db.text)
  f5.db.text <- getURL("https://github.com/davilahg/ma-populationdynamics/blob/master/understory.csv")
  f5.db <- read.csv(text = f5.db.text)
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
	gam.f4 <- gam(f4.a~s(x, k = 3), family = betar(link = "logit"))
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
# kernel plots
k.ages <- c(1,10,20,30,40,50,60,80,100)
zlim.k = c(0, max(k.i.j.a[,,]))
#png(filename="kernel.atesis.png",width=8,height=8,units="in",res=150)
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
zlim.g = c(0, max(g.i.j.a[,,]))
#png(filename="G.a-tesis.png",width=8,height=8,units="in",res=150)
par(mfrow=c(3,3), tcl=-0.5, family="serif", mai=c(0.3,0.3,0.3,0.3))
for (a in k.ages) {
	if (a != 100) {
		image(ex.pr, ex.pr, G[a,,], zlim = zlim.k, xlab = "", ylab = "", main = paste0("Successional age =  ", a))
	 } else {
		image(ex.pr, ex.pr, G[a,,], zlim = zlim.k, xlab = "", ylab = "", main = "Mature forest")
	 } 
	}
mtext(substitute(paste("Size ", italic(t))), side=1, outer=T, at=0.5)
mtext(substitute(paste("Size ", italic(t), " + 1")), side=2, outer=T, at=0.5)
#image.plot(legend.only=TRUE, zlim= zlim.g, col =  heat.colors(12),horizontal = F) 
#dev.off()
# p(x,y,t) = g(x,y,t)*s(x,t)
zlim.p = c(0, max(p.i.j.a[,,]))
#png(filename="P.a.png",width=8,height=8,units="in",res=150)
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
#png(filename="F.a.png",width=4,height=4,units="in",res=150)
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
	size.v.a.NM <- list(NA) # create list for size structure change, NM = no migration
	n.list <- c() # create population size vector & setting first value
	for (a in 1:Age.mature) {
       		n.a.v <- k.i.j.a[a,,]%*%init.n.a.v
       		lam.a <- sum(n.a.v)/sum(init.n.a.v) # add c individuals
       		size.v.a.NM[[a]] <- n.a.v
        	n.list <- c(n.list, sum(n.a.v))
       		lam.list <- c(lam.list, lam.a)
       		init.n.a.v <- n.a.v			# add size structure
		}
	tot.lam.pred <- lam.list.NM <- lam.list # lam.list is the transitory lambda vector
	prd.N.total <- data.frame(N = n.list, Age = Age.pred)
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
	  n.0 <- which(s1.db$Age == 0)	# get row numbers with age = 0
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
	        dist.N <- mean(dist.N)
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
 # Try different estimators (for population size)
 {
    opt.brent <- optim(0, distance.N, method = "Brent", lower = 0, upper = 100)
    opt.brent$par # c = 9.35491439069071, dist = 81.8803591350505
    opt.lbf <- optim(0, distance.N, method = "L-BFGS-B")
    c <- opt.lbf$par # c = 9.35391449733551, dist = 81.8803596331657
    opt.bfgs <- optim(0, distance.N, method = "BFGS")
    opt.bfgs$par # c = 9.35491911006061, dist = 81.8803591350612
 }
# Try different estimators (for lambda)
 {
    opt.brent <- optim(0, distance.L, method = "Brent", lower = 0, upper = 100000000)
    opt.brent$par # c = 2.47619435505125, dist = 0.530052570070713
    opt.lbf <- optim(0, distance.L, method = "L-BFGS-B")
    c <- opt.lbf$par # c = 2.47519458486379, dist = 0.530052570638488
    opt.bfgs <- optim(0, distance.L, method = "BFGS")
    opt.bfgs$par # c = 15.7093246924741, dist = 0.535953969391775
 }
}

# plots
# total, no migration
{
     lambda.df.NM <- as.data.frame(list(lambda = lam.list.NM, Age = 1:Age.mature))
     lambda.edad <- qplot(x = Age, y = log(lambda), data = lambda.df.NM, color = "red", geom = "line", xlab = "Succesional age (years)", ylab = expression(italic(r))) +
     	theme_minimal() +
	theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15), legend.text = element_text(size = 12), legend.title = element_text(size = 15), legend.position = "none") + 
	geom_point(data = ob.lam.df, mapping = aes(x = Age, y = log(ob.lambda)))
     	lambda.edad
     lambda.edad
     #ggsave("lambda-edad.pdf", lambda.edad, device = "pdf", width = 9, height = 6, units = "in", dpi = 180*2)
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
	 scale_y_continuous(limits = c(0, 5.8))
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
 {
     	lam.list <- c()			# create lambda vector
     	n.0 <- which(s1.db$Age == 0)	# get row numbers with age = 0
     	n.0.h <- log(s1.db[n.0,]$h2)	#
     	n.0.v <- hist(n.0.h, breaks = e.pred, plot = FALSE)$counts
     	init.n.a.v <- n.0.v+c*F5
	size.v.a.WM <- list(NA) # size structure vector by year, WM = with migration
     	n.list <- c()
     	for (a in 1:Age.mature) {
       		n.a.v <- k.i.j.a[a,,]%*%init.n.a.v
       		n.a.v <- n.a.v+c*F5
       		lam.a <- sum(n.a.v)/sum(init.n.a.v) # add c individuals
       		size.v.a.WM[[a]] <- n.a.v
        	n.list <- c(n.list, sum(n.a.v))
       		lam.list <- c(lam.list, lam.a)
       		init.n.a.v <- n.a.v			# add size structure
       		}
            	tot.lam.pred <- lam.list
	lambda.df <- as.data.frame(list(lambda = lam.list, Age = 1:Age.mature))
   	lambda.df.2 <- transform(lambda.df, plot = as.factor("total"))

	plot.l.graf.c <- ggplot(lambda.df.2, aes(x = Age, y = lambda)) +
				theme_minimal() +
				xlab(expression(paste("Abandonment time ", italic(t), " (years)"))) +
				ylab(expression(lambda)) +
				theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 15)) + 
				geom_line(size = 2, alpha = 1/3, col = "red") + 
				geom_point(data = obs.lam.p, aes(x = Age, y = lambda), alpha = 1/2, size = 1) +
				geom_line(data = ob.lam.df, aes(x = Age, y = ob.lambda), col = "black", size = 2, alpha = 1/3, show.legend = FALSE) +
	 			scale_alpha(guide = "none") +
				scale_y_continuous(limits = c(0, 5.8))
	plot.l.graf.c
	ggsave("with-migration.pdf", plot.l.graf.c, device = "pdf", width = 9, height = 6, units = "in", dpi = 180*2)
}
# size structure change over time
{
	pop.size.NM <- c()
	pop.size.WM <- c()
	size.v.a.NM.s <- list(NA) # size structure vector by age, NO MIGRATION, standardized
	size.v.a.WM.s <- list(NA) # size structure vector by age, WITH MIGRATION, standardized
	for (a in 1:Age.mature) {
		pop.size.NM <- c(pop.size.NM, sum(size.v.a.NM[[a]]))			# get population size 
		size.v.a.NM.s[[a]] <- size.v.a.NM[[a]]/sum(size.v.a.NM[[a]])		# standardized sized vector, NO MIGRATION
		pop.size.WM <- c(pop.size.WM, sum(size.v.a.WM[[a]]))
		size.v.a.WM.s[[a]] <- size.v.a.WM[[a]]/sum(size.v.a.WM[[a]])	
		}
	ps.df <- data.frame(Age = Age.pred, N.nm = pop.size.NM, N.wm = pop.size.WM)
	ps.plot <- ggplot(ps.df, aes(x = Age) +
			theme_minimal() +
			xlab("Succesional age (years)") +
			ylab("Projected population size") +
			theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 15)) + 
			geom_line(data = ps.df, aes(x = Age, y = N.nm), col = "red", size = 2, alpha = 1/3) +
			geom_line(data = ps.df, aes(x = Age, y = N.wm), col = "blue", size = 2, alpha = 1/3) +


     	ps.nm.plot
	ps.wm.df <- data.frame(Age = Age.pred, N = pop.size.WM)
	ps.wm.plot <- qplot(xlab = "Succesional age (years)", ylab = "Projected population size") +
		theme_minimal() +
		theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15), legend.text = element_text(size = 12), legend.title = element_text(size = 15), legend.position = "none")
     	ps.nm.plot

	plot(1:Age.mature, pop.size.NM, main = "Population vector, no migration", type = "l", las = 1, bty = "l") # plot population size
	plot(1:Age.mature, pop.size.WM, main = "Population vector, with migration", type = "l", las = 1, bty = "l") 
	size.str.mat.NM.s <- matrix(unlist(size.v.a.NM.s), ncol = 100, byrow = TRUE) #    ### !! QUITAR PRIMERA ESTRUCTURA OBSERVADA
	size.str.mat.WM.s <- matrix(unlist(size.v.a.WM.s), ncol = 100, byrow = TRUE) 
	zlim <- max(max(size.str.mat.NM.s), max(size.str.mat.WM.s))
	hist3D(y = exp(x.pred), x = 1:Age.mature, z = size.str.mat.NM.s, col = "grey", border = "black", xlab = "Age", ylab = "Size", zlab = "Probability", main = "Size structure change without migration", zlim = c(0, zlim) 
	       ,theta = -90
	       )	
	hist3D(y = exp(x.pred), x = 1:Age.mature, z = size.str.mat.WM.s, col = "grey", border = "black", xlab = "Age", ylab = "Size", zlab = "Probability", main = "Size structure change with migration", zlim = c(0, zlim)
	       ,theta = -90
	       )
}

# esto ya no está limpio

# plots F5
  f5.p.df <- data.frame(f5.db$Age, f5.db$h1)
  f5.p.gg <- ggplot(f5.p.df, aes(x=f5.db.h1)) + 
    theme_minimal() +
    geom_histogram(aes(y=..density..), colour="black", fill="white")+
    geom_density(alpha=.2, fill="#FF6666") +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
    xlab(expression(paste("Recruits height (m)"))) +
    ylab(expression("Probability"))
   #ggsave("art-f5.png", f5.p.gg, device = "png", width = 10, height = 7, units = "in", dpi = 180*2)

	#plots F4
	
 	#png("f4Plot.png", width = 900, height = 600, units = "px", pointsize = 20)
 	f4.pred <- seq(0, Age.mature)
 	plot(x, f4.a, xlim = c(0, Age.mature), xlab = "Edad sucesional (años)", ylab = "Probabilidad de establecimiento")
 	#lines(f4.pred, predict(glm.f4, newdata = data.frame(x = f4.pred), type = "response"), col = "green")
    #lines(f4.pred, predict(gam.f4, newdata = data.frame(x = f4.pred), type = "response"), col = "blue")
    lines(f4.pred, predict(gam.f4, newdata = data.frame(x = f4.pred), type = "response"), col = "red")
 	#dev.off()
}

