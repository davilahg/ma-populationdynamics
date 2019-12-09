# Population dynamics of Mimosa acantholoba
# open r terminal in desired work directory or substitute dir by directory path
dir <- getwd()
setwd(dir) # setting work directory
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
}
# modelling vital rates
{
 glmm.s <- glmer(sup~ln.h1*Age+(0+ln.h1| Census)+(0+Age+ln.h1|plot), s.db, binomial) # survival: s
 lmm.g <- lmer(ln.h2~ln.h1*Age+(1+ln.h1| Census)+(1+ln.h1|plot), g.db) # growth: g
 glmm.f1 <- glmer(Rep~Age*ln.h1 + (1| plot), f1.db, binomial) # reproduction probability: f1
 gam.f2 <- gamm4(TotFrut~t2(ln.h1, Age, k = 4), random = ~ (0 + Age + ln.h1|plot), family = negbin(1.179556), data = f2.db) # fruit number per individual: f2
 glmm.f3 <- glmer(N.seed~ln.h1+Age + (1| plot), f3.db, poisson) # seed number per fruit: f3
 den.f5 <- density(x = f5.db$h1, n = m, na.rm = TRUE, from = min(f5.db$h1, na.rm = TRUE), to = max(exp(max.lh1), exp(max.lh2)))
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
}
# creating kernel
{
 k.i.j.a <- s.i.j.a <- g.i.j.a <- p.i.j.a <- f.i.j.a <- F5.i.j.a <- array(NA, dim = c(Age.mature, m, m)) 
 s.i.a <- S
 f.i.a <- matrix(NA, nrow = m, ncol = Age.mature)
 for (a in 1:Age.mature)
  f.i.a[, a] <- F1[, a]*F2[, a]*F3[, a]*F4[a]
 for (a in 1:Age.mature)
  for (i in 1:m) {
   s.i.j.a[a, i, ] <- s.i.a[i, a]
   g.i.j.a[a, , ] <- G[a, , ]
   p.i.j.a[a, i, ] <- p.i <- s.i.j.a[a, i, ]*G[a, i, ]
   f.i.j.a[a, i, ] <- f.i.a[i, a]
   F5.i.j.a[a, i, ] <- F5
   f.i <- f.i.a[i, a]*F5
   f.i.j.a[a, i, ] <- f.i.j.a[a, i, ]*F5
   k.i.j.a[a, i, ] <- p.i+f.i
  }
}
# total predicted lambda
{
	lam.list <- c()
	n.0 <- which(s.db$Age == 0)
	n.0.h <- log(s.db[n.0,]$h2) # since this dataframe has a new estimated h1, h2 is the observed first height for the first year
	n.0.v <- hist(n.0.h, breaks = e.pred, plot = FALSE)$counts
	init.n.a.v <- n.0.v
	n.list <- c(sum(init.n.a.v))
	for (a in 1:Age.mature) {
		 n.a.v <- k.i.j.a[a,,]%*%init.n.a.v
		 lam.a <- sum(n.a.v)/sum(init.n.a.v)
		 init.n.a.v <- n.a.v
		 n.list <- c(n.list, sum(n.a.v))
		 lam.list <- c(lam.list, lam.a)
		 }
		tot.lam.pred <- lam.list # lam.list is the transitory lambda vector
} 
# total observed lambda
{
   ob.lam.list <- c()
 	   N.0 <- subset(s1.db, Age == 0)
 	   N.0 <- subset(N.0, sup == 1)
 	   N.init <- nrow(N.0)
 	   N.list <- c(N.init)
 	   for (i in 1:Age.mature) {
 		   ns.i <- subset(s1.db, Age == i)
 		   ns.i <- subset(ns.i, sup == 1)
 		   N.ns.i <- nrow(ns.i)
 		   if (N.ns.i != 0 & N.init != 0) {
 			   lambda <- N.ns.i/N.init
 			 } else {
 				 lambda <- NA
 			 }
 		   N.init <- N.ns.i
 		   N.list <- c(N.list, N.init)
 		   ob.lam.list <- c(ob.lam.list, lambda)
 	    }
     ns.mature <- subset(s1.db, Age == Age.mature) # avoid population overestimating due to individuals repeat at Age.mature
     ns.mature <- subset(ns.mature, sup == 1)
     ns.mature <- transform(ns.mature, Census = as.factor(as.character(Census)))
     ns.mature <- droplevels(ns.mature)
     N.mat.list <- c()
     for (i in levels(ns.mature$Census)) {
     	ns.mat.yr <- subset(ns.mature, Census == i)
     	ns.mat.N <- nrow(ns.mat.yr)
     	N.mat.list <- c(N.mat.list, ns.mat.N)
     }
     lam.m.list <- c()
     init.n.mat <- N.mat.list[1]
     for (j in 2:length(N.mat.list)) {
     	lam.j <- N.mat.list[j]/init.n.mat
     	lam.m.list <- c(lam.m.list,lam.j)
     	init.n.mat <- N.mat.list[j]
     }
     lam.mature <- geoMean(lam.m.list)
 	   ob.lam.list <- c(ob.lam.list, lam.mature)
 	   ob.lam.df <- as.data.frame(matrix(ncol = 2, nrow = Age.mature+1)) # creo base de datos para graficar
 	   names(ob.lam.df) <- c("Age", "ob.lambda")
 	   ob.lam.df$Age <- 1:(Age.mature+1)
 	   ob.lam.df$ob.lambda[1:(Age.mature+1)] <- ob.lam.list
     # create lambda data frame for plotting
     lambda.df <- as.data.frame(list(lambda = lam.list, Age = 1:Age.mature))
}
# observed lambda by plot
{   
           nplot <- nlevels(s1.db$PLOT)
           s1.db <- transform(s1.db, PLOT = as.factor(as.character(plot)))
     	   plot.lam.list <- as.data.frame(matrix(nrow = Age.mature, ncol = length(levels(s1.db$PLOT)), NA))
 	   names(plot.lam.list) <- 1:nplot
 	   for (i in levels(s1.db$PLOT)) {
  	 	 Mim.s.p <- subset(s1.db, PLOT == i)
		 if (Mim.s.p$PLOT[1] != 5) {
		  if (all(!is.na(Mim.s.p$Age))) {
		     lam.list <- c()
		     min.a <- min(Mim.s.p$Age, na.rm = TRUE)
		     max.a <- max(Mim.s.p$Age, na.rm = TRUE)
		     n.1.i <- which(Mim.s.p$Age == min.a)
		     n.1.h.i <- log(Mim.s.p[n.1.i,]$h2)[order(log(Mim.s.p[n.1.i,]$h2))]
		     init.n.a.v.i <- n.1.v.i <- hist(n.1.h.i, breaks = e.pred, plot = FALSE)$counts
		     for (a in ((min.a+1):max.a)) {
			  n.a.db <- droplevels(subset(Mim.s.p, Age == a))
			  n.a.v.i <- nlevels(n.a.db$id)
			  lam.a.i <- n.a.v.i/sum(init.n.a.v.i)
			  init.n.a.v.i <- n.a.v.i
			  lam.list <- c(lam.list, lam.a.i)
		      }
		     plot.lam.list[(min.a+1):max.a, as.numeric(i)] <- lam.list
		}
	       }
	   }
 	   plot.lam.plot <- list(NA)
 	   for (i in 1:nplot) {
      		if (i != 5) {
        		plot.l.i <- as.data.frame(matrix(nrow = length(which(!is.na(plot.lam.list[,i]) == TRUE)), ncol = 2))
        		names(plot.l.i) <- c("Age", "lambda")
       			plot.l.i$Age <- which(!is.na(plot.lam.list[,i]))
        		plot.l.i$lambda <- plot.lam.list[,i][which(!is.na(plot.lam.list[,i]) == TRUE)]
        		plot.lam.plot[[i]] <- plot.l.i 
      			}
      		else if (i == 5) {
        		plot.l.5 <- as.data.frame(matrix(nrow = 1, ncol = 2))
        		names(plot.l.5) <- c("Age", "lambda")
        		plot.l.5$Age[1] <- Age.mature
        		plot.l.5$lambda[1] <- lam.mature
        		plot.lam.plot[[5]] <- plot.l.5
      			}
 	   	       }
 }
# plots
# total, no migration
{
     lambda.df <- as.data.frame(list(lambda = lam.list, Age = 1:Age.mature))
     lambda.edad <- qplot(x = Age, y = log(lambda), data = lambda.df, color = "red", geom = "line", xlab = "Edad sucesional (años)", ylab = expression(italic(r))) +
     	scale_fill_discrete(guide=FALSE) +
     	geom_point(data = ob.lam.df, mapping = aes(x = Age, y = log(ob.lambda)))
     	lambda.edad + 
     	theme_minimal() +
     	theme(legend.position = "none")
     	lambda.edad
     #ggsave("lambda-edad.png", lambda.edad, device = "png", width = 9, height = 7, units = "in", dpi = 180*2)
}
# by plot, no migration
 {
 	 lambda.df.2 <- transform(lambda.df)
 	 p.lam.p.2 <- as.data.frame(matrix(nrow = 0, ncol = 2))
 	 names(p.lam.p.2) <- c("Age", "lambda")
 	 for (i in 1:nplot) {
 		 p.lam.p.i <-transform(plot.lam.plot[[i]])
 		 p.lam.p.2  <- rbind(p.lam.p.2, p.lam.p.i)
 	 }
 	 p.lam.ob <- ggplot(p.lam.p.2, aes(x = Age, y = lambda, alpha = 1/2), show.legend = FALSE) + 
 	 theme_minimal() +
 	 geom_line(size = 1) + 
     	 theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 15)) + 
 	 geom_line(data = lambda.df.2, aes(x = Age, y = lambda), col = "red", size = 1, alpha = 1/3, show.legend = FALSE) +
 	 labs(x = expression(paste("Abandonment time ", italic(t), " (years)")), y = expression(lambda))+
 	 scale_alpha(guide = "none")
 	 p.lam.ob
 	 #ggsave("no-migration.png", p.lam.ob, device = "png", width = 18, height = 12, units = "in", dpi = 180*2)
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

