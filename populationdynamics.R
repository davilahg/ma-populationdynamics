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
   G[a, m, i] <- 1-pnorm(e.pred[m], MUG[i, a], sdg)
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
 # f5 vector
 ff5 <- den.f5$y*h
 FF5 <- matrix(NA, nrow = m, ncol = m)
 for (i in 1:m) {
  FF5[i,] <- ff5
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






# plots F5
  f5.p.df <- data.frame(offspring$Age, offspring$h1)
  f5.p.gg <- ggplot(f5.p.df, aes(x=offspring.h1)) + 
    geom_histogram(aes(y=..density..), colour="black", fill="white")+
    geom_density(alpha=.2, fill="#FF6666") +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
    xlab(expression(paste("Altura de reclutas (m)"))) +
    ylab(expression("Probabilidad"))
   # ggsave("/Users/gdh/Documents/MEGA/Tesis/Gráficas/fc5-gg.png", f5.p.gg, device = "png", width = 10, height = 7, units = "in", dpi = 180*2)

	#plots F4
	
 	#png("f4Plot.png", width = 900, height = 600, units = "px", pointsize = 20)
 	f4.pred <- seq(0, Age.mature)
 	plot(x, f4.a, xlim = c(0, Age.mature), xlab = "Edad sucesional (años)", ylab = "Probabilidad de establecimiento")
 	#lines(f4.pred, predict(glm.f4, newdata = data.frame(x = f4.pred), type = "response"), col = "green")
    #lines(f4.pred, predict(gam.f4, newdata = data.frame(x = f4.pred), type = "response"), col = "blue")
    lines(f4.pred, predict(gam.f4, newdata = data.frame(x = f4.pred), type = "response"), col = "red")
 	#dev.off()
}

# CONSTRUCCIÓN DE MATRICES






# CONSTRUCCIÓN DE MATRICES
{
 load("Mlist.RData")
 G <- Mlist$G
 S <- MUS <- Mlist$MUS
 F1 <- MUF1 <- Mlist$MUF1
 F2 <- MUF2 <- Mlist$MUF2
 F3 <- MUF3 <- Mlist$MUF3
 F4 <- predict(gam.f4, newdata = data.frame(x = f4.pred), type = "response")
 F5 <- den.f5$y
}

# LAMBDAS (Sin migración)
{
 # Lambda total predicha
 {
     lam.list <- c()
  	 n.0 <- which(Survival$Age == 0)
     n.0.h <- log(Survival[n.0,]$h2)
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
  	  tot.lam.pred <- lam.list
 } 

 # Lambda total observdada
 {
     ob.lam.list <- c()
 	   N.0 <- subset(Survival, Age == 0)
 	   N.0 <- subset(N.0, sup == 1)
 	   N.init <- nrow(N.0)
 	   N.list <- c(N.init)
 	   for (i in 1:Age.mature) {
 		   ns.i <- subset(Survival, Age == i)
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
     ###
 	   ns.mature <- subset(Survival, Age == Age.mature) #aparto para Age.mature porque como no avanzan los años se repiten los individuos para edad 100 y se infla la población
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
 }

 # Plot: lambda total
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
 
 # Lambda por parcela (PREDICHA)
 {   
     Survival <- transform(Survival, PLOT = as.factor(as.character(plot)))
 	   nplot <- nlevels(Survival$PLOT)
     Survival <- transform(Survival, PLOT = as.factor(as.character(plot)))
 	   plot.lam.list <- as.data.frame(matrix(nrow = Age.mature, ncol = length(levels(Survival$PLOT)), NA))
 	   names(plot.lam.list) <- 1:nplot
 	   for (i in levels(Survival$PLOT)) {
  	 	 Mim.s.p <- subset(Survival, PLOT == i)
         if (Mim.s.p$PLOT[1] != 5) {
          if (all(!is.na(Mim.s.p$Age))) {
             lam.list <- c()
             min.a <- min(Mim.s.p$Age, na.rm = TRUE)
             max.a <- max(Mim.s.p$Age, na.rm = TRUE)
             n.1.i <- which(Mim.s.p$Age == min.a)
             n.1.h.i <- log(Mim.s.p[n.1.i,]$h2)[order(log(Mim.s.p[n.1.i,]$h2))]
             init.n.a.v.i <- n.1.v.i <- hist(n.1.h.i, breaks = e.pred, plot = FALSE)$counts
             for (a in (min.a:(max.a-1))) { ## CHECAR!!
                  n.a.v.i <- k.i.j.a[a+1,,]%*%init.n.a.v.i
                  lam.a.i <- sum(n.a.v.i)/sum(init.n.a.v.i)
                  init.n.a.v.i <- n.a.v.i
                  lam.list <- c(lam.list, lam.a.i)
              }
             plot.lam.list[(min.a+1):max.a, as.numeric(i)] <- lam.list
        }
       }
  	 }
 	   plot.est.plot <- list(NA)
 	   for (i in 1:nplot) {
      if (i != 5) {
        plot.l.i <- as.data.frame(matrix(nrow = length(which(!is.na(plot.lam.list[,i]) == TRUE)), ncol = 2))
        names(plot.l.i) <- c("Age", "lambda")
        plot.l.i$Age <- which(!is.na(plot.lam.list[,i]))
        plot.l.i$lambda <- plot.lam.list[,i][which(!is.na(plot.lam.list[,i]) == TRUE)]
        plot.est.plot[[i]] <- plot.l.i 
      }
      else if (i == 5) {
        plot.l.5 <- as.data.frame(matrix(nrow = 1, ncol = 2))
        names(plot.l.5) <- c("Age", "lambda")
        plot.l.5$Age[1] <- Age.mature
        plot.l.5$lambda[1] <- lam.mature
        plot.est.plot[[5]] <- plot.l.5
      }
 	   }
 }

  # Plot: lambda por parcela (PREDICHA)
 {
 	 lambda.df.1 <- transform(lambda.df, plot = as.factor("total"))
 	 p.lam.p <- as.data.frame(matrix(nrow = 0, ncol = 3))
 	 names(p.lam.p) <- c("Age", "lambda", "plot")
 	 for (i in 1:nplot) {
 		 p.lam.p.i <-transform(plot.est.plot[[i]], plot = as.factor(i))
 		 p.lam.p  <- rbind(p.lam.p, p.lam.p.i)
 	 }
 	 p.lam <- ggplot(p.lam.p, aes(x = Age, y = lambda, col = plot)) + 
 	 geom_line(size = 1) + 
 	 geom_line(data = lambda.df.1, aes(x = Age, y = lambda, col = plot), size = 2, alpha = 1/5) +
 	 labs(x = "Edad sucesional", y = expression(lambda), title = "Dinámica poblacional")
 	 p.lam
 	 #ggsave("lambda-parcela.png", p.lam, device = "png", width = 9, height = 7, units = "in", dpi = 180*2)
 }

 # Lambda por parcela (OBSERVADA)
 {   
     nplot <- nlevels(Survival$PLOT)
     Survival <- transform(Survival, PLOT = as.factor(as.character(plot)))
 	   plot.lam.list <- as.data.frame(matrix(nrow = Age.mature, ncol = length(levels(Survival$PLOT)), NA))
 	   names(plot.lam.list) <- 1:nplot
 	   for (i in levels(Survival$PLOT)) {
  	 	 Mim.s.p <- subset(Survival, PLOT == i)
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

 # Plot: lambda por parcela (OBSERVADA)
 {
 	 lambda.df.2 <- transform(lambda.df, plot = as.factor("total"))
 	 p.lam.p.2 <- as.data.frame(matrix(nrow = 0, ncol = 3))
 	 names(p.lam.p.2) <- c("Age", "lambda", "plot")
 	 for (i in 1:nplot) {
 		 p.lam.p.i <-transform(plot.lam.plot[[i]], plot = as.factor(i))
 		 p.lam.p.2  <- rbind(p.lam.p.2, p.lam.p.i)
 	 }
 	 p.lam.ob <- ggplot(p.lam.p.2, aes(x = Age, y = lambda, col = plot, alpha = 1/3)) + 
 	 geom_line(size = 1) + 
 	 geom_line(data = lambda.df.2, aes(x = Age, y = lambda, col = plot), size = 2, alpha = 1/5) +
 	 labs(x = "Edad sucesional", y = expression(lambda), title = "Dinámica poblacional")
 	 p.lam.ob
 	 #ggsave("lambda-parcela-line.png", p.lam.ob, device = "png", width = 18, height = 12, units = "in", dpi = 180*2)
 }

 # Plot: SIN MIGRACIÓN PARA ARTÍCULO
 {
 	 lambda.df.2 <- transform(lambda.df)
 	 p.lam.p.2 <- as.data.frame(matrix(nrow = 0, ncol = 2))
 	 names(p.lam.p.2) <- c("Age", "lambda")
 	 for (i in 1:nplot) {
 		 p.lam.p.i <-transform(plot.lam.plot[[i]])
 		 p.lam.p.2  <- rbind(p.lam.p.2, p.lam.p.i)
 	 }
 	 p.lam.p.2 <- p.lam.p.2[-c(53),]


 	 p.lam.ob <- ggplot(p.lam.p.2, aes(x = Age, y = lambda, alpha = 1/2), show.legend = FALSE) + 
 	 theme_minimal() +
 	 geom_line(size = 1) + 
     theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 15)) + 
 	 geom_line(data = lambda.df.2, aes(x = Age, y = lambda), col = "red", size = 1, alpha = 1/3, show.legend = FALSE) +
 	 labs(x = expression(paste("Abandonment time ", italic(t), " (years)")), y = expression(lambda))+
 	 scale_alpha(guide = "none")
 	 p.lam.ob
 	 ggsave("no-migration.png", p.lam.ob, device = "png", width = 18, height = 12, units = "in", dpi = 180*2)
 }


 # Lambda por parcela (PREDICHA por PARCELA) !!!
 {   
     load("k-p-list.RData")
     Survival <- transform(Survival, PLOT = as.factor(as.character(plot)))
     nplot <- nlevels(Survival$PLOT)
     est.lam.df <- as.data.frame(matrix(nrow = Age.mature, ncol = nlevels(Survival$PLOT), NA))
     names(est.lam.df) <- 1:nplot
     for (i in levels(Survival$PLOT)) {
       Mim.s.p <- subset(Survival, PLOT == i & !is.na(ln.h2))
         if (Mim.s.p$PLOT[1] != 5) {
          if (all(!is.na(Mim.s.p$Age))) {
             lam.list <- c()
             min.a <- min(Mim.s.p$Age, na.rm = TRUE)
             max.a <- max(Mim.s.p$Age, na.rm = TRUE)
             n.1.i <- which(Mim.s.p$Age == min.a)
             n.1.h.i <- log(Mim.s.p[n.1.i,]$h2)#[order(log(Mim.s.p[n.1.i,]$h2))]
             init.n.a.v.i <- n.1.v.i <- hist(n.1.h.i, breaks = e.pred, plot = FALSE)$counts
             for (a in (min.a+1):max.a) {
                  n.a.v.i <- k.p.list[[as.integer(i)]][(a-min.a),,]%*%init.n.a.v.i
                  lam.a.i <- sum(n.a.v.i)/sum(init.n.a.v.i)
                  init.n.a.v.i <- n.a.v.i
                  lam.list <- c(lam.list, lam.a.i)
              }
             est.lam.df[(min.a+1):max.a, as.numeric(i)] <- lam.list
        }
       }
     }
     est.lam.list <- list(NA)
     for (i in 1:nplot) {
      if (i != 5) {
        plot.l.i <- as.data.frame(matrix(nrow = length(which(!is.na(est.lam.df[,i]) == TRUE)), ncol = 2))
        names(plot.l.i) <- c("Age", "lambda")
        plot.l.i$Age <- which(!is.na(est.lam.df[,i]))
        plot.l.i$lambda <- est.lam.df[,i][which(!is.na(est.lam.df[,i]) == TRUE)]
        est.lam.list[[i]] <- plot.l.i 
      } else {
        plot.l.5 <- as.data.frame(matrix(nrow = 1, ncol = 2))
        names(plot.l.5) <- c("Age", "lambda")
        plot.l.5$Age[1] <- Age.mature
        plot.l.5$lambda[1] <- lam.mature
        est.lam.list[[5]] <- plot.l.5
      }
     }
 }

 # Lambda por parcela (OBSERVADA por PARCELA) !!!
 {   
     nplot <- nlevels(Survival$PLOT)
     Survival <- transform(Survival, PLOT = as.factor(as.character(plot)))
     obs.lam.df <- as.data.frame(matrix(nrow = Age.mature, ncol = length(levels(Survival$PLOT)), NA))
     names(obs.lam.df) <- 1:nplot
     for (i in levels(Survival$PLOT)) {
       Mim.s.p <- subset(Survival, PLOT == i & !is.na(h2))
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
             obs.lam.df[(min.a+1):max.a, as.numeric(i)] <- lam.list
        }
       }
     }
     obs.lam.list <- list(NA)
     for (i in 1:nplot) {
      if (i != 5) {
        plot.l.i <- as.data.frame(matrix(nrow = length(which(!is.na(obs.lam.df[,i]) == TRUE)), ncol = 2))
        names(plot.l.i) <- c("Age", "lambda")
        plot.l.i$Age <- which(!is.na(obs.lam.df[,i]))
        plot.l.i$lambda <- obs.lam.df[,i][which(!is.na(obs.lam.df[,i]) == TRUE)]
        obs.lam.list[[i]] <- plot.l.i 
      }
      else if (i == 5) {
        plot.l.5 <- as.data.frame(matrix(nrow = 1, ncol = 2))
        names(plot.l.5) <- c("Age", "lambda")
        plot.l.5$Age[1] <- Age.mature
        plot.l.5$lambda[1] <- lam.mature
        obs.lam.list[[5]] <- plot.l.5
      }
     }
 }

 # Plot: lambda por parcela (GUSANITOS) !!!
 {
   		lambda.df.2 <- transform(lambda.df, plot = as.factor("total"))
   		est.lam.p <- as.data.frame(matrix(nrow = 0, ncol = 3))
   		names(est.lam.p) <- c("Age", "lambda", "plot")
   		for (i in 1:nplot) {
     		est.lam.p.i <-transform(est.lam.list[[i]], plot = as.factor(i))
     		est.lam.p  <- rbind(est.lam.p, est.lam.p.i)
   		}
   		obs.lam.p <- as.data.frame(matrix(nrow = 0, ncol = 3))
   		names(est.lam.p) <- c("Age", "lambda", "plot")
   		for (i in 1:nplot) {
     		obs.lam.p.i <-transform(obs.lam.list[[i]], plot = as.factor(i))
     		obs.lam.p  <- rbind(obs.lam.p, obs.lam.p.i)
   		}
   		plot.l.graf <- ggplot(lambda.df.2, aes(x = Age, y = lambda, col = plot)) +
   		theme_minimal() +
   		xlab(expression(paste("Edad sucesionl ", italic(t), " (años)"))) +
   		ylab(expression(lambda)) +
	 	theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 15)) +
   		geom_line(size = 2, alpha = 0.3) + 
   		geom_line(data = est.lam.p, aes(x = Age, y = lambda, col = plot, alpha = 1/3), size = 1) +
   		geom_line(data = obs.lam.p, linetype = "dotted", aes(x = Age, y = lambda, col = plot), size = 1) +
   		scale_alpha(guide = "none") +
   		guides(col = guide_legend(title = "Parcela"))
   		plot.l.graf
   		ggsave("lambda-parcela-line.png", plot.l.graf, device = "png", width = 10, height = 7, units = "in", dpi = 180*2)
 }
}

# ESTIMACIÓN DE INMIGRACIÓN
{ 
  load("k-p-list.RData") #cargar k.p.list
  lambda.obs <- plot.lam.plot
  lambda.est <- plot.est.plot
  distance.l <- function(c) {
    dist.l <- 0
    lam.est <- c()
    for (p in 1:nplot)
        if (p != 5) {
            S.p <- subset(Survival, plot == p & !is.na(ln.h2))
            S.p <- droplevels(S.p)
            min.a <- min(S.p$Age)
            max.a <- max(S.p$Age)
            n.1 <- which(S.p$Age == min.a)
            n.1.h <- log(S.p[n.1,]$h2)
            n.1.v <- hist(n.1.h, breaks = e.pred, plot = FALSE)$counts
            init.n.a.v <- n.1.v+c*F5
            for (a in (min.a+1):max.a) {
                n.a.v <- k.p.list[[p]][(a-min.a),,]%*%init.n.a.v
                lam.a <- (sum(n.a.v)+c)/sum(init.n.a.v)
                init.n.a.v <- n.a.v+c*F5
                lam.est <- c(lam.est, lam.a)
            }
            lambda.obs.plot <- rep(NA, nrow(lambda.obs[[p]]))
            for (j in (1:nrow(lambda.obs[[p]])))
                dist.l <- dist.l + (lam.est[j]-lambda.obs[[p]]$lambda[j])^2
        dist.l <- sqrt(dist.l)
        cat(paste0("c = ", c, ", dist = ", dist.l, "\n"))
       }
    return(dist.l)
  } 

  # Probar diferentes estimadores
  {
    opt.brent <- optim(0, distance.l, method = "Brent", lower = 0, upper = 100)
    opt.brent$par # 10.79303, dist = 1.702837
    opt.lbf <- optim(0, distance.l, method = "L-BFGS-B")
    c <- opt.lbf$par # 10.79292, dist = 1.702837
    opt.bfgs <- optim(0, distance.l, method = "BFGS")
    opt.bfgs$par # 10.64919, dist = 1.702837
  }

  # Plot: lambda por parcela (GUSANITOS) !!!
  {
   # Lambda total predicha
   {
     lam.list <- c()
     n.0 <- which(Survival$Age == 0)
     n.0.h <- log(Survival[n.0,]$h2)
     n.0.v <- hist(n.0.h, breaks = e.pred, plot = FALSE)$counts
     init.n.a.v <- n.0.v
     n.list <- c(sum(init.n.a.v))
     for (a in 1:Age.mature) {
       n.a.v <- k.i.j.a[a,,]%*%init.n.a.v
       lam.a <- (sum(n.a.v)+c)/sum(init.n.a.v)
       init.n.a.v <- n.a.v+c*F5
       n.list <- c(n.list, sum(n.a.v))
       lam.list <- c(lam.list, lam.a)
       }
      tot.lam.pred <- lam.list
   } 

   lambda.df <- as.data.frame(list(lambda = lam.list, Age = 1:Age.mature))
   lambda.df.2 <- transform(lambda.df, plot = as.factor("total"))
   est.lam.p <- as.data.frame(matrix(nrow = 0, ncol = 3))
   names(est.lam.p) <- c("Age", "lambda", "plot")
   for (i in 1:nplot) {
     est.lam.p.i <-transform(est.lam.list[[i]], plot = as.factor(i))
     est.lam.p  <- rbind(est.lam.p, est.lam.p.i)
   }
   obs.lam.p <- as.data.frame(matrix(nrow = 0, ncol = 3))
   names(est.lam.p) <- c("Age", "lambda", "plot")
   for (i in 1:nplot) {
     obs.lam.p.i <-transform(obs.lam.list[[i]], plot = as.factor(i))
     obs.lam.p  <- rbind(obs.lam.p, obs.lam.p.i)
   }

   plot.l.graf.c <- ggplot(lambda.df.2, aes(x = Age, y = lambda, col = plot)) +
   theme_minimal() +
   xlab(expression(paste("Edad sucesionl ", italic(t), " (años)"))) +
   ylab(expression(lambda)) +
   theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 15)) + 
   geom_line(size = 2, alpha = 0.3) + 
   geom_line(data = est.lam.p, aes(x = Age, y = lambda, col = plot, alpha = 1/3), size = 1) +
   geom_line(data = obs.lam.p, linetype = "dotted", aes(x = Age, y = lambda, col = plot), size = 1) +
   scale_alpha(guide = "none") +
   guides(col = guide_legend(title = "Parcela"))
   plot.l.graf.c
   ggsave("gusanitos.png", plot.l.graf.c, device = "png", width = 10, height = 7, units = "in", dpi = 180*2)
  }
}

 # Plot: lambda por parcela (GUSANITOS) !!!
  {
   # Lambda total predicha
   {
     lam.list <- c()
     n.0 <- which(Survival$Age == 0)
     n.0.h <- log(Survival[n.0,]$h2)
     n.0.v <- hist(n.0.h, breaks = e.pred, plot = FALSE)$counts
     init.n.a.v <- n.0.v
     n.list <- c(sum(init.n.a.v))
     for (a in 1:Age.mature) {
       n.a.v <- k.i.j.a[a,,]%*%init.n.a.v
       lam.a <- (sum(n.a.v)+c)/sum(init.n.a.v)
       init.n.a.v <- n.a.v+c*F5
       n.list <- c(n.list, sum(n.a.v))
       lam.list <- c(lam.list, lam.a)
       }
      tot.lam.pred <- lam.list
   } 

   ## CON MIGRACIÓN PARA ARTÍCULO

   plot.l.graf.c <- ggplot(lambda.df.2, aes(x = Age, y = lambda)) +
   theme_minimal() +
   xlab(expression(paste("Abandonment time ", italic(t), " (years)"))) +
   ylab(expression(lambda)) +
   theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 15)) + 
   geom_line(size = 2, alpha = 1/3, col = "red") + 
   geom_line(data = obs.lam.p, aes(x = Age, y = lambda), alpha = 1/2, size = 1) +
   scale_alpha(guide = "none")
   plot.l.graf.c
   ggsave("w-migration.png", plot.l.graf.c, device = "png", width = 10, height = 7, units = "in", dpi = 180*2)
  }
}


# PLOTS
{	
 # Definiciones necesarias
 {
	 bin.size <- 0.1
	 n.sim <- 1000
	 h1.pred <- seq(min(newSurvival$h1, na.rm = TRUE), ceiling(max(newSurvival$h1, na.rm = TRUE)), length.out = m+1)
	 h2.pred <- seq(min(newSurvival$h2, na.rm = TRUE), ceiling(max(newSurvival$h2, na.rm = TRUE)), length.out = m+1)
	 mids.h2 <- hist(newSurvival$h2, breaks = h2.pred, plot = FALSE)$mids
	 ln.mids.h2 <- log(mids.h2)
	 lh1.pred <- log(h1.pred)
	 data.s.j <- subset(newSurvival, !is.na(h1) & !is.na(sup))
	 data.s.j <- droplevels(data.s.j)
	 data.g.j <- subset(data.s.j, !is.na(h2))
	 data.g.j <- droplevels(data.g.j)
	 mids <- hist(newSurvival$h1, breaks = h1.pred, plot = FALSE)$mids
	 ln.mids <- log(mids)
	 age.pred.f1 <- 1:Age.mature
	 age.pred.j <- 1:Age.mature
 }
  
 # Survival
 {
	 data.s.j.1 <- subset(data.s.j, sup == 1)
	 count.s.j.1 <- hist(data.s.j.1$h1, breaks = h1.pred, plot = FALSE)$counts
	 count.s.j <- hist(data.s.j$h1, breaks = h1.pred, plot = FALSE)$counts
	 prop.s.j <- count.s.j.1/count.s.j
	 plot.data.s.j <- as.data.frame(matrix(ncol = 4, nrow = 4))
	 names(plot.data.s.j) <- c("age", "h1", "prop", "n")
	 for (a in 1:length(age.pred.j)) {
	 	data.a <- subset(data.s.j, Age == age.pred.j[a])
	 	data.a.1 <- subset(data.a, sup == 1)
	 	count.a.1 <- hist(data.a.1$h1, breaks = h1.pred, plot = FALSE)$counts
	 	count.a <- hist(data.a$h1, breaks = h1.pred, plot = FALSE)$counts
	 	plot.data.s.j[(a-1)*length(mids)+1:length(mids), ] <- cbind(rep(age.pred.j[a], length(mids)), mids, count.a.1/count.a, count.a)
	 }
	 plot.data.s.j <- subset(plot.data.s.j, !is.na(prop))
	 pred.s.j <- as.data.frame(matrix(ncol = 3, nrow = length(age.pred.j)*length(mids)))
	 names(pred.s.j) <- c("age.pred", "h1.pred", "s.pred")
	 for (a in 1:length(age.pred.j)) {
	 	s.pred <- as.vector(predict(glmm.s, newdata = data.frame(ln.h1 = ln.mids, Age = rep(age.pred.j[a], length(mids))), type = "response", re.form = NA))
	 	pred.s.j[(a-1)*length(mids)+1:length(mids), ] <- cbind(rep(age.pred.j[a], length(mids)), mids, s.pred)
	 }
	 g.s.j <- ggplot(pred.s.j) +
	 theme_minimal() +
	 theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
	 xlab(expression(paste("Edad sucesionl ", italic(t), " (años)"))) +
	 ylab(expression("Altura (m)")) +
	 scale_x_continuous(expand = c(0,0)) +
	 scale_y_continuous(expand = c(0,0)) +
	 scale_fill_gradient(low = "#1569f1", high = "#f21515", guide = "none") +
	 scale_alpha(guide = "none") +
	 geom_raster(aes(age.pred, h1.pred, fill = s.pred, alpha = 0.35)) +
	 scale_colour_gradient(low = "#1569f1", high = "#f21515") +
	 geom_point(data = plot.data.s.j, aes(age, h1, color = prop, size = n), shape = 21) +
	 guides(colour = guide_colourbar(title = "Probabilidad de\nsupervivencia", order = 1), size = guide_legend(title = "Tamaño de\nmuestra", order = 2))
	 dev.new(width = 10, height = 8)
	 g.s.j
	 #ggsave("/home/gerardo/MEGA/Tesis/Gráficas/sur-gg.png", g.s.j, device = "png", width = 10, height = 7, units = "in", dpi = 180*2)	
 }
 
 # Growth
 {
	 # ggplot 2 version
	 sdg <- sigma(lmm.g)
	 setwd("/Users/gdh/Google\ Drive/L\ 2017\ Gerardo\ Dávila/Gráficas/2 Growth")
	 for (a in 1:Age.mature) {
	  pred.g.j <- as.data.frame(matrix(ncol = 3, nrow = m*m))
	  names(pred.g.j) <- c("x", "y", "prob")
	  lh2.predicted <- as.vector(predict(lmm.g, newdata = data.frame(ln.h1 = ln.mids, Age = rep(a,m)), type = "response", re.form = NA))
	  for (i in 1:m) {
	   pred.g.j[(i-1)*m+1, ] <- cbind(mids[i], mids[1], pnorm(lh1.pred[2], lh2.predicted[i], sdg))
	   pred.g.j[(i-1)*m+m, ] <- cbind(mids[i], mids[m], 1-pnorm(lh1.pred[m], lh2.predicted[i], sdg))
	   for (j in 2:(m-1))
	    pred.g.j[(i-1)*m+j, ] <- cbind(mids[i], mids[j], pnorm(lh1.pred[j+1],lh2.predicted[i], sdg)-pnorm(lh1.pred[j],lh2.predicted[i], sdg))
	  }
	  data.g.j.t <- subset(data.g.j, Age == a)
	  gg.j.t <- ggplot(pred.g.j) +
	  theme_minimal() +
	  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
	  xlab(expression(paste("Altura en el tiempo ", italic(t), " (m)"))) +
	  ylab(expression(paste("Altura en el tiempo ", italic(t), "+1 (m)"))) +
	  scale_x_continuous(expand = c(0,0)) +
	  scale_y_continuous(expand = c(0,0)) +
	  scale_fill_gradient(low = "#1569f1", high = "#f21515", guide = guide_legend(title = "Probabilidad")) +
	  scale_alpha(guide = "none") +
	  geom_raster(aes(x, y, fill = prob, alpha = 0.35)) +
	  scale_colour_gradient(low = "#1569f1", high = "#f21515") +
	  geom_point(data = data.g.j.t, aes(h1, h2), shape = 21) +
	  guides(colour = guide_colourbar(title = "Crecimiento (m)")) +
	  dev.new(width = 10, heigth = 8)
	  #gg.j.t
	  #ggsave(paste0("/home/gerardo/MEGA/Tesis/Gráficas/grt-gg-", a, ".png"), gg.j.t, device = "png", width = 10, height = 7, units = "in", dpi = 180*2)
	 }

	 #setwd("/Users/gdh/Google\ Drive/L\ 2017\ Gerardo\ Dávila/Gráficas/2 Growth/GIF1")
	 #system("convert -delay 0.25 *.png gg.j.gif")
	 
	 #  Plot 3D version
	 setwd("/Users/gdh/Google\ Drive/L\ 2017\ Gerardo\ Dávila/Gráficas/2 Growth/GIF1")
	 for (a in 1:Age.mature) {
	  png(paste0("G.", a, ".png"))
	  hist3D(x.pred, x.pred, G[a,,], scale = TRUE, col = colfunc(100))
	  dev.off()
	 }
	 #system("convert -delay .4 *.png G.gif")	 
 }
 
 # Fecundity 1
 {
	 data.rep <- subset(Mim.rep, !is.na(h1) & !is.na(sup))
	 data.rep <- droplevels(data.rep)
	 data.rep1 <- subset(Mim.rep, Rep == 1)
	 count.rep1 <- hist(data.rep1$h1, breaks = h1.pred, plot = FALSE)$counts
	 count.rep <- hist(data.rep$h1, breaks = h1.pred, plot = FALSE)$counts
	 prop.rep <- count.rep1/count.rep
	 plot.data.rep <- as.data.frame(matrix(ncol = 4, nrow = 4))
	 names(plot.data.rep) <- c("age", "h1", "prop.f1", "n")
	 for (a in 1:length(age.pred.f1)) { # EJG: poner los años que tienen sentido
	     data.rep.a <- subset(data.rep, Age == age.pred.f1[a])
	     data.rep.a.1 <- subset(data.rep.a, Rep == 1)
	     count.rep.a.1 <- hist(data.rep.a.1$h1, breaks = h1.pred, plot = FALSE)$counts
	     count.rep.a <- hist(data.rep.a$h1, breaks = h1.pred, plot = FALSE)$counts
	     plot.data.rep[(a-1)*length(mids)+1:length(mids), ] <- cbind(rep(age.pred.f1[a], length(mids)), mids, count.rep.a.1/count.rep.a, count.rep.a)
	     }
	 plot.data.rep <- subset(plot.data.rep, !is.na(prop.f1))
	 pred.f1.j <- as.data.frame(matrix(ncol = 3, nrow = length(age.pred.f1)*length(mids)))
	 names(pred.f1.j) <- c("age.pred", "h1.pred", "f1.pred")
	 for (a in 1:length(age.pred.f1)) {
	  pred.f1 <- as.vector(predict(glmm.f1, newdata = data.frame(ln.h1 = ln.mids, Age = rep(age.pred.f1[a], length(mids))), type = "response", re.form = NA))
	  pred.f1.j[(a-1)*length(mids)+1:length(mids), ] <- cbind(rep(age.pred.f1[a], length(mids)), mids, pred.f1)
	  }
	 g.f1.j <- ggplot(pred.f1.j) +
	 theme_minimal() +
	 theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
	 xlab(expression(paste("Edad sucesional ", italic(t), " (años)"))) +
	 ylab(expression("Altura (m)")) +
	 scale_x_continuous(expand = c(0,0)) +
	 scale_y_continuous(expand = c(0,0)) +
	 scale_fill_gradient(low = "#1569f1", high = "#f21515", guide = "none") +
	 scale_alpha(guide = "none") +
	 geom_raster(aes(age.pred, h1.pred, fill = f1.pred, alpha = 0.35)) +
	 scale_colour_gradient(low = "#1569f1", high = "#f21515") +
	 geom_point(data = plot.data.rep, aes(age, h1, color = prop.f1, size = n), alpha = 5/10) +
	 guides(colour = guide_colourbar(title = "Probabilidad de\nreproducción", order = 1), size = guide_legend(title = "Tamaño de\nmuestra", order = 2))
	 dev.new(width = 10, height = 8)
	 g.f1.j
	 #ggsave("/home/gerardo/MEGA/Tesis/Gráficas/fc1-gg.png", g.f1.j, device = "png", width = 10, height = 7, units = "in", dpi = 180*2)	
 }
 
 # Fecundity 2
 {
	 plot.data.f2 <- as.data.frame(matrix(ncol = 3, nrow = m*length(age.pred.j)))
	 names(plot.data.f2) <- c("age.pred.f2", "h1.pred.f2", "frt.n.pred")
	 for (a in 1:length(age.pred.j)) {
	     frt.n.pred <- as.vector(predict(gam.f2$gam, newdata = data.frame(ln.h1 = x.pred, Age = rep(age.pred.j[a], m)), type = "response", re.form = NA))
	     plot.data.f2[(a-1)*m+1:m, ] <- cbind(rep(age.pred.j[a], m), mids, frt.n.pred)
	     }
	 real.rep.data <- subset(Mim.frut, Rep == 1)
	 my_varsf2 <- c("Age", "h1", "TotFrut")
	 real.rep.data <- real.rep.data[my_varsf2]
	 g.f2.j <- ggplot(plot.data.f2) +
	 theme_minimal() +
	 theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
	 xlab(expression(paste("Edad sucesional ", italic(t), " (años)"))) +
	 ylab(expression("Altura (m)")) +
	 scale_x_continuous(expand = c(0,0)) +
	 scale_y_continuous(expand = c(0,0)) +
	 scale_fill_gradient(low = "#1569f1", high = "#f21515") +
	 scale_alpha(guide = "none") +
	 geom_raster(aes(age.pred.f2, h1.pred.f2, fill = frt.n.pred, alpha = 0.35)) +
	 geom_point(data = real.rep.data, aes(Age, h1, color = TotFrut), size = 10, alpha = 4/10) +
	 scale_colour_gradient(low = "#1569f1", high = "#f21515") +
	 guides(color = guide_legend(title = "Número de frutos\n observado"), fill = guide_legend("Número de\nfrutos predicho"))
	 dev.new(width = 10, height = 8)
	 g.f2.j
	 # ggsave("/Users/gdh/Documents/MEGA/Tesis/Gráficas/fc2-gg.png", g.f2.j, device = "png", width = 10, height = 7, units = "in", dpi = 180*2)
 }
 
 # Fecundity 3
 {
	 #Plot f3
	 plot.data.f3 <- as.data.frame(matrix(ncol = 3, nrow = m*length(age.pred.j)))
	 names(plot.data.f3) <- c("age.pred.f3", "h1.pred.f3", "sd.n.pred")
	 for (a in 1:length(age.pred.j)) {
	     sd.n.pred <- as.vector(predict(glmm.f3, newdata = data.frame(ln.h1 = ln.mids, Age = rep(age.pred.j[a], m)), type = "response", re.form = NA))
	     plot.data.f3[(a-1)*length(mids)+1:length(mids), ] <- cbind(rep(age.pred.j[a], length(mids)), mids, sd.n.pred)
	     }
	 real.rep.data <- Seeds
	 my_varsf3 <- c("Age", "h1", "N.seed")
	 real.rep.data <- real.rep.data[my_varsf3]
	 g.f3.j <- ggplot(plot.data.f3) +
	 theme_minimal() +
	 theme(axis.text = element_text(size = 12), axis.title = element_text(size = 20, face = "bold"), legend.text = element_text(size = 12), legend.title = element_text(size = 20)) +
	 xlab(expression(paste("Edad sucesional ", italic(t), " (años)"))) +
	 ylab(expression("Altura (m)")) +
	 scale_x_continuous(expand = c(0,0)) +
	 scale_y_continuous(expand = c(0,0)) +
	 scale_fill_gradient(low = "#1569f1", high = "#f21515", guide = guide_legend(title = "Número de semillas\npredicho")) +
	 scale_alpha(guide = "none") +
	 geom_raster(aes(age.pred.f3, h1.pred.f3, fill = sd.n.pred, alpha = 0.35)) +
	 geom_point(data = real.rep.data, aes(Age, h1, color = N.seed), size = 10, alpha = 4/10) +
	 scale_colour_gradient(low = "#1569f1", high = "#f21515") +
	 guides(color = guide_legend(title = "Número de semillas\nobservado"), size = guide_legend("Número de semillas\nobservado")) + 
	 #scale_fill_continuous(guide = guide_legend(title = "Número de\nsemillas predicho"))
	 dev.new(width = 10, height = 8)
	 g.f3.j
	 #ggsave("/Users/gdh/Documents/MEGA/Tesis/Gráficas/fc3-gg.png", g.f3.j, device = "png", width = 10, height = 7, units = "in", dpi = 180*2)
 }
 
 # F1234
 {
	 # Plot f12 y  f123
	 my_varsf1 <- c("age.pred", "h1.pred", "f1.pred")
	 my_varsf2 <- c("age.pred.f2", "h1.pred.f2", "frt.n.pred")
	 my_varsf3 <- c("sd.n.pred")
	 plot.data.f12 <- cbind(pred.f1.j[my_varsf1], plot.data.f2[my_varsf2], plot.data.f3[my_varsf3])
	 plot.data.f12 <- transform(plot.data.f12, f12.pred = f1.pred*frt.n.pred)
	 plot.data.f123 <- transform(plot.data.f12, f123.pred = f12.pred*sd.n.pred)
	 stb.pred <- as.data.frame(rep(predict(gam.f4, newdata = data.frame(x = 1:100), type = "response"), Age.mature))
	 colnames(stb.pred) <- c("stb.pred")
	 plot.data.f1234 <- cbind(plot.data.f123, stb.pred)
	 plot.data.f1234 <- transform(plot.data.f1234, f1234.pred = f123.pred*stb.pred)

	 # plot f12
	 real.rep.data <- subset(Mim.frut, Rep == 1)
	 my_varsf12 <- c("Age", "h1", "TotFrut")
	 real.rep.data <- real.rep.data[my_varsf12]
	 g.f12.j <- ggplot(plot.data.f123) +
	 theme_minimal() +
	 theme(axis.text = element_text(size = 10), axis.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 10), legend.title = element_text(size = 15)) +
	 ggtitle(expression(italic(f)[12])) +
     xlab(expression(paste("Edad sucesional ", italic(t), " (años)"))) +
	 ylab(expression("Altura en el tiempo (m)")) +
	 scale_x_continuous(expand = c(0,0)) +
	 scale_y_continuous(expand = c(0,0)) +
	 scale_fill_gradientn(colours = c("#f22a0c", "#f7f445"), guide = guide_legend(title = "Número de frutos\npredicho")) +
	 scale_alpha(guide = "none") +
	 geom_raster(aes(age.pred, h1.pred, fill = f12.pred, alpha = 0.35)) +
	 geom_point(data = real.rep.data, aes(Age, h1, color = TotFrut), size = 10, alpha = 4/10) +
	 scale_colour_gradient(low = "#f22a0c", high = "#f7f445") +
	 guides(color = guide_legend(title = "Número de frutos\nobservado"), size = guide_legend("Número de frutos\nobservado")) 
	 #scale_fill_continuous(guide = guide_legend(title = "Número de\nsemillas predicho"))
	 dev.new(width = 10, height = 8)
	 g.f12.j
	 #ggsave("/home/gerardo/MEGA/Tesis/Gráficas/f12Plot.pdf", g.f12.j, device = "pdf", width = 14, height = 10, units = "in")


	 # plot f123
	 real.rep.data <- Seeds
	 my_varsf3 <- c("Age", "h1", "N.seed")
	 real.rep.data <- real.rep.data[my_varsf3]
	 g.f123.j <- ggplot(plot.data.f123) +
	 theme_linedraw() +
	 theme(axis.text = element_text(size = 15), axis.title = element_text(size = 25, face = "bold"), legend.text = element_text(size = 15), legend.title = element_text(size = 25)) +
	 ggtitle("F123") +
     xlab("Edad de abandono (años)") +
	 ylab("Altura en el tiempo t (m)") +
	 scale_x_continuous(expand = c(0,0)) +
	 scale_y_continuous(expand = c(0,0)) +
	 scale_fill_gradient(low = "#f22a0c", high = "#f7f445", guide = guide_legend(title = "Número de semillas\npredicho")) +
	 scale_alpha(guide = "none") +
	 geom_raster(aes(age.pred, h1.pred, fill = f123.pred, alpha = 0.35)) +
	 geom_point(data = real.rep.data, aes(Age, h1, color = N.seed), size = 10, alpha = 4/10) +
	 scale_colour_gradient(low = "#f22a0c", high = "#f7f445") +
	 guides(color = guide_legend(title = "Número de semillas\nobservado"), size = guide_legend("Número de semillas\nobservado")) + 
	 #scale_fill_continuous(guide = guide_legend(title = "Número de\nsemillas predicho"))
	 dev.new(width = 10, height = 8)
	 g.f123.j
	 #ggsave("/home/gerardo/MEGA/Tesis/Gráficas/f123Plot.pdf", g.f123.j, device = "pdf", width = 14, height = 10, units = "in")


	 #plot f1234
	 #real.rep.data <- Seeds
	 #my_varsf3 <- c("Age", "h1", "N.seed")
	 #real.rep.data <- real.rep.data[my_varsf3]
	 g.f1234.j <- ggplot(plot.data.f1234) +
	 theme_linedraw() +
	 theme(axis.text = element_text(size = 15), axis.title = element_text(size = 25, face = "bold"), legend.text = element_text(size = 15), legend.title = element_text(size = 25)) +
	 ggtitle("F1234") +
     xlab("Edad de abandono (años)") +
	 ylab("Altura en el tiempo t (m)") +
	 scale_x_continuous(expand = c(0,0)) +
	 scale_y_continuous(expand = c(0,0)) +
	 scale_fill_gradient(low = "#f22a0c", high = "#f7f445", guide = guide_legend(title = "Número de reclutas\npredicho")) +
	 scale_alpha(guide = "none") +
	 geom_raster(aes(age.pred, h1.pred, fill = f1234.pred)) +
	 #geom_point(data = real.rep.data, aes(Age, h1, color = N.seed), size = 10, alpha = 4/10) +
	 scale_colour_gradient(low = "#f22a0c", high = "#f7f445") +
	 #guides(color = guide_legend(title = "Número de semillas\nobservado"), size = guide_legend("Número de semillas\nobservado")) + 
	 #scale_fill_continuous(guide = guide_legend(title = "Número de\nsemillas predicho"))
	 dev.new(width = 10, height = 8)
	 g.f1234.j
	 #ggsave("/home/gerardo/MEGA/Tesis/Gráficas/f123Plot.pdf", g.f123.j, device = "pdf", width = 14, height = 10, units = "in")
 }

 # All Fecundity plots
 {
   	 g.f1.j <- g.f1.j + theme(legend.position="none") +  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13, face = "bold"), legend.text = element_text(size = 10), legend.title = element_text(size = 2)) 
   	 g.f2.j <- g.f2.j + theme(legend.position="none") +  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13, face = "bold"), legend.text = element_text(size = 10), legend.title = element_text(size = 2)) 
   	 g.f3.j <- g.f3.j + theme(legend.position="none")+  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13, face = "bold"), legend.text = element_text(size = 10), legend.title = element_text(size = 2))  
   	 g.f12.j <- g.f12.j + theme(legend.position="none") +  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13, face = "bold"), legend.text = element_text(size = 10), legend.title = element_text(size = 2)) 
   	 g.f123.j <- g.f123.j + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 13, face = "bold"), legend.text = element_text(size = 10), legend.title = element_text(size = 2))  
   	 all.f.plots <-  ggarrange(g.f1.j, g.f2.j, g.f3.j, g.f12.j, g.f123.j, 
           labels = c("A", "B", "C", "D", "E"),
           ncol = 2, nrow = 3)
   	 #ggsave("/Users/gdh/Google\ Drive/L\ 2017\ Gerardo\ Dávila/Gráficas/3\ Fecundity/allFplots.pdf", all.f.plots, device = "pdf")
 }

 # Kernel
 {
     # getting maximum value for plot limits
     max.v.list <- list(NA)
     for (a in 1:(Age.mature-1)) {
      long_k.a <- melt(k.i.j.a[a,,])
      max.v <- max(long_k.a$value)
      max.v.list[[a]] <- max.v
     }
     max.value <- max(unlist(max.v.list))
     # ggplot2 version 
     for (a in 1:(Age.mature-1)) {
      long_k.a <- melt(k.i.j.a[a,,])
      k.t.plot <- ggplot(long_k.a, aes(x = Var1, y = Var2)) + 
        geom_raster(aes(fill=value)) + 
        scale_fill_gradient(low="#FFFFFF", high="#000000", breaks = c(0, max.value), labels = c(0, round(max.value, 4)), limits = c(0, max.value)) +
        labs(x = "Altura 1", y = "Altura 2", title = paste0("Año: ", a)) +
        theme_bw() + theme(text = element_text(size = 20*2.5),
                           axis.text.x=element_text(size=20*2.5, angle=0, vjust=0.3),
                           axis.text.y=element_text(size=20*2.5),
                           plot.title=element_text(size=22*2.5)) + theme(legend.position="none")
      #ggsave(paste0("/home/gerardo/MEGA/Tesis/Gráficas/Kernel/",a,"k-1.png"), k.t.plot, device = "png", width = 18, height = 18)
      #ggsave(paste0("/home/gerardo/MEGA/Tesis/Gráficas/Kernel1/",a,"k-1.png"), k.t.plot, device = "png")
      #ggsave(paste0("/Users/gdh/Google\ Drive/L\ 2017\ Gerardo\ Dávila/Gráficas/4 Kernel/",a,"k-1.png"), k.t.plot, device = "png")
     }
     #setwd("/Users/gdh/Google\ Drive/L\ 2017\ Gerardo\ Dávila/Gráficas/4 Kernel/GIF")
     #system("convert -delay 0.5 *.png G.gif")
  
     #Plot 3D version
     #setwd("/Users/gdh/Google\ Drive/L\ 2017\ Gerardo\ Dávila/Gráficas/4 Kernel/GIF2")
     setwd("/home/gerardo/MEGA/Tesis/Gráficas/Kernel2")
     #setwd("/Users/Edgar/Google\ Drive/Edgar/Trabajo/FC/Dirección\ de\ tesis/L\ 2017\ Gerardo\ Dávila/Gráficas/4 Kernel/GIF2")
     for (a in 1:(Age.mature-1)) {
      png(paste0("K.", a, ".png"))
      hist3D(x.pred, x.pred, k.i.j.a[a,,], zlim = c(0, max.value))
      dev.off()
     }
     #setwd("/Users/gdh/Google\ Drive/L\ 2017\ Gerardo\ Dávila/Gráficas/4 Kernel/GIF2")
     #system("convert -delay 0.5 *.png G.gif")
 }
 
 # Lambdas
 {
     lambda.df <- as.dR
     ata.frame(list(lambda = lam.list, Age = seq(0,Age.mature,1)))
     lambda.plot <- qplot(x = Age, y = lambda, data = lambda.df, color = "red", geom = "line", xlab = "Age", ylab = "Lambda", main = "Population dynamics") + scale_fill_discrete(guide=FALSE) + geom_point(data = ob.lam.df, mapping = aes(x = Age, y = ob.lambda))
     lambda.plot
     ggsave("/Users/gdh/Google\ Drive/L\ 2017\ Gerardo\ Dávila/Gráficas/Lambda.plot.png", lambda.plot, device = "png",width = 18, height = 12, units = "in")
     #ggsave("/home/gerardo/MEGA/Tesis/Gráficas/Lambda/lambda.plot.1.png", lambda.plot, device = "png")
 }

 # Plot lambdas
 {
   for (i in 1:nrow(plot.lam.list)) {
    lambda.df <- as.data.frame(list(lambda = plot.lam.list[i,], Age = seq(min(Mim.s$PLOT(,82,1)))
   }
 }
}

# Graficas f y p con  plot3D
{
  png("f-plot3d.png", width = 6, height = 6, units = 'in', res = 300)
  f.3d <- persp3D(z = f.i.a, theta = 220, xlab = "Altura", ylab ="Edad sucesional", zlab = "Reclutas")
  dev.off()
  for(i in c(1,25,40,60,75,99)) {
    png(paste0("p-plot3d-", i, ".png"), width = 4, height = 4 , units = 'in', res = 300)
    persp3D(z = p.i.j.a[i,,], theta = 220, xlab = "Altura t", ylab = "Altura t+1", zlab = "Probabilidad", title = i)
    dev.off()
  }
  for(i in c(1,25,40,60,75,99)) {
    png(paste0("k-plot3d-", i, ".png"), width = 6, height = 6, units = 'in', res = 300)
    persp3D(z = k.i.j.a[i,,], theta = 220, xlab = "Altura t", ylab = "Altura t+1", zlab = "Probabilidad")
    dev.off()
  }
} 

# Graficas f y p con  ggplot
{
     # f
     long_f.a <- melt(f.i.a)
     for(i in 1:m)
     long_f.a$Var1 <- rep(x.pred, m)
      f.t.plot <- ggplot(long_f.a, aes(x = Var1, y = Var2)) + 
        geom_raster(aes(fill=value)) + 
        scale_fill_gradient(low="#FFFFFF", high="#000000", breaks = c(0, max.value), labels = c(0, round(max.value, 4)), limits = c(0, max.value)) +
        labs(x = expression("Altura (m)"), y = expression(paste("Edad sucesional ", italic(t)))) +
        theme_bw() + theme(text = element_text(size = 20),
                           axis.text.x=element_text(size=20, angle=0, vjust=0.3),
                           axis.text.y=element_text(size=20),
                           plot.title=element_text(size=22)) +
        scale_fill_continuous(guide = "colourbar") +
        guides(fill = guide_legend(title = "Prob."))
      f.t.plot
      ggsave("/Users/gdh/Documents/MEGA/Tesis/Gráficas/f-gg.png", f.t.plot, device = "png", width = 5, height = 5)
        

     # k
     max.v.list <- list(NA)
     for (a in 1:(Age.mature-1)) {
      long_k.a <- melt(k.i.j.a[a,,])
      max.v <- max(long_k.a$value)
      max.v.list[[a]] <- max.v
     }
     max.value <- max(unlist(max.v.list))
     # ggplot2 version 
     for (a in c(1,25,40,60,75,99)) {
      long_k.a <- melt(k.i.j.a[a,,])
      k.t.plot <- ggplot(long_k.a, aes(x = Var1, y = Var2)) + 
        geom_raster(aes(fill=value)) + 
        scale_fill_gradient(low="#FFFFFF", high="#000000", breaks = c(0, max.value), labels = c(0, round(max.value, 4)), limits = c(0, max.value)) +
        labs(x = expression(paste("Altura ", italic(t))), y = expression(paste("Altura ", italic(t), "+1")), title = paste0("año = ", a)) +
        theme_bw() + theme(text = element_text(size = 20),
                           axis.text.x=element_text(size=20, angle=0, vjust=0.3),
                           axis.text.y=element_text(size=20),
                           plot.title=element_text(size=22)) + theme(legend.position="none") #+
        #scale_fill_continuous(guide = "colourbar")
      k.t.plot
      ggsave(paste0("/Users/gdh/Documents/MEGA/Tesis/Gráficas/",a,"k.png"), k.t.plot, device = "png", width = 5, height = 5)
     }

     # p
     plist = list()
     max.v.list <- list(NA)
     for (a in 1:(Age.mature-1)) {
      long_p.a <- melt(p.i.j.a[a,,])
      max.v <- max(long_p.a$value)
      max.v.list[[a]] <- max.v
     }
     max.value <- max(unlist(max.v.list))
     # ggplot2 version 
     for (a in c(1,25,40,60,75,99)) {
      long_p.a <- melt(p.i.j.a[a,,])
      for(i in 1:m) {
       long_k.a$Var1 <- rep(x.pred, m)
       long_k.a$Var2 <- rep(x.pred, m)
      }
      p.t.plot <- ggplot(long_p.a, aes(x = Var1, y = Var2)) + 
        geom_raster(aes(fill=value)) + 
        scale_fill_gradient(low="#FFFFFF", high="#000000", breaks = c(0, max.value), labels = c(0, round(max.value, 4)), limits = c(0, max.value)) +
        labs(x = expression(paste("Altura ", italic(t))), y = expression(paste("Altura ", italic(t), "+1")), title = paste0("año = ", a)) +
        theme_bw() + theme(text = element_text(size = 20),
                           axis.text.x=element_text(size=20, angle=0, vjust=0.3),
                           axis.text.y=element_text(size=20),
                           plot.title=element_text(size=22)) + theme(legend.position="none")# +
        #scale_fill_continuous(guide = "colourbar")
      p.t.plot
      ggsave(paste0("/Users/gdh/Documents/MEGA/Tesis/Gráficas/",a,"p.png"), p.t.plot, device = "png", width = 5, height = 5)
     } 
}

# GRÁFICAS CON Plotly
{
  #G promedio
  {
  	G.m.m <- matrix(NA, ncol = 100, nrow = 100)
    for (a in 1:Age.mature) {
      g.m <- predict(lmm.g, newdata = data.frame(ln.h1 = x.pred, Age = rep(a, m)), type = "response", re.form = NA)
      G.m.m[a, ] <- exp(g.m)-exp(x.pred)
    }

    x <- 1:100
    y <- exp(x.pred)
    z <- G.m.m
    plot_ly(x=x,y=y,z=z, type="surface") #%>% layout(
      #title = "Crecimiento promedio",
      #scene = list(
      #xaxis = list(title="Edad sucesional (años)"),
      #yaxis = list(title="Altura (m)"),
      #zaxis = list(title="x-y")
    ))
  }

  #S
  {
     x <- 1:100
     y <- exp(x.pred)
     z <- S
     plot_ly(x=x,y=y,z=z, type="surface") %>% layout(
       title = "Supervivencia",
       scene = list(
       xaxis = list(title="Edad sucesional (años)"),
       yaxis = list(title="Altura (m)"),
       zaxis = list(title="Prob. de supervivencia")
       ))
  }

  #f(x,t)
  {
  	 x <- 1:100
     y <- exp(x.pred)
     z <- f.i.a
     plot_ly(x=x,y=y,z=z, type="surface") %>% layout(
       title = "Fecundidad",
       scene = list(
       xaxis = list(title="Edad sucesional (años)"),
       yaxis = list(title="Altura (m)"),
       zaxis = list(title="Reclutas")
       ))
  }

  #k(x,y,t) solo años 1 25 50 75 y 99
  {
  	 for (i in c(1,25,50,75,99)) {
  	 	 x <- x.pred
     	 y <- x.pred
     	 z <- k.i.j.a[i,,]
     	 plot_ly(x=x,y=y,z=z, type="surface") %>% layout(
       		 title = "k(x,y,t)",
       		 scene = list(
       		 xaxis = list(title="Altura t (m)"),
       	 	 yaxis = list(title="Altura t+1 (m)"),
       		 zaxis = list(title="Probabilidad")
       		))
  	 }
  }
}

# NO CORRER
# FUNCIÓN PARA CALCULAR h ANTES DEL PRIMER REGISTRO
{
 ## ln.h1 = altura antes de ser registrado por primera vez, ln.h2 = altura en su primer registro
 geth1 <- function(ln.h2, Age, Plot, Census) { ## Census? plot? id?
  b0 <- fx.coef[1] + pt.coef$int[which(pt.coef$plot == Plot)] + ifelse(Census != 0, cs.coef$int[which(cs.coef$Census == Census)], 0)
  b1 <- fx.coef[2] + pt.coef$h1[which(pt.coef$plot == Plot)] + ifelse(Census != 0, cs.coef$h1[which(cs.coef$Census == Census)], 0)
  b2 <- fx.coef[3]
  b3 <- fx.coef[4]
	 ln.h1 <- (ln.h2-b0-b2*Age)/(b1+b3*Age) #calcular el modelo general
	 return(ln.h1)
	 }
}

# CALCULAR TAMAÑO ANTERIOR PARA PRIMER REGISTRO
{
	fx.coef <- fixef(lmm.g)
	pt.coef <- ranef(lmm.g)$plot
	pt.coef <- transform(pt.coef, plot = as.factor(row.names(pt.coef)))
	colnames(pt.coef) <- c("int", "h1", "plot")
	cs.coef <- ranef(lmm.g)$Census
	cs.coef <- transform(cs.coef, Census = as.factor(row.names(cs.coef)))
	colnames(cs.coef) <- c("int", "h1", "Census")
	censuses <- as.factor(row.names(cs.coef))
	newSurvival <- Survival
	for (i in levels(Survival$id)) {
		S.i <- subset(Survival, id == i)
		ln.h2 <- S.i$ln.h2[1]
		Age <- S.i$Age[1]
 		Plot <- S.i$plot[1]
 		Census <- ifelse(any(censuses == S.i$Census[1]-1) == TRUE, S.i$Census[1]-1, 0)
 		if(Census != 0) {
 		newSurvival[which(newSurvival$id == i)[1],]$ln.h1[1] <- geth1(ln.h2, Age, Plot, Census)
		newSurvival[which(newSurvival$id == i)[1],]$h1[1] <- exp(newSurvival[which(newSurvival$id == i)[1],]$ln.h1[1])
 		}
 		if(Census == 0) {
  		ln.h1.new <- geth1(ln.h2, Age, Plot, Census)
  		new.row <- as.data.frame(c(S.i$Age[1]-1, ln.h1.new, ln.h2, exp(ln.h1.new), exp(ln.h2), Plot, 1, 0, S.i$Census[1]-1))
  		rownames(new.row) <- colnames(newSurvival)
  		colnames(new.row) <- "1"
  		new.row[which(rownames(new.row) == "id"), 1] <- i
  		row <- as.data.frame(matrix(rep(NA, length(rownames(new.row))), nrow = 1, ncol = length(rownames(new.row))))
  		for (j in 1:length(row)) {
   			if (rownames(new.row)[j] != "id") {
    		row[1, j] <- as.numeric(new.row[j, 1])
   			} else {
    		row[1, j] <- new.row[j, 1]
   			}
  		}
  		colnames(row) <- rownames(new.row)
  		newSurvival <- rbind(newSurvival, row)
 		}
 	#write.table(newSurvival, file = "nSur.csv", sep = ",", row.names = FALSE)
}
